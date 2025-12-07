library(tradeSeq)
library(dplyr)

file_path <- "D:\\wsl\\elsa\\data\\aging-hsc-young_kowalczyk.rds"
dyno_object <- readRDS(file_path)

# 1. 原始 counts
counts_raw <- dyno_object$counts          # 493 x 2227 (cells x genes)
cell_ids <- dyno_object$cell_ids
rownames(counts_raw) <- cell_ids

# 2. 基本检查
cat("原始维度(细胞×基因):", dim(counts_raw), "\n")
cat("是否存在 NA:", any(is.na(counts_raw)), "\n")
cat("最小/最大值:", range(counts_raw), "\n")

# 3. 转为整数（如果不是） - 通常 counts 就是整数
if (!all(counts_raw == floor(counts_raw))) {
  counts_raw <- round(counts_raw)
  cat("已将 counts 四舍五入转换为整数。\n")
}

# 4. 去除全零基因 & 全零细胞
gene_totals <- colSums(counts_raw)
cell_totals <- rowSums(counts_raw)

zero_genes <- gene_totals == 0
zero_cells <- cell_totals == 0

if (any(zero_genes)) {
  cat("移除全零基因数:", sum(zero_genes), "\n")
  counts_raw <- counts_raw[, !zero_genes]
}
if (any(zero_cells)) {
  cat("移除全零细胞数:", sum(zero_cells), "\n")
  counts_raw <- counts_raw[!zero_cells, ]
}

cat("过滤后维度:", dim(counts_raw), "\n")

# 5. 简单过滤：只保留在至少 5% 细胞中 count > 0 的基因
presence_pct <- colSums(counts_raw > 0) / nrow(counts_raw) * 100
keep_genes <- presence_pct >= 5
counts_filt <- counts_raw[, keep_genes]
cat("基因存在率过滤后维度:", dim(counts_filt), "\n")

# 6. 减少基因数做最小测试（例如取前 200 个）
if (ncol(counts_filt) > 200) {
  counts_filt <- counts_filt[, 1:200]
  cat("子集化到前200个基因测试。\n")
}

# 7. 提取伪时间并按其排序
pseudotime_df <- dyno_object$progressions %>%
  select(cell_id, percentage) %>%
  distinct(cell_id, .keep_all = TRUE)   # 防止重复
pseudotime_df <- pseudotime_df[order(pseudotime_df$percentage), ] # **关键：显式排序**

pseudotime_vec <- pseudotime_df$percentage
names(pseudotime_vec) <- pseudotime_df$cell_id

# 对齐到过滤后的细胞
common_cells <- intersect(rownames(counts_filt), names(pseudotime_vec))
counts_filt <- counts_filt[common_cells, ]
pseudotime_vec <- pseudotime_vec[common_cells]
# 再次按伪时间排序 counts 和 pseudotime
counts_filt <- counts_filt[names(pseudotime_vec), ] # **再次确保顺序一致**

cat("最终用于拟合的细胞数:", length(pseudotime_vec), "\n")
cat("伪时间范围:", range(pseudotime_vec), "\n")

# 8. 转换为 genes x cells (注意：tradeSeq 期望的是 genes x cells)
counts_for_tradeSeq <- t(counts_filt)   # genes x cells

# 9. cellWeights 向量（单一路径）
cellWeights_vec <- rep(1, length(pseudotime_vec))

# 10. **关键修改：对 counts 矩阵添加微小噪声以提高鲁棒性**
# 这是解决 NULL 问题的常用技巧
cat("对计数矩阵添加微小噪声以提高稳定性...\n")
set.seed(123) # 为了结果可复现
counts_for_tradeSeq <- counts_for_tradeSeq + abs(rnorm(n = length(counts_for_tradeSeq), mean = 0, sd = 1e-6))

# 11. 核心拟合（最简形式）
cat("开始最小 fitGAM 测试...\n")
gam_min <- tryCatch(
  tradeSeq::fitGAM(
    counts = counts_for_tradeSeq,
    pseudotime = pseudotime_vec,      # 向量
    cellWeights = cellWeights_vec,    # 向量
    nknots = 6,
    verbose = TRUE,
    parallel = FALSE
  ),
  error = function(e) e
)

if (inherits(gam_min, "error")) {
  cat("添加噪声后，最小测试仍失败。错误：", gam_min$message, "\n")
} else {
  cat("添加噪声后，最小测试成功。返回对象类型:", class(gam_min), "\n")
  # 统计成功的拟合（非NULL）
  success_gam_list <- gam_min[!sapply(gam_min, is.null)]
  success_ct <- length(success_gam_list)
  cat("成功拟合的基因数:", success_ct, "\n")
  if (success_ct > 0) {
    cat("尝试对第1个成功基因做平滑预测。\n")
    tryCatch({
      test_pred <- tradeSeq::predictSmooth(
        models = success_gam_list[1], # 传入列表中的第一个模型
        nPoints = 50,
        tidy = FALSE
      )
      print(paste("预测矩阵维度 (基因, 伪时间点):", dim(test_pred)))
    }, error = function(e) cat("预测失败：", e$message, "\n"))
  } else {
    cat("所有基因的拟合都返回了 NULL。\n")
  }
}