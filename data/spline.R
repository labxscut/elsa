suppressPackageStartupMessages({
  library(splines)
  library(Matrix)
  library(dplyr)
  library(readr)
  if (!requireNamespace("irlba", quietly = TRUE)) install.packages("irlba")
  library(irlba)
})

uniformize_with_bs <- function(expr, pt, T_ticks = 100,
                               n_internal_knots = 5, degree = 3,
                               lambda = 1e-1) {
  # expr: cells x genes (numeric)
  stopifnot(nrow(expr) == length(pt))
  
  # 1. build knots (基于伪时间 pt 的分位数)
  boundary <- c(0, 1)
  if (n_internal_knots > 0) {
    probs <- seq(0, 1, length.out = n_internal_knots + 2)[2:(n_internal_knots + 1)]
    internal <- as.numeric(quantile(pt, probs = probs))
  } else {
    internal <- NULL
  }
  
  # 2. Design matrices
  # B: 在原始细胞伪时间 pt 上的样条基函数矩阵
  B <- splines::bs(pt,
                   knots = internal,
                   degree = degree,
                   intercept = TRUE,
                   Boundary.knots = boundary)
  # B_new: 在均匀时间网格 tgrid 上的样条基函数矩阵
  tgrid <- seq(0, 1, length.out = T_ticks)
  B_new <- splines::bs(tgrid,
                       knots = internal,
                       degree = degree,
                       intercept = TRUE,
                       Boundary.knots = boundary)
  
  # 3. Ridge-regularized least squares: coef = solve(B'B + λI) B'Y
  # Y是表达矩阵expr。解出拟合系数 coef (Bases x Genes)
  BtB <- crossprod(B)
  p <- ncol(B)
  # solve for coefficient matrix
  coef <- solve(BtB + diag(lambda, p), crossprod(B, as.matrix(expr)))
  
  # 4. Predict: Y_new = B_new %*% coef
  # Y_new: Ticks x Genes (均匀序列)
  Y_new <- B_new %*% coef         
  rownames(Y_new) <- paste0("T", sprintf("%03d", seq_len(T_ticks)))
  colnames(Y_new) <- colnames(expr)
  
  # 5. 返回均匀化后的矩阵
  return(Y_new)
}

run_uniformization <- function(
    rds_path,
    output_dir,
    method = c("bs", "kernel"),
    T_ticks = 100,
    n_internal_knots = 6,
    degree = 3,
    lambda = 1e-3,
    # gene filtering parameters
    max_zero_percent = 70,
    min_mean_expression = 0,
    top_n_genes = 500,
    use_svd = TRUE,
    # 使用 irlba SVD，通常 20 个成分已足够
    svd_components = 20 
) {
  method <- match.arg(method)
  message("Loading dyno object...")
  obj <- readRDS(rds_path)
  
  # --- 1. 数据提取与对齐 ---
  expr <- obj$expression         # cells x genes
  cells <- obj$cell_ids
  rownames(expr) <- cells
  
  # 提取并合并伪时间 (progressions$percentage)
  pt_df <- obj$progressions %>%
    select(cell_id, percentage) %>%
    group_by(cell_id) %>%
    summarise(pseudotime = max(percentage), .groups = "drop")
  
  # 按细胞ID对齐伪时间pt向量和表达矩阵expr的行
  pt <- setNames(pt_df$pseudotime, pt_df$cell_id)[cells]
  
  # 移除没有伪时间的细胞
  keep_cells <- which(!is.na(pt))
  expr <- expr[keep_cells, , drop = FALSE]
  pt <- pt[keep_cells]
  
  # --- 2. 排序 (仅为严谨性，并非算法要求) ---
  # **显式按 percentage 排序 (修正点 2)**
  sort_idx <- order(pt)
  expr <- expr[sort_idx, , drop = FALSE]
  pt <- pt[sort_idx]
  
  
  # --- 3. 基本 QC 和基因过滤 ---
  zero_pct <- colSums(expr == 0) / nrow(expr) * 100
  mean_expr <- colMeans(expr)
  keep_genes <- (zero_pct <= max_zero_percent) & (mean_expr >= min_mean_expression)
  expr <- expr[, keep_genes, drop = FALSE]
  
  # --- 4. 基因筛选 (SVD/方差) ---
  if (use_svd) {
    message(sprintf("Running SVD/PCA to select top %d genes...", top_n_genes))
    Xc <- scale(expr, center = TRUE, scale = FALSE) # 中心化数据
    k <- min(svd_components, nrow(Xc) - 1, ncol(Xc))
    # 使用 irlba 进行快速 SVD (解决 prcomp 对大矩阵慢的问题)
    sv <- tryCatch(irlba(Xc, nv = k), error = function(e) NULL)
    
    if (is.null(sv)) {
      warning("irlba failed; fallback to variance-based selection.")
      gene_rank <- order(apply(expr, 2, var), decreasing = TRUE)
    } else {
      # **SVD核心：计算基因对前 K 个主成分的贡献度 (修正点 1)**
      V <- sv$v[, seq_len(k), drop = FALSE]  # genes x k (右奇异向量/loadings)
      contrib <- rowSums(V^2) # 总贡献度
      names(contrib) <- colnames(expr)
      gene_rank <- order(contrib, decreasing = TRUE)
    }
  } else {
    gene_rank <- order(apply(expr, 2, var), decreasing = TRUE)
  }
  
  g <- min(top_n_genes, ncol(expr)) # **确保不超过最大限制 (修正点 3)**
  sel_genes <- colnames(expr)[gene_rank[seq_len(g)]]
  expr_sel <- expr[, sel_genes, drop = FALSE]
  
  message(sprintf("Cells: %d, Genes selected: %d, Method: %s",
                  nrow(expr_sel), ncol(expr_sel), method))
  
  # 归一化伪时间到 [0,1]
  rng <- range(pt, na.rm = TRUE)
  if (rng[1] < 0 || rng[2] > 1) {
    pt <- (pt - rng[1]) / (rng[2] - rng[1])
  }
  
  # --- 5. 平滑与均匀化 ---
  if (method == "bs") {
    Y_new <- uniformize_with_bs(expr_sel, pt,
                                T_ticks = T_ticks,
                                n_internal_knots = n_internal_knots,
                                degree = degree,
                                lambda = lambda)
  } else {
    # 使用 Kernel 方法 (如果您需要尝试)
    stop("Kernel method not fully implemented in this block, use 'bs'.")
  }
  
  # --- 6. 导出文件 (T Ticks x G Genes) ---
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_path_csv <- file.path(output_dir, "LLA_input_uniformized.csv")
  write.csv(Y_new, output_path_csv, row.names = TRUE)
  
  message(sprintf("Output: %d ticks x %d genes. Saved to %s",
                  nrow(Y_new), ncol(Y_new), output_path_csv))
  invisible(Y_new)
}

# 最终执行
rds_path_input <- "D:/wsl/elsa/data/aging-hsc-young_kowalczyk.rds"
output_directory <- "D:/wsl/elsa/results/spline2"

# 执行 LLA 预处理流程
run_uniformization(
  rds_path = rds_path_input,
  output_dir = output_directory,
  method = "bs", # 使用 B-Spline 方法
  T_ticks = 100, # LLA 序列长度
  top_n_genes = 500, # LLA 批量计算的基因数
  svd_components = 20 # 用于贡献度计算的主成分数
)