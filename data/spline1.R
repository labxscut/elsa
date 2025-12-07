suppressPackageStartupMessages({
  library(splines)
  library(dplyr)
  library(irlba)    # SVD
  library(ggplot2)  # 可选，用于可视化
})

# ========== 辅助函数 ==========

#' 从 dyno progressions 和 milestone_network 计算全局伪时间
#' 
#' @param progressions_df dyno 的 progressions 数据框
#' @param milestone_network dyno 的 milestone_network
#' @return data.frame(cell_id, global_pseudotime)
compute_global_pseudotime <- function(progressions_df, milestone_network) {
  stopifnot(all(c("cell_id", "from", "to", "percentage") %in% colnames(progressions_df)))
  stopifnot(all(c("from", "to", "length") %in% colnames(milestone_network)))
  
  # 定义正确的发育顺序（根据您的数据）
  correct_order <- c("LT-HSC", "ST-HSC", "MPP")
  
  # 按发育顺序排列 milestone_network
  milestone_network <- milestone_network %>%
    arrange(match(from, correct_order))
  
  # 计算每段路径的累积起始伪时间
  milestone_network <- milestone_network %>%
    mutate(
      segment_length = length,
      pseudotime_start = cumsum(c(0, segment_length[-n()]))
    )
  
  # 合并 progressions 与路径信息，计算全局伪时间
  global_pt <- progressions_df %>%
    left_join(milestone_network, by = c("from", "to")) %>%
    mutate(
      global_pseudotime = pseudotime_start + percentage * segment_length
    ) %>%
    select(cell_id, global_pseudotime)
  
  # 处理无效值
  global_pt$global_pseudotime[!is.finite(global_pt$global_pseudotime)] <- NA_real_
  
  return(global_pt)
}


#' B-spline 自适应均匀化（核心函数）
#' 
#' @param expr 表达矩阵（细胞 × 基因）
#' @param pt 伪时间向量（与 expr 行对应）
#' @param T_ticks 输出的均匀时间点数
#' @param degree 样条次数（默认3=cubic）
#' @param lambda 岭正则化参数
#' @return 均匀化矩阵（T_ticks × 基因），带 metadata 属性
uniformize_with_bs_adaptive <- function(expr, pt, T_ticks = 100,
                                        degree = 3, lambda = 1e-2) {
  stopifnot(nrow(expr) == length(pt))
  
  unique_pt <- sort(unique(pt))
  n_unique <- length(unique_pt)
  n_cells <- length(pt)
  
  message(sprintf("  Data: %d cells at %d unique pseudotime points", 
                  n_cells, n_unique))
  
  # --- 1. 标准化伪时间到 [0,1] ---
  pt_min <- min(pt)
  pt_max <- max(pt)
  
  if (pt_max - pt_min <= 0) {
    stop("All pseudotime values are identical; cannot uniformize.")
  }
  
  pt_norm <- (pt - pt_min) / (pt_max - pt_min)
  unique_pt_norm <- (unique_pt - pt_min) / (pt_max - pt_min)
  
  message(sprintf("  Normalized pseudotime to [0, 1] (original: [%.3f, %.3f])", 
                  pt_min, pt_max))
  
  # --- 2. 自适应选择内部节点 ---
  if (n_unique <= 2) {
    # 只有1-2个时间点：无内部节点，全局样条
    internal <- NULL
    message("  Strategy: 0 internal knots (≤2 distinct time points)")
    
  } else if (n_unique == 3) {
    # 3个时间点：1个节点在中位数
    internal <- median(unique_pt_norm)
    message(sprintf("  Strategy: 1 internal knot at %.3f (median of 3 points)", 
                    internal))
    
  } else if (n_unique <= 6) {
    # 4-6个时间点：段间中点策略
    internal <- c()
    for (i in 1:(n_unique - 1)) {
      midpoint <- (unique_pt_norm[i] + unique_pt_norm[i+1]) / 2
      internal <- c(internal, midpoint)
    }
    # 移除边界点（确保在开区间内）
    internal <- internal[internal > 0 & internal < 1]
    internal <- sort(unique(internal))
    message(sprintf("  Strategy: %d internal knots (segment midpoints)", 
                    length(internal)))
    
  } else {
    # 多于6个时间点：基于分位数
    n_knots <- min(ceiling(n_unique / 3), 10)
    probs <- seq(0, 1, length.out = n_knots + 2)[2:(n_knots + 1)]
    internal <- as.numeric(quantile(unique_pt_norm, probs = probs))
    internal <- unique(internal)
    internal <- internal[internal > 0 & internal < 1]
    message(sprintf("  Strategy: %d internal knots (quantile-based)", 
                    length(internal)))
  }
  
  # 确保内部节点在开区间 (0, 1) 内
  if (!is.null(internal)) {
    internal <- internal[internal > 0 & internal < 1]
    internal <- sort(unique(internal))
  }
  
  # --- 3. 构建 B-spline 基矩阵 ---
  boundary <- c(0, 1)
  
  B <- splines::bs(pt_norm,
                   knots = internal,
                   degree = degree,
                   intercept = TRUE,
                   Boundary.knots = boundary)
  
  p <- ncol(B)  # 基函数数量
  n <- nrow(B)  # 样本数
  
  message(sprintf("  Basis functions: %d (degree=%d, %d internal knots)", 
                  p, degree, length(internal)))
  
  # --- 4. 检查并处理退化情况 ---
  if (p >= n - 1) {
    warning(sprintf("Too many basis functions (p=%d) for n=%d samples; reducing knots", 
                    p, n))
    # 减少内部节点
    if (length(internal) > 1) {
      internal <- internal[seq(1, length(internal), length.out = ceiling(length(internal)/2))]
    } else {
      internal <- NULL
    }
    # 重新构建
    B <- splines::bs(pt_norm, knots = internal, degree = degree, 
                     intercept = TRUE, Boundary.knots = boundary)
    p <- ncol(B)
    message(sprintf("  Reduced to %d basis functions", p))
  }
  
  # --- 5. 构建预测网格 ---
  tgrid_norm <- seq(0, 1, length.out = T_ticks)
  B_new <- splines::bs(tgrid_norm,
                       knots = internal,
                       degree = degree,
                       intercept = TRUE,
                       Boundary.knots = boundary)
  
  # --- 6. 岭回归拟合系数 ---
  BtB <- crossprod(B)
  ridge_mat <- BtB + diag(lambda, p)
  
  # 添加条件数检查
  cond_num <- kappa(ridge_mat, exact = FALSE)
  if (cond_num > 1e12) {
    warning(sprintf("Matrix condition number is high: %.2e", cond_num))
  }
  
  coef <- tryCatch({
    solve(ridge_mat, crossprod(B, as.matrix(expr)))
  }, error = function(e) {
    warning("solve() failed; using Moore-Penrose pseudoinverse (MASS::ginv)")
    if (!requireNamespace("MASS", quietly = TRUE)) {
      stop("Need MASS package. Install with: install.packages('MASS')")
    }
    MASS::ginv(ridge_mat) %*% crossprod(B, as.matrix(expr))
  })
  
  # --- 7. 预测均匀时间点上的表达量 ---
  Y_new <- B_new %*% coef
  
  rownames(Y_new) <- paste0("T", sprintf("%03d", seq_len(T_ticks)))
  colnames(Y_new) <- colnames(expr)
  
  # --- 8. 保存元数据 ---
  attr(Y_new, "metadata") <- list(
    internal_knots_normalized = internal,
    internal_knots_original = if (!is.null(internal)) {
      internal * (pt_max - pt_min) + pt_min
    } else {
      NULL
    },
    n_basis_functions = p,
    degree = degree,
    lambda = lambda,
    pseudotime_range_original = c(pt_min, pt_max),
    n_unique_timepoints = n_unique
  )
  
  return(Y_new)
}


#' 主均匀化流程
#' 
#' @param rds_path dyno RDS 文件路径
#' @param output_dir 输出目录
#' @param method 方法（目前只支持 "bs"）
#' @param T_ticks 输出时间点数
#' @param degree 样条次数
#' @param lambda 正则化参数
#' @param max_zero_percent 基因过滤：最大零值比例
#' @param min_mean_expression 基因过滤：最小平均表达量
#' @param top_n_genes 筛选的基因数
#' @param use_svd 是否使用 SVD 筛选基因
#' @param svd_components SVD 保留的主成分数
#' @return 均匀化矩阵（隐式返回）
run_uniformization <- function(
    rds_path,
    output_dir,
    method = c("bs"),
    T_ticks = 100,
    degree = 3,
    lambda = 1e-2,
    max_zero_percent = 70,
    min_mean_expression = 0,
    top_n_genes = 500,
    use_svd = TRUE,
    svd_components = 20
) {
  method <- match.arg(method)
  
  message("========================================")
  message("伪时间均匀化流程")
  message("========================================\n")
  
  # --- 1. 加载数据 ---
  message("步骤 1/6: 加载 dyno 对象...")
  obj <- readRDS(rds_path)
  
  expr <- obj$expression
  cells <- obj$cell_ids
  rownames(expr) <- cells
  
  message(sprintf("  原始数据: %d 细胞 × %d 基因", nrow(expr), ncol(expr)))
  
  # --- 2. 计算全局伪时间 ---
  message("\n步骤 2/6: 计算全局伪时间...")
  pt_df <- compute_global_pseudotime(obj$progressions, obj$milestone_network)
  pt <- setNames(pt_df$global_pseudotime, pt_df$cell_id)[cells]
  
  # 移除无伪时间的细胞
  keep_cells <- which(!is.na(pt))
  expr <- expr[keep_cells, , drop = FALSE]
  pt <- pt[keep_cells]
  
  # 按伪时间排序
  sort_idx <- order(pt)
  expr <- expr[sort_idx, , drop = FALSE]
  pt <- pt[sort_idx]
  
  message(sprintf("  伪时间范围: [%.3f, %.3f]", min(pt), max(pt)))
  message(sprintf("  唯一时间点: %s", 
                  paste(sort(unique(pt)), collapse = ", ")))
  message(sprintf("  细胞数: %d", length(pt)))
  
  # --- 3. 基因质量过滤 ---
  message("\n步骤 3/6: 基因质量控制...")
  
  zero_pct <- colSums(expr == 0) / nrow(expr) * 100
  mean_expr <- colMeans(expr)
  
  keep_genes <- (zero_pct <= max_zero_percent) & 
    (mean_expr >= min_mean_expression)
  
  n_removed <- sum(!keep_genes)
  expr <- expr[, keep_genes, drop = FALSE]
  
  message(sprintf("  移除 %d 个低质量基因", n_removed))
  message(sprintf("  保留 %d 个基因", ncol(expr)))
  
  # --- 4. 基因选择（SVD 或方差） ---
  message("\n步骤 4/6: 基因筛选...")
  
  if (use_svd && ncol(expr) > top_n_genes) {
    message("  使用 SVD 方法...")
    
    Xc <- scale(expr, center = TRUE, scale = FALSE)
    k <- min(svd_components, nrow(Xc) - 1, ncol(Xc))
    
    sv <- tryCatch({
      irlba::irlba(Xc, nv = k)
    }, error = function(e) {
      warning("SVD 失败: ", e$message)
      NULL
    })
    
    if (!is.null(sv)) {
      # 计算每个基因对前k个主成分的贡献度
      V <- sv$v[, seq_len(k), drop = FALSE]
      contrib <- rowSums(V^2)
      names(contrib) <- colnames(expr)
      gene_rank <- order(contrib, decreasing = TRUE)
      
      var_explained <- sum(sv$d[1:k]^2) / sum(Xc^2) * 100
      message(sprintf("  SVD 前 %d 个主成分解释 %.1f%% 方差", k, var_explained))
    } else {
      message("  回退到基于方差的排序...")
      gene_rank <- order(apply(expr, 2, var), decreasing = TRUE)
    }
  } else {
    message("  使用方差排序...")
    gene_rank <- order(apply(expr, 2, var), decreasing = TRUE)
  }
  
  g <- min(top_n_genes, ncol(expr))
  sel_genes <- colnames(expr)[gene_rank[seq_len(g)]]
  expr_sel <- expr[, sel_genes, drop = FALSE]
  
  message(sprintf("  最终选择 %d 个基因", g))
  
  # --- 5. B-spline 平滑与均匀化 ---
  message("\n步骤 5/6: B-spline 平滑与均匀化...")
  
  if (method == "bs") {
    Y_new <- uniformize_with_bs_adaptive(
      expr_sel, pt,
      T_ticks = T_ticks,
      degree = degree,
      lambda = lambda
    )
  } else {
    stop("当前仅支持 'bs' 方法")
  }
  
  # --- 6. 保存结果 ---
  message("\n步骤 6/6: 保存结果...")
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # 保存均匀化矩阵
  saveRDS(Y_new, file.path(output_dir, "LLA_input_uniformized.rds"))
  write.csv(Y_new, file.path(output_dir, "LLA_input_uniformized.csv"), 
            row.names = TRUE)
  
  # 保存元数据
  meta <- attr(Y_new, "metadata")
  metadata <- list(
    timestamp = Sys.time(),
    user = Sys.info()["user"],
    input_file = rds_path,
    n_cells = nrow(expr_sel),
    n_genes_selected = ncol(expr_sel),
    T_ticks = T_ticks,
    method = method,
    spline_degree = degree,
    lambda = lambda,
    internal_knots_original = meta$internal_knots_original,
    n_basis_functions = meta$n_basis_functions,
    pseudotime_range = meta$pseudotime_range_original,
    n_unique_timepoints = meta$n_unique_timepoints,
    selected_genes = sel_genes,
    svd_components = if (use_svd) svd_components else NULL
  )
  
  saveRDS(metadata, file.path(output_dir, "metadata.rds"))
  write.csv(as.data.frame(lapply(metadata, function(x) {
    if (is.null(x)) "NULL" else paste(x, collapse = ", ")
  })), file.path(output_dir, "metadata.csv"), row.names = FALSE)
  
  # 保存伪时间信息
  pt_info <- data.frame(
    cell_id = rownames(expr_sel),
    pseudotime_original = pt,
    stringsAsFactors = FALSE
  )
  write.csv(pt_info, file.path(output_dir, "cell_pseudotime.csv"), 
            row.names = FALSE)
  
  message(sprintf("\n✓ 均匀化完成！"))
  message(sprintf("  输出矩阵: %d 时间点 × %d 基因", nrow(Y_new), ncol(Y_new)))
  message(sprintf("  保存路径: %s", output_dir))
  message("========================================\n")
  
  invisible(Y_new)
}

# ========== 使用示例 ==========

# 运行主流程
if (interactive()) {
  Y_result <- run_uniformization(
    rds_path = "D:/wsl/elsa/data/aging-hsc-young_kowalczyk.rds",
    output_dir = "D:/wsl/elsa/results/final",
    method = "bs",
    T_ticks = 100,
    degree = 3,
    lambda = 1e-2,
    max_zero_percent = 70,
    min_mean_expression = 0,
    top_n_genes = 500,
    use_svd = TRUE,
    svd_components = 20
  )
}