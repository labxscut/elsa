library(dplyr)
library(ggplot2)
library(tibble)

# 1. 提取具有有效伪时间的细胞ID和它们的伪时间值
#    为了方便，将其转换为一个 tibble
cell_pseudotime_map <- pseudotime_valid %>%
  rownames_to_column("cell") %>%
  select(cell, pseudotime)

# 2. 提取这些细胞的原始基因表达矩阵
#    假设你的 cds 对象中，表达矩阵存储在 counts(cds) 中
#    并且，你只想分析高可变基因 (HVG) 以减少计算量
#    (如果你想分析所有基因，可以跳过 HVG 筛选)

# 前面已经通过 preprocess_cds 计算了 HVG，并且num_dim=50
# 如果没有，运行:
# cds <- preprocess_cds(cds, num_dim = 50)

# 提取 HVG 列表
hvg_genes <- fData(cds)$gene_short_name[fData(cds)$is_hvg]

# 提取表达矩阵，并只保留有效细胞和 HVG
# 注意：counts(cds) 是一个稀疏矩阵，行是基因，列是细胞
expression_matrix <- counts(cds)[hvg_genes, cell_pseudotime_map$cell_id]

# 将稀疏矩阵转换为普通矩阵（如果数据量不大）
expression_matrix <- as.matrix(expression_matrix)