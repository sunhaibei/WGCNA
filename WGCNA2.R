library(Seurat)
library(tidyverse)

# 选取常见T细胞耗竭基因
exhaustion_genes <- c("PDCD1", "CTLA4", "LAG3", "TIGIT", "HAVCR2")

# 计算基因表达相关性
cor_matrix <- cor(t(GetAssayData(seurat_cd8, slot = "data")[exhaustion_genes, ]), 
                  t(GetAssayData(seurat_cd8, slot = "data")), 
                  method = "pearson")

# 提取与T细胞耗竭基因显著相关的基因
cor_threshold <- 0.6  # 相关性阈值
p_threshold <- 0.05   # p值阈值
cor_results <- as.data.frame(cor_matrix) %>%
  rownames_to_column("Gene") %>%
  pivot_longer(-Gene, names_to = "Exhaustion_Gene", values_to = "Correlation") %>%
  filter(abs(Correlation) > cor_threshold)

# 显示前10个最相关的基因
head(cor_results, 10)

immune_checkpoints <- c("PDCD1", "CTLA4", "LAG3", "TIGIT", "HAVCR2", "CD274") # 已知免疫检查点
novel_candidates <- rownames(deg_up)[rownames(deg_up) %in% immune_checkpoints]
novel_candidates

power <- 6  # 选择合适的power
net <- blockwiseModules(datExpr, power = power, TOMType = "unsigned", minModuleSize = 30, mergeCutHeight = 0.25)
module_genes <- names(net$colors)[net$colors == "blue"]  # 选一个显著模块
potential_checkpoints <- intersect(module_genes, rownames(deg_up))



# 使用模块的基因表达矩阵
module_genes <- names(net$colors)[net$colors == "blue"]  # 假设选择 "blue" 模块
module_expression <- datExpr[, module_genes]

# 绘制热图
library(pheatmap)
pheatmap(module_expression, cluster_rows = TRUE, cluster_cols = TRUE, 
         scale = "row", show_rownames = FALSE, show_colnames = FALSE)

# 检查是否存在非有限值
sum(!is.finite(module_expression))

str(module_expression)
class(module_expression)


str(module_expression)

# 或者
sapply(module_expression, class)

# 将非数值型列转换为数值型
module_expression <- data.frame(lapply(module_expression, function(x) as.numeric(as.character(x))))
# 检查是否存在非有限值
sum(!is.finite(as.matrix(module_expression)))

# 检查是否存在全 0 的行
zero_rows <- rowSums(module_expression == 0) == ncol(module_expression)
print(sum(zero_rows))  # 输出全 0 的行数

# 检查是否存在全 0 的列
zero_cols <- colSums(module_expression == 0) == nrow(module_expression)
print(sum(zero_cols))  # 输出全 0 的列数

# 检查是否存在常数的行
constant_rows <- apply(module_expression, 1, function(x) length(unique(x)) == 1)
print(sum(constant_rows))  # 输出常数的行数

# 检查是否存在常数的列
constant_cols <- apply(module_expression, 2, function(x) length(unique(x)) == 1)
print(sum(constant_cols))  # 输出常数的列数

# 删除全 0 或常数的行
module_expression <- module_expression[rowSums(module_expression != 0) > 0, ]
module_expression <- module_expression[apply(module_expression, 1, function(x) length(unique(x)) > 1), ]

# 删除全 0 或常数的列
module_expression <- module_expression[, colSums(module_expression != 0) > 0]
module_expression <- module_expression[, apply(module_expression, 2, function(x) length(unique(x)) > 1)]












# 处理非有限值
module_expression[!is.finite(as.matrix(module_expression))] <- 0
library(pheatmap)
pheatmap(module_expression, cluster_rows = TRUE, cluster_cols = TRUE, 
         scale = "row", show_rownames = FALSE, show_colnames = FALSE)

# 筛选方差最大的前 100 个基因
gene_var <- apply(module_expression, 1, var)
top_genes <- names(sort(gene_var, decreasing = TRUE)[1:100])
module_expression_filtered <- module_expression[top_genes, ]

# 绘制热图
pheatmap(module_expression_filtered, cluster_rows = TRUE, cluster_cols = TRUE, 
         scale = "row", show_rownames = TRUE, show_colnames = FALSE, 
         fontsize_row = 8, fontsize_col = 8)

# 保存为高分辨率图片
png("heatmap.png", width = 2000, height = 2000, res = 300)
pheatmap(module_expression_filtered, cluster_rows = TRUE, cluster_cols = TRUE, 
         scale = "row", show_rownames = TRUE, show_colnames = FALSE)
dev.off()
# 对数据进行行标准化
module_expression_scaled <- t(scale(t(module_expression)))

# 计算每行表达量的绝对值最大值
max_abs_value <- apply(module_expression_scaled, 1, function(x) max(abs(x)))

# 筛选绝对值最大值较大的基因（例如 > 2）
filtered_genes <- names(max_abs_value[max_abs_value > 2])
module_expression_filtered <- module_expression[filtered_genes, ]
