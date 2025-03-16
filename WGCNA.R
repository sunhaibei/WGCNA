library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(limma)
library(DESeq2)
library(WGCNA)

library(Seurat)

# 指定路径
normal_path <- "D:/GSE175453_RAW (1)/GSE533785"
sepsis_path <- "D:/GSE175453_RAW (1)/GSM5333788" # 请替换为正确的脓毒症数据路径

# 读取数据
normal_data <- Read10X(data.dir = normal_path)
sepsis_data <- Read10X(data.dir = sepsis_path)

# 检查数据层
print(names(normal_data))  # 需要是 "Gene Expression"
print(names(sepsis_data))

# 选择 "Gene Expression" 作为表达矩阵
normal_expr <- normal_data$`Gene Expression`
sepsis_expr <- sepsis_data$`Gene Expression`

# 统一基因名格式（替换 "_" 为 "-"）
rownames(normal_expr) <- gsub("_", "-", rownames(normal_expr))
rownames(sepsis_expr) <- gsub("_", "-", rownames(sepsis_expr))

# 确保两个数据集的基因名一致（取交集）
common_genes <- intersect(rownames(normal_expr), rownames(sepsis_expr))
normal_expr <- normal_expr[common_genes, ]
sepsis_expr <- sepsis_expr[common_genes, ]

# 创建Seurat对象
normal_seurat <- CreateSeuratObject(counts = normal_expr, project = "Normal", min.cells = 3, min.features = 200)
sepsis_seurat <- CreateSeuratObject(counts = sepsis_expr, project = "Sepsis", min.cells = 3, min.features = 200)

# 检查对象是否成功创建
print(normal_seurat)
print(sepsis_seurat)

# 合并两个数据集
seurat_combined <- merge(normal_seurat, y = sepsis_seurat, add.cell.ids = c("Normal", "Sepsis"))

# 计算线粒体基因比例
seurat_combined[["percent.mt"]] <- PercentageFeatureSet(seurat_combined, pattern = "^MT-")

# 质量控制：去除低质量细胞
seurat_combined <- subset(seurat_combined, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)



# 数据标准化
seurat_combined <- NormalizeData(seurat_combined)
seurat_combined <- FindVariableFeatures(seurat_combined, selection.method = "vst", nfeatures = 2000)



# 线性缩放
seurat_combined <- ScaleData(seurat_combined)

# PCA降维
seurat_combined <- RunPCA(seurat_combined)

# UMAP/tSNE 聚类
seurat_combined <- FindNeighbors(seurat_combined, dims = 1:20)
seurat_combined <- FindClusters(seurat_combined, resolution = 0.5)
seurat_combined <- RunUMAP(seurat_combined, dims = 1:20)

# 可视化
DimPlot(seurat_combined, reduction = "umap", label = TRUE, group.by = "orig.ident")


# 典型CD8+ T细胞标志基因
FeaturePlot(seurat_combined, features = c("CD8A", "CD8B"))
VlnPlot(seurat_combined, features = c("CD8A", "CD8B"))

# 提取CD8+ T细胞群
Idents(seurat_combined) <- "seurat_clusters"
cd8_cluster <- WhichCells(seurat_combined, expression = CD8A > 1)
seurat_cd8 <- subset(seurat_combined, cells = cd8_cluster)




head(seurat_cd8@meta.data)

seurat_cd8$orig.ident <- factor(seurat_cd8$orig.ident, levels = c("Normal", "Sepsis"))

Idents(seurat_cd8) <- "orig.ident"


levels(seurat_cd8)  # 确保包含 "Normal" 和 "Sepsis"

table(seurat_cd8$orig.ident)  # 统计分组情况


print(Assays(seurat_cd8))  # 查看数据包含哪些层

DefaultAssay(seurat_cd8) <- "RNA"
seurat_cd8 <- JoinLayers(seurat_cd8)
deg_cd8 <- FindMarkers(seurat_cd8, ident.1 = "Sepsis", ident.2 = "Normal",
                       min.pct = 0.25, logfc.threshold = 0.25)


# 找出脓毒症 vs 正常 的差异基因(无效)
deg_cd8 <- FindMarkers(seurat_cd8, ident.1 = "Sepsis", ident.2 = "Normal", min.pct = 0.25, logfc.threshold = 0.25)

# 只保留上调基因
deg_cd8_up <- deg_cd8 %>% filter(avg_log2FC > 1 & p_val_adj < 0.05)



immune_checkpoints <- c("PDCD1", "CTLA4", "LAG3", "TIGIT", "HAVCR2", "CD274") # 经典免疫检查点
novel_candidates <- rownames(deg_cd8_up)[!(rownames(deg_cd8_up) %in% immune_checkpoints)]

# WGCNA
library(WGCNA)
datExpr <- as.data.frame(t(seurat_cd8@assays$RNA@layers$counts))

slotNames(seurat_cd8@assays$RNA)
datExpr <- as.data.frame(t(GetAssayData(seurat_cd8, slot = "scale.data")))
names(seurat_cd8@assays$RNA@layers)
datExpr <- as.data.frame(t(seurat_cd8@assays$RNA@layers$counts))

library(WGCNA)

# 允许并行计算（提高运行速度）
enableWGCNAThreads()

# 确保基因表达数据符合 WGCNA 要求
goodGenes <- goodSamplesGenes(datExpr, verbose = 3)
datExpr <- datExpr[goodGenes$goodSamples, goodGenes$goodGenes]

# 选择最佳 soft-threshold power（差点跑蹦了）
powers <- c(1:20)
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
softPower <- sft$powerEstimate

# 构建 WGCNA 网络
net <- blockwiseModules(datExpr, power = softPower, TOMType = "unsigned",
                        minModuleSize = 30, mergeCutHeight = 0.25)

# 获取关键模块的基因
module_genes <- names(net$colors)[net$colors == "blue"]  # 选择显著模块


library(pheatmap)
top_deg_genes <- rownames(deg_cd8_up)[1:50]  # 取前50个上调基因
pheatmap(datExpr[top_deg_genes, ], cluster_rows = TRUE, cluster_cols = TRUE, 
         scale = "row", show_rownames = TRUE, 
         main = "Top DEG Genes Heatmap")




moduleColors <- net$colors
table(moduleColors)  # 查看模块分布

# 可视化基因模块
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]], 
                    "Module colors", main = "Gene dendrogram and module colors")





library(VennDiagram)
venn.plot <- draw.pairwise.venn(area1 = length(immune_checkpoints), 
                                area2 = length(novel_candidates),
                                cross.area = length(intersect(immune_checkpoints, novel_candidates)), 
                                category = c("Classic ICPs", "Novel Candidates"), 
                                fill = c("red", "blue"))
grid.draw(venn.plot)





power <- 6  # 选择合适的power
net <- blockwiseModules(datExpr, power = power, TOMType = "unsigned", minModuleSize = 30, mergeCutHeight = 0.25)
module_genes <- names(net$colors)[net$colors == "blue"]  # 选显著模块


potential_checkpoints <- intersect(module_genes, rownames(deg_cd8_up))

VlnPlot(seurat_cd8, features = potential_checkpoints, group.by = "orig.ident")
DotPlot(seurat_cd8, features = potential_checkpoints, group.by = "orig.ident")


install.packages("VennDiagram")
library(VennDiagram)
venn.plot <- draw.pairwise.venn(area1 = length(immune_checkpoints), 
                                area2 = length(novel_candidates),
                                cross.area = length(intersect(immune_checkpoints, novel_candidates)), 
                                category = c("Classic ICPs", "Novel Candidates"), 
                                fill = c("red", "blue"))
grid.draw(venn.plot)

install.packages("ggVennDiagram")
library(ggVennDiagram)

sets <- list(
  "Classic ICPs" = immune_checkpoints,
  "Novel Candidates" = novel_candidates
)

ggVennDiagram(sets)





















# 归一化 & 标准化
seurat_combined <- merge(normal_seurat, y = sepsis_seurat, add.cell.ids = c("Normal", "Sepsis"))
seurat_combined <- NormalizeData(seurat_combined)
seurat_combined <- FindVariableFeatures(seurat_combined, selection.method = "vst", nfeatures = 2000)

# 线性缩放
seurat_combined <- ScaleData(seurat_combined)

# PCA降维
seurat_combined <- RunPCA(seurat_combined)

# 进行UMAP/tSNE聚类
seurat_combined <- FindNeighbors(seurat_combined, dims = 1:20)
seurat_combined <- FindClusters(seurat_combined, resolution = 0.5)
seurat_combined <- RunUMAP(seurat_combined, dims = 1:20)


# 使用CD8 T细胞特异性基因进行标注
FeaturePlot(seurat_combined, features = c("CD8A", "CD8B"))

Idents(seurat_combined) <- "seurat_clusters"
cd8_cluster <- WhichCells(seurat_combined, expression = CD8A > 1)
seurat_cd8 <- subset(seurat_combined, cells = cd8_cluster)


# 查看 Idents
table(Idents(seurat_cd8))
head(seurat_cd8@meta.data)
table(seurat_cd8@meta.data$orig.ident)
# 假设 sepsis_seurat 是脓毒症样本的 Seurat 对象
combined_seurat <- merge(seurat_cd8, sepsis_seurat, add.cell.ids = c("Normal", "Sepsis"))

# 查看合并后的对象
combined_seurat
# 合并正常样本和脓毒症样本
combined_seurat <- merge(seurat_cd8, sepsis_seurat, add.cell.ids = c("Normal", "Sepsis"))

# 查看合并后的对象
combined_seurat

# 根据 orig.ident 设置分组信息
combined_seurat@meta.data$group <- ifelse(combined_seurat@meta.data$orig.ident == "Normal", "Normal", "Sepsis")

# 设置 Idents
Idents(combined_seurat) <- combined_seurat@meta.data$group

# 检查 Idents
table(Idents(combined_seurat))

# 运行 FindMarkers
markers <- FindMarkers(combined_seurat, ident.1 = "Sepsis", ident.2 = "Normal", min.pct = 0.25, logfc.threshold = 0.25)

# 查看结果
head(markers)

names(seurat_cd8@assays$RNA@layers)

markers <- FindMarkers(seurat_cd8, ident.1 = "Sepsis", ident.2 = "Normal", min.pct = 0.25, logfc.threshold = 0.25)
deg_up <- markers %>% filter(avg_log2FC > 1 & p_val_adj < 0.05)

immune_checkpoints <- c("PDCD1", "CTLA4", "LAG3", "TIGIT", "HAVCR2", "CD274") # 已知免疫检查点
novel_candidates <- rownames(deg_up)[!(rownames(deg_up) %in% immune_checkpoints)]
datExpr <- as.data.frame(t(seurat_cd8@assays$RNA@data))
power <- 6  # 选择合适的power
net <- blockwiseModules(datExpr, power = power, TOMType = "unsigned", minModuleSize = 30, mergeCutHeight = 0.25)
module_genes <- names(net$colors)[net$colors == "blue"]  # 选一个显著模块
potential_checkpoints <- intersect(module_genes, rownames(deg_up))
