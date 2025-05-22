### Seurat
library(Seurat)
library(Signac)
library(Matrix)

# 1. 加载你已经保存好的整合输入文件
load("D:/HuaweiMoveData/Users/hya/Desktop/CMML/ICA2/Results/scAI_input.Rda")

rna_counts  <- input_list$data$RNA
atac_counts <- input_list$data$ATAC

# 2. 创建 Seurat 对象
seurat_obj <- CreateSeuratObject(counts = rna_counts)
seurat_obj[["ATAC"]] <- CreateAssayObject(counts = atac_counts)

# 3. RNA 模态处理
seurat_obj <- NormalizeData(seurat_obj, assay = "RNA")
seurat_obj <- FindVariableFeatures(seurat_obj, assay = "RNA")
seurat_obj <- ScaleData(seurat_obj, assay = "RNA")
seurat_obj <- RunPCA(seurat_obj, assay = "RNA", reduction.name = "pca")

# 4. ATAC 模态处理
seurat_obj <- RunTFIDF(seurat_obj, assay = "ATAC")
seurat_obj <- FindTopFeatures(seurat_obj, assay = "ATAC")
seurat_obj <- RunSVD(seurat_obj, assay = "ATAC", reduction.name = "lsi")

# 5. 多模态整合（WNN）
seurat_obj <- FindMultiModalNeighbors(
  seurat_obj,
  reduction.list = list("pca", "lsi"),
  dims.list      = list(1:20, 2:20)
)

seurat_obj <- RunUMAP(seurat_obj, nn.name = "weighted.nn", reduction.name = "wnn.umap")
seurat_obj <- FindClusters(seurat_obj, graph.name = "wsnn", resolution = 0.5)

# 6. 提取聚类标签以供后续 ARI/Silhouette 评估
seurat_clusters <- seurat_obj$seurat_clusters

library(ggplot2)
DimPlot(seurat_obj, reduction = "wnn.umap", group.by = "seurat_clusters", label = TRUE) +
  ggtitle("Seurat WNN: UMAP visualization")

embedding_seurat <- Embeddings(seurat_obj, "wnn.umap")  # matrix: cells × 2
labels_seurat <- seurat_obj$seurat_clusters

library(cluster)
library(fpc)

dmat_seurat <- dist(embedding_seurat)
sil_seurat <- silhouette(as.numeric(factor(labels_seurat)), dmat_seurat)
mean_sil_seurat <- mean(sil_seurat[, 3])

cs_seurat <- cluster.stats(as.matrix(dmat_seurat), as.numeric(labels_seurat))
dunn_seurat <- cs_seurat$dunn

cat("Silhouette (Seurat) =", round(mean_sil_seurat, 3), "\n")
cat("Dunn Index (Seurat) =", round(dunn_seurat, 3), "\n")

# 写入 CSV 方便后续外部可视化比对
umap_df <- as.data.frame(embedding_seurat)
umap_df$cluster <- labels_seurat
write.csv(umap_df, file = "seurat_wnn_umap_clusters.csv")


### PCA
# 0. 载入必要包
library(Matrix)
library(ggplot2)
library(cluster)      # silhouette
library(Rtsne)        # optional: for tSNE visualization
library(uwot)         # for UMAP
library(stats)

# 1. 加载你准备好的矩阵数据
load("D:/HuaweiMoveData/Users/hya/Desktop/CMML/ICA2/Results/scAI_input.Rda")

rna_mat  <- as.matrix(input_list$data$RNA)
atac_mat <- as.matrix(input_list$data$ATAC)

# 2.1 过滤掉全为 0 或常数的特征
filter_constant_rows <- function(mat) {
  row_sd <- apply(mat, 1, sd)
  mat[row_sd > 0, , drop = FALSE]
}

rna_mat_filtered  <- filter_constant_rows(rna_mat)
atac_mat_filtered <- filter_constant_rows(atac_mat)

# 2.2 然后再 PCA
rna_pca  <- prcomp(t(rna_mat_filtered), center = TRUE, scale. = TRUE)
atac_pca <- prcomp(t(atac_mat_filtered), center = TRUE, scale. = TRUE)


# 3. 取前 20 主成分，拼接为细胞 × 特征 矩阵
rna_pc  <- rna_pca$x[, 1:20]
atac_pc <- atac_pca$x[, 1:20]
joint_pc <- cbind(rna_pc, atac_pc)

# 4. 在联合空间进行 K-means 聚类（保持与 scAI 一致，设为 K=6）
set.seed(42)
clus <- kmeans(joint_pc, centers = 6)$cluster

# 5. UMAP 可视化
umap_coords <- umap(joint_pc)

df_umap <- data.frame(
  UMAP1 = umap_coords[,1],
  UMAP2 = umap_coords[,2],
  Cluster = factor(clus)
)

ggplot(df_umap, aes(x = UMAP1, y = UMAP2, color = Cluster)) +
  geom_point(size = 0.5) +
  ggtitle("PCA-based Integration: UMAP Visualization") +
  theme_minimal()

# 6. Silhouette 评估聚类质量
dist_mat <- dist(joint_pc)
sil <- silhouette(clus, dist_mat)
mean_sil <- mean(sil[, "sil_width"])
cat("Average Silhouette width =", round(mean_sil, 3), "\n")

# 7. 可选：Dunn index 评估
if (!requireNamespace("fpc", quietly = TRUE)) install.packages("fpc")
library(fpc)
cs <- cluster.stats(d = as.matrix(dist_mat), clustering = clus)
cat("Dunn Index =", round(cs$dunn, 4), "\n")

