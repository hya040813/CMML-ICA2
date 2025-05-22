# 04_prepare_scAI_input.R

# 1) 加载 MatrixMarket 格式的稀疏矩阵
library(Matrix)

rna_counts <- readMM("D:/HuaweiMoveData/Users/hya/Desktop/CMML/ICA2/Results/scAI_rna_counts.mtx")
atac_counts <- readMM("D:/HuaweiMoveData/Users/hya/Desktop/CMML/ICA2/Results/scAI_atac_counts.mtx")

# 2) 读取行／列标签
rna_genes  <- read.delim("D:/HuaweiMoveData/Users/hya/Desktop/CMML/ICA2/Results/scAI_rna_genes.tsv",  header=TRUE, stringsAsFactors=FALSE)$gene
rna_cells  <- read.delim("D:/HuaweiMoveData/Users/hya/Desktop/CMML/ICA2/Results/scAI_rna_barcodes.tsv", header=TRUE, stringsAsFactors=FALSE)$cell
atac_peaks <- read.delim("D:/HuaweiMoveData/Users/hya/Desktop/CMML/ICA2/Results/scAI_atac_peaks.tsv", header=TRUE, stringsAsFactors=FALSE)$peak
atac_cells <- read.delim("D:/HuaweiMoveData/Users/hya/Desktop/CMML/ICA2/Results/scAI_atac_barcodes.tsv", header=TRUE, stringsAsFactors=FALSE)$cell

# 3) 给矩阵打行列名
rownames(rna_counts) <- rna_genes
colnames(rna_counts) <- rna_cells
rownames(atac_counts) <- atac_peaks
colnames(atac_counts) <- atac_cells

# 4) （可选）加载 labels
if (file.exists("D:/HuaweiMoveData/Users/hya/Desktop/CMML/ICA2/Results/scAI_cell_labels.tsv")) {
  labels <- read.delim("Results/scAI_cell_labels.tsv", header=TRUE, stringsAsFactors=FALSE)
} else {
  labels <- NULL
}

# 5) 对齐细胞：只取两端共有的
common_cells <- intersect(colnames(rna_counts), colnames(atac_counts))
rna_counts <- rna_counts[, common_cells]
atac_counts <- atac_counts[, common_cells]
if (!is.null(labels)) {
  labels <- labels[labels$cell %in% common_cells, , drop=FALSE]
}

# 6) 打包进一个 list
input_list <- list(
  data   = list(RNA = rna_counts, ATAC = atac_counts),
  labels = labels
)

# 7) 保存为 Rda 文件
save(input_list, file = "D:/HuaweiMoveData/Users/hya/Desktop/CMML/ICA2/Results/scAI_input.Rda")

message("✔ scAI input written to Results/scAI_input.Rda")
