

library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(GenomicRanges)
library(dplyr)

# 1. 直接读取你的 markers_ATAC_top.rds
markers_top <- readRDS("markers_ATAC_top.rds")

# 2. 构造 GRanges（峰坐标形式："chr:start-end"）
peak_strings <- markers_top$features
gr_peaks <- GRanges(
  seqnames = sub("^([^:]+):.*$", "\\1", peak_strings),
  ranges   = IRanges(
    start = as.integer(sub("^.*:(\\d+)-.*$", "\\1", peak_strings)),
    end   = as.integer(sub("^.*-(\\d+)$",        "\\1", peak_strings))
  )
)

# 3. 用 annotatePeak() 注释到基因
peak_anno <- annotatePeak(
  gr_peaks,
  TxDb      = TxDb.Hsapiens.UCSC.hg38.knownGene,
  tssRegion = c(-2000, 500),
  annoDb    = "org.Hs.eg.db"
)

# 4. 提取 SYMBOL 并合并
anno_df <- as.data.frame(peak_anno)[, c("seqnames","start","end","SYMBOL")]
anno_df$feature <- paste0(anno_df$seqnames, ":", anno_df$start, "-", anno_df$end)

markers_annotated <- markers_top %>%
  ungroup() %>%
  left_join(anno_df, by = c("features" = "feature"))

# 5. 查看带有 SYMBOL 注释的结果
head(markers_annotated)

# 6. 保存带基因注释的 ATAC 峰表
saveRDS(
  markers_annotated,
  file = "markers_ATAC_annotated.rds"
)
# 也可以导出为 csv 方便 Excel 查看
write.csv(
  markers_annotated,
  file = "markers_ATAC_annotated.csv",
  row.names = FALSE
)
