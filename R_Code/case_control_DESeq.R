library(DESeq2)
library(ggplot2)
library(pheatmap)
library(dplyr)

setwd('~/PPMI_IP/R_Code/')

# Load Normalised Dataset -------------------------------------------------
load('./../Data/preprocessedData/preprocessed_gPD_BL.RData')

vsd <- vst(dds)
res <- results(dds)
summary(res)
hist(res$pvalue)
hist(res$pvalue[!is.na(res$padj)])

# Volcano Plot -------------------
de <- genes_info
de$diffexpressed <- "NO"
# if log2Foldchange > 0.1 and adjusted pvalue < 0.05, set as "UP" 
de$diffexpressed[de$log2FoldChange > 0.1 & de$padj < 0.05] <- "UP"
# if log2Foldchange < -0.1 and adjusted pvalue < 0.05, set as "DOWN"
de$diffexpressed[de$log2FoldChange < -0.1 & de$padj < 0.05] <- "DOWN"
de$delabel <- NA
de$delabel[de$diffexpressed != "NO"] <- rownames(de)[de$diffexpressed != "NO"]

ggplot(data=de, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed , label=delabel)) +
  geom_point() + 
  theme_minimal() +
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=c(-0.1, 0.1), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red") + 
  ggtitle('gPD-HC BL Volcano Plot DESeq') + theme(plot.title = element_text(hjust = 0.5)) + 
  xlim(c(-0.75,0.5))

sum(de$diffexpressed != 'NO')

# Heatmap -----------------------------------------------------------------
library(pheatmap)
stopifnot(rownames(vsd) == rownames(res))
mat <- assay(vsd)
rownames(mat) <- ifelse(!is.na(rowData(vsd)$SYMBOL), rowData(vsd)$SYMBOL, rownames(vsd))
mat <- mat[head(order(res$padj), 30), ]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(vsd)[, c("PD_status"), drop = FALSE])
df <- df %>% arrange(desc(PD_status))
pheatmap(mat[,rownames(df)], annotation_col = df , show_rownames = FALSE , show_colnames = FALSE)


plot(metadata(res)$filterNumRej, 
     type="b", ylab="number of rejections",
     xlab="quantiles of filter")
lines(metadata(res)$lo.fit, col="red")
abline(v=metadata(res)$filterTheta)

write.csv(de , file = './../Data/output_data/case_control_deseq_filt.csv')
