library(DESeq2)
library(ggplot2)
library(pheatmap)
library(dplyr)

setwd('~/PPMI_IP/R_Code/')

# Load Normalised Dataset -------------------------------------------------
TS <- 'V08'
load(paste0('./../gPD_HC/preprocessedData/',TS,'/preprocessed2.RData'))

genes_of_interest <- c('ENSG00000188906' , 'ENSG00000145335' , 'ENSG00000177628')
genes_info$ID <- rownames(genes_info)
genes_info$ID <- gsub("\\..*","", genes_info$ID)
genes_info[which(genes_info$ID %in% genes_of_interest),]

vsd <- vst(dds)
res <- results(dds)
summary(res)
hist(res$pvalue)
hist(res$pvalue[!is.na(res$padj)])

# Volcano Plot -------------------
de <- genes_info
de$diffexpressed <- "NO"
# if log2Foldchange > 0.1 and adjusted pvalue < 0.05, set as "UP" 
de$diffexpressed[de$log2FoldChange > 0.1 & de$padj < 0.1] <- "UP"
# if log2Foldchange < -0.1 and adjusted pvalue < 0.05, set as "DOWN"
de$diffexpressed[de$log2FoldChange < -0.1 & de$padj < 0.1] <- "DOWN"
de$delabel <- NA
de$delabel[de$diffexpressed != "NO"] <- rownames(de)[de$diffexpressed != "NO"]

ggplot(data=de, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed , label=delabel)) +
  geom_point() + 
  theme_minimal() +
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=c(-0.1, 0.1), col="red") +
  geom_hline(yintercept=-log10(0.1), col="red") + 
  ggtitle(paste0('iPD-HC ',TS,' Volcano Plot DESeq')) + theme(plot.title = element_text(hjust = 0.5))  
  #xlim(c(-0.75,0.5))

sum(de$diffexpressed != 'NO')

# Heatmap -----------------------------------------------------------------
library(pheatmap)
stopifnot(rownames(vsd) == rownames(res))
mat <- datExpr
mat <- mat[head(order(res$padj), 50), ]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(vsd)[, c("PD_status"), drop = FALSE])
pheatmap(mat[,rownames(df)], annotation_col = df , show_colnames = FALSE , show_rownames = FALSE)


DEGs <- res %>% data.frame %>% arrange(padj) %>% mutate(de = (abs(log2FoldChange) > 0.1) & (padj < 0.1)) 

write.csv(DEGs , file = paste0('./../gPD_HC/outputData/DEG_',TS,'.csv'))
write.csv(de , file = './../Data/output_data/case_control_deseq_filt.csv')
