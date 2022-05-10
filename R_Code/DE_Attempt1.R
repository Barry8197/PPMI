library(dplyr)
library(DESeq2)

gene_ids_paper <- read.csv('PPMI_IP/Data/gene_ids_prodromal_HC.csv')[,-1]

count_mtx_PD <- read.csv('~/PPMI_IP/Data/Idiopathic_PD/BL/count_matrix.csv' , check.names = FALSE)
count_mtx_HC <- cbind(read.csv('~/PPMI_IP/Data/Healthy_Control/BL/count_matrix.csv' , check.names = FALSE) , read.csv('~/PPMI_IP/Data/Genetic_Unaffected/BL/count_matrix.csv' , check.names = FALSE))
count_mtx <- cbind(count_mtx_PD[sample(colnames(count_mtx_PD) , 212)],count_mtx_HC[sample(colnames(count_mtx_HC) , 212)])
count_mtx <- count_mtx[which(count_mtx$Geneid %in% gene_ids_paper) , ]

# Convert the dataframe to matrix, saving the gene IDs in another vector
gene_ids = count_mtx$Geneid
count_mtx = count_mtx[,-1] 

meta_df <- read.csv('./PPMI_IP/Data/meta_data.11192021.csv')
coldata <- meta_df[which(meta_df$HudAlphaSampleName %in% colnames(count_mtx)) , ]
coldata <- coldata[,c('HudAlphaSampleName','Disease.Status', 'Race')]
rownames(coldata) <- coldata$HudAlphaSampleName
colnames(coldata)[2] <- 'condition'
coldata$condition <- gsub(" ", "_", coldata$condition, fixed = TRUE)
coldata$condition <- gsub("Genetic_Unaffected", "Healthy_Control", coldata$condition, fixed = TRUE)
coldata$condition <- factor(coldata$condition , levels = c('Idiopathic_PD' , 'Healthy_Control'))
count_mtx <- count_mtx[, rownames(coldata)]

dds = DESeqDataSetFromMatrix(countData = count_mtx, colData = coldata,  design = ~ condition)
rownames(dds) <- gene_ids


# To do the differential expression analysis (DEA) with this dataset you only do:
dds = DESeq(dds)

# And to extract the results from the DEA:

# If you want regular log-fold change values
DE_info = results(dds, lfcThreshold = 0, altHypothesis = 'greaterAbs')

library(ggplot2)

de <- as.data.frame(DE_info)
de$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
de$diffexpressed[de$log2FoldChange > 0.6 & de$pvalue < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
de$diffexpressed[de$log2FoldChange < -0.6 & de$pvalue < 0.05] <- "DOWN"
de$delabel <- NA
de$delabel[de$diffexpressed != "NO"] <- rownames(de)[de$diffexpressed != "NO"]

ggplot(data=de, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed , label=delabel)) +
  geom_point() + 
  theme_minimal() +
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red") + 
  xlim(-1.5,1) + ylim(-0.1,8)

library(EnhancedVolcano)
EnhancedVolcano(DE_info,
                lab = rownames(DE_info),
                x = 'log2FoldChange',
                y = 'pvalue')
