library(DESeq2)
library(variancePartition)
library(ggplot2)

count_mtx_PD <- NA
for (DS in c('Idiopathic_PD' , 'Genetic_PD')) {
  for (i in c('BL' , 'V02' , 'V04' , 'V06' , 'V08')) {
    if (length(colnames(count_mtx_PD)) == 0) {
      count_mtx_PD <- read.csv(paste0('~/PPMI_IP/Data/',DS,'/',i,'/count_matrix.csv') , check.names = FALSE)
      gene_ids = count_mtx_PD$Geneid
      count_mtx_PD = count_mtx_PD[,-1]
    }
    else {
      count_mtx_PD <- cbind(read.csv(paste0('~/PPMI_IP/Data/',DS,'/',i,'/count_matrix.csv') , check.names = FALSE) , count_mtx_PD) 
      count_mtx_PD = count_mtx_PD[,-1]
    }
  }
}

count_mtx_HC <- NA
for (HC in c('Healthy_Control', 'Genetic_Unaffected')) {
  for (i in c('BL' , 'V02' , 'V04' , 'V06' , 'V08')) {
    if (length(colnames(count_mtx_PD)) == 0) {
      count_mtx_HC <- read.csv(paste0('~/PPMI_IP/Data/',HC,'/',i,'/count_matrix.csv') , check.names = FALSE)
      count_mtx_HC = count_mtx_HC[,-1]
    }
    else {
      count_mtx_HC <- cbind(read.csv(paste0('~/PPMI_IP/Data/',HC,'/',i,'/count_matrix.csv') , check.names = FALSE) , count_mtx_HC)
      count_mtx_HC = count_mtx_HC[,-1]
    }
  }
}
count_mtx <- cbind(count_mtx_PD , count_mtx_HC)
count_mtx <- cbind(count_mtx_PD[sample(colnames(count_mtx_PD) , 275)],count_mtx_HC[sample(colnames(count_mtx_HC) , 275)])
rownames(count_mtx) <- gene_ids

genes_removed <- read.csv('PPMI_IP/Data/genes_removed.csv' , header = FALSE)$V1
genes_subset <- setdiff(gene_ids , genes_removed)
count_mtx <- count_mtx[genes_subset , ]

#Subset to genes from paper
paper_df <- read.csv('PPMI_IP/Data/case_control.csv' , skip = 1)
gene_ids_paper <- paper_df$gene_id
count_mtx <- count_mtx[which(count_mtx$Geneid %in% gene_ids_paper) , ]


meta_df <- read.csv('./PPMI_IP/Data/meta_data.11192021.csv')

i <- match(meta_df$HudAlphaSampleName, names(count_mtx), nomatch = 0)
count_mtx <- count_mtx[,i]
meta_df <- meta_df[i,]

identical(names(count_mtx) , meta_df$HudAlphaSampleName)

coldata <- meta_df[which(meta_df$HudAlphaSampleName %in% colnames(count_mtx)) , ]
coldata <- coldata[,c('HudAlphaSampleName','Disease.Status', 'Race' , 'Age..Bin.' , 'Sex' , 'Plate' , 'QC')]
rownames(coldata) <- coldata$HudAlphaSampleName
colnames(coldata)[2] <- 'condition'
coldata$condition <- gsub(" ", "_", coldata$condition, fixed = TRUE)
coldata$condition <- factor(coldata$condition , levels = c("Healthy_Control", "Idiopathic_PD", "Genetic_PD", "Genetic_Unaffected"))
coldata$Sex <- factor(coldata$Sex, levels = c("Female", "Male"))
coldata$age <- factor(coldata$Age..Bin., levels = c("under_55", "55_to_65", "over_65"))
coldata$Plate[coldata$Plate == ""] = 'NA'
coldata$Plate <- factor(coldata$Plate)
coldata$QC <- factor(coldata$QC)
coldata <- coldata[complete.cases(coldata),]
count_mtx <- count_mtx[, rownames(coldata)]

coldata$condition <- gsub("Healthy_Control", "conGroup", coldata$condition, fixed = TRUE)
coldata$condition <- gsub("Genetic_Unaffected", "conGroup", coldata$condition, fixed = TRUE)
coldata$condition <- gsub("Idiopathic_PD", "PDGroup", coldata$condition, fixed = TRUE)
coldata$condition <- gsub("Genetic_PD", "PDGroup", coldata$condition, fixed = TRUE)
coldata$condition <- factor(coldata$condition , levels = c("conGroup","PDGroup"))

keep <- filterByExpr(count_mtx , coldata$condition)
count_mtx <- count_mtx[keep , ]

design <- model.matrix( ~ 0 + coldata$condition + coldata$Sex + coldata$age + coldata$Plate + coldata$QC)

dds = DESeqDataSetFromMatrix(countData = count_mtx, colData = coldata,  design = ~ 0 + condition + Sex + age + Plate + QC)

# To do the differential expression analysis (DEA) with this dataset you only do:
dds = DESeq(dds)

# And to extract the results from the DEA:

# If you want regular log-fold change values
#DE_info = results(dds , lfcThreshold = 0 , altHypothesis = 'greaterAbs')
DE_info = results(dds)



de <- as.data.frame(DE_info)
de$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
de$diffexpressed[de$log2FoldChange > 0.1 & de$padj < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
de$diffexpressed[de$log2FoldChange < -0.1 & de$padj < 0.05] <- "DOWN"
de$delabel <- NA
de$delabel[de$diffexpressed != "NO"] <- rownames(de)[de$diffexpressed != "NO"]

ggplot(data=de, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed , label=delabel)) +
  geom_point() + 
  theme_minimal() +
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=c(-0.1, 0.1), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red") 

de <- as.data.frame(paper_df)
de$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
de$diffexpressed[de$logFC > 0.1 & de$adj.P.Val < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
de$diffexpressed[de$logFC < -0.1 & de$adj.P.Val < 0.05] <- "DOWN"
de$delabel <- NA
de$delabel[de$diffexpressed != "NO"] <- rownames(de)[de$diffexpressed != "NO"]

ggplot(data=de, aes(x=logFC, y=-log10(adj.P.Val), col=diffexpressed , label=delabel)) +
  geom_point() + 
  theme_minimal() +
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=c(-0.1, 0.1), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red") 

diff_expr_paper <- de$gene_id[de$diffexpressed != 'NO']
diff_expr_me <- rownames(de)[de$diffexpressed != 'NO']
intersect(diff_expr_paper,diff_expr_me)
setdiff(diff_expr_paper , rownames(count_mtx))

write.csv(DE_info , file = 'PPMI_IP/Data/iPD_de_info.csv')
