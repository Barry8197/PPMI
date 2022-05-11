library(DESeq2)
library(variancePartition)
library(ggplot2)

count_mtx <- NA
for (DS in c('Idiopathic_PD' , 'Genetic_PD' ,'Healthy_Control', 'Genetic_Unaffected')) {
  for (i in list.files(paste0('PPMI_IP/Data_2/',DS,'/'))[-c(3)]) {
    count_mtx_tmp <- read.csv(paste0('~/PPMI_IP/Data_2/',DS,'/',i,'/count_matrix.csv') , check.names = FALSE)
    count_mtx <- cbind(count_mtx_tmp[sample(colnames(count_mtx_tmp[,2:length(count_mtx_tmp)]) , 272)] , count_mtx) 
  }
}
rownames(count_mtx) <- count_mtx_tmp$Geneid

genes_removed <- read.csv('PPMI_IP/Data/genes_removed.csv' , header = FALSE)$V1
genes_subset <- setdiff(rownames(count_mtx) , genes_removed)
count_mtx <- count_mtx[genes_subset , ]
count_mtx <- count_mtx[,1:1632]

meta_df <- read.csv('./PPMI_IP/Data/meta_filtered.csv')

i <- match(meta_df$HudAlphaSampleName, names(count_mtx), nomatch = 0)
count_mtx <- count_mtx[,i]
meta_df <- meta_df[which(meta_df$HudAlphaSampleName %in% names(count_mtx)),]

identical(names(count_mtx) , meta_df$HudAlphaSampleName)

visit <- factor(meta_df$Clinical_Event, levels = c("BL", "V02", "V04", "V06", "V08"))
subject <- meta_df$PATNO
disease <- factor(gsub(" ", "", meta_df$Disease_Status), levels = c("HealthyControl", "IdiopathicPD", "GeneticPD", "GeneticUnaffected"))
PD_status <- factor(disease, levels = c("Unaffected", "Affected"))
PD_status[disease %in% c("HealthyControl", "GeneticUnaffected")] <- "Unaffected"
PD_status[disease %in% c("IdiopathicPD", "GeneticPD")] <- "Affected"
group <- factor(paste(disease, visit, sep = "_"))
age <- factor(meta_df$ageBin, levels = c("under_55", "55_to_65", "over_65"))
sex <- factor(meta_df$Sex,
              levels = c("Female", "Male"))
plate <- factor(meta_df$Plate)
usableBases <- meta_df$Usable_Bases

coldata <- data.frame( sex , age , PD_status , plate , usableBases, row.names = names(count_mtx))

keep <- filterByExpr(count_mtx, group = PD_status)
sum(keep) #18996 genes

dds = DESeqDataSetFromMatrix(countData = count_mtx[keep,], colData = coldata , design = ~0  + sex + plate + age + usableBases + PD_status)

dds = DESeq(dds)


DE_info = results(dds)

de <- as.data.frame(DE_info)
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
  ggtitle('Case-Control Volcano Plot DESeq') + theme(plot.title = element_text(hjust = 0.5))

