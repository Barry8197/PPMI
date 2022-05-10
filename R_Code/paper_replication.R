library(DESeq2)
library(variancePartition)
library(tidyverse) #42; base for data/table manipulation
library(tidyr) #for pivot tables
library(edgeR) #expression normalization
library(limma) #more expression normalization
library(corrplot) #correlation plot matrix
library(beepr) #beep alerts when job is complete
library(ggpubr) #paper level visuals
library(seqsetvis) #ggplot2 venn

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

design <- model.matrix(~ 0 + PD_status + sex + plate + age + usableBases)
colnames(design) <- gsub("PD_status", "", colnames(design))
colnames(design)

keep <- filterByExpr(count_mtx, group = PD_status)
sum(keep) #18996 genes


dge <- DGEList(count_mtx[keep,])
dge  <- calcNormFactors(dge)

logCPM <- cpm(dge, log = TRUE, prior.count = 3)


#make table of mean expression for each group
PDGroup <- meta_df$HudAlphaSampleName[meta_df$Disease.Status %in% c("Idiopathic PD", "Genetic PD")]
conGroup <- meta_df$HudAlphaSampleName[meta_df$Disease.Status %in% c("Healthy Control", "Genetic Unaffected")]


PDCPM <- logCPM[,PDGroup]
conCPM <- logCPM[,conGroup]

meanCPM <- data.frame("gene_id" = row.names(logCPM),
                      "AvgExpr_allPD" = rowMeans(PDCPM),
                      "AvgExpr_allControls" = rowMeans(conCPM))

v <- voom(dge, design, plot=TRUE) #check voom plot for curve to see if we need to do more filtering
vfit <- lmFit(v, design)

contr.matrix <- makeContrasts(
  PD = Affected - Unaffected,
  levels = colnames(design)
)

vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit)
summary(decideTests(efit))


dt <- decideTests(efit)
contrName <- colnames(contr.matrix)[1]
print(contrName)
df <- rownames_to_column(topTable(efit, coef = 1, n = Inf), var = "gene_id")
with(df, plot(logFC, -log10(P.Value), pch = 20, main = paste("Volcano Plot: ", contrName, sep = "")))
with(subset(df, adj.P.Val < 0.05 & abs(logFC) > 0.1), points(logFC, -log10(P.Value), pch=20, col="red"))
dev.off()


