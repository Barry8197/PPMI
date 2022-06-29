library(edgeR) #expression normalization
library(devtools)
install_github("vqv/ggbiplot")
library(ggbiplot)
library(DESeq2)

# Pull in Count Matrices --------------------------------------------------
count_mtx <- NA
for (DS in c('Idiopathic_PD' , 'Genetic_PD' ,'Healthy_Control', 'Genetic_Unaffected','Prodromal' , 'SWEDD')) {
  for (i in list.files(paste0('PPMI_IP/Data_2/',DS,'/'))) {
    print(DS)
    count_mtx_tmp <- read.csv(paste0('~/PPMI_IP/Data_2/',DS,'/',i,'/count_matrix.csv') , check.names = FALSE)
    count_mtx <- cbind(count_mtx_tmp , count_mtx) 
  }
}
rownames(count_mtx) <- count_mtx_tmp$Geneid #add gene ids to rownames

genes_removed <- read.csv('PPMI_IP/meta_files/genes_removed.csv' , header = FALSE)$V1 #remove genes as per paper
genes_subset <- setdiff(rownames(count_mtx) , genes_removed)
count_mtx <- count_mtx[genes_subset , ]
count_mtx <- count_mtx[,1:length(colnames(count_mtx))-1]


# Create design matrix from meta data -------------------------------------
meta_df <- read.csv('./PPMI_IP/meta_files/metadata_final.csv')

i <- match(meta_df$HudAlphaSampleName, names(count_mtx), nomatch = 0)
count_mtx <- count_mtx[,i]
meta_df <- meta_df[which(meta_df$HudAlphaSampleName %in% names(count_mtx)),]

identical(names(count_mtx) , meta_df$HudAlphaSampleName) #double check ids and genes match

visit <- factor(meta_df$Clinical_Event, levels = c("BL", "V02", "V04", "V06", "V08"))
subject <- meta_df$PATNO
disease <- factor(gsub(" ", "", meta_df$Disease_Status), levels = c("HealthyControl", "IdiopathicPD", "GeneticPD", "GeneticUnaffected" , "Prodromal" , "SWEDD"))
PD_status <- factor(disease, levels = c("Unaffected", "Affected"))
PD_status[disease %in% c("HealthyControl", "GeneticUnaffected")] <- "Unaffected"
PD_status[disease %in% c("IdiopathicPD", "GeneticPD" , "Prodromal" , "SWEDD")] <- "Affected"
group <- factor(paste(disease, visit, sep = "_"))
age <- factor(meta_df$ageBin, levels = c("under_55", "55_to_65", "over_65"))
sex <- factor(meta_df$Sex,
              levels = c("Female", "Male"))
plate <- factor(meta_df$Plate)
usableBases <- meta_df$Usable_Bases

coldata <- data.frame( meta_df$Neutrophils, meta_df$Lymphocytes, PD_status , meta_df$Genetic_Status , plate ,
                       visit, meta_df$Position , meta_df$PD_medication , meta_df$QC ,
                       sex , meta_df$Race , usableBases, age , row.names = names(count_mtx))

# Filter out non-expressed genes ------------------------------------------
keep <- filterByExpr(count_mtx)
sum(keep) #18996 genes


# Obtain R Log  ------------------------------------
dds = DESeqDataSetFromMatrix(countData = count_mtx[keep_tst,] ,  colData = coldata , design = ~0 + sex + plate + age + usableBases + PD_status)
rld <- vst(dds , nsub = 2)

(pca_data <- plotPCA(rld , intgroup = c('sex','plate','age','usableBases','PD_status') , returnData = TRUE))

plot(pca_data$PC1 , pca_data$PC2)

ggplot(pca_data) + geom_point(aes(x = PC1 , y = PC2) , colour = as.numeric(sex))

geneIDs <- c("ENSG00000229807.11",
             "ENSG00000129824.15",
             "ENSG00000280969.1",
             "ENSG00000012817.15",
             "ENSG00000067048.16",
             "ENSG00000114374.12")

tst <- count_mtx[which(rownames(count_mtx) %in% geneIDs) , ]
dds = DESeqDataSetFromMatrix(countData = tst ,  colData = coldata , design = ~0 + sex + plate + age + usableBases + PD_status)
pca_tst <- prcomp(t(dds@assays@data$counts))

library(ggbiplot)
ggbiplot(tst)
plot(pca_tst$x[,1] , pca_tst$x[,2])

tst <- estimateSizeFactors(count_mtx[keep_tst,])
pca_tst <- 

plot(tst$x[,1] , tst$x[,2])
