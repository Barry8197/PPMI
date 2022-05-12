library(edgeR) #expression normalization
library(devtools)
install_github("vqv/ggbiplot")
library(ggbiplot)
library(DESeq2)

# Pull in Count Matrices --------------------------------------------------
count_mtx <- NA
for (DS in c('Idiopathic_PD' , 'Genetic_PD' ,'Healthy_Control', 'Genetic_Unaffected')) {
  for (i in list.files(paste0('PPMI_IP/Data_2/',DS,'/'))[-c(3)]) {
    count_mtx_tmp <- read.csv(paste0('~/PPMI_IP/Data_2/',DS,'/',i,'/count_matrix.csv') , check.names = FALSE)
    count_mtx <- cbind(count_mtx_tmp[sample(colnames(count_mtx_tmp[,2:length(count_mtx_tmp)]) , 272)] , count_mtx) 
  }
}
rownames(count_mtx) <- count_mtx_tmp$Geneid #add gene ids to rownames

genes_removed <- read.csv('PPMI_IP/Data/genes_removed.csv' , header = FALSE)$V1 #remove genes as per paper
genes_subset <- setdiff(rownames(count_mtx) , genes_removed)
count_mtx <- count_mtx[genes_subset , ]
count_mtx <- count_mtx[,1:1632]


# Create coldata and condition table --------------------------------------
meta_df <- read.csv('./PPMI_IP/meta_files/metadata_fianl.csv')

i <- match(meta_df$HudAlphaSampleName, names(count_mtx), nomatch = 0)
count_mtx <- count_mtx[,i]
meta_df <- meta_df[which(meta_df$HudAlphaSampleName %in% names(count_mtx)),]

identical(names(count_mtx) , meta_df$HudAlphaSampleName) #doube check that ids and genes match

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


# Filter out non-expressed genes ------------------------------------------
keep <- filterByExpr(count_mtx, group = PD_status)
sum(keep) #18996 genes


# Obtain R Log  ------------------------------------
dds = DESeqDataSetFromMatrix(countData = count_mtx[keep,], colData = coldata , design = ~0  + sex + plate + age + usableBases + PD_status)
rld <- rlog(dds , blind = FALSE)


# Perform PCA -------------------------------------------------------------
coldata.pca <- prcomp(count_mtx[keep,] , scale = TRUE , center = TRUE)
summary(coldata.pca)


# Visualise PCA -----------------------------------------------------------
ggbiplot(pca)

