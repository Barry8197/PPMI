library(DESeq2)
library(variancePartition)
library(edgeR) #expression normalization

# Pull in Count Matrices --------------------------------------------------
count_mtx <- NA
for (DS in c('Idiopathic_PD' , 'Genetic_PD' ,'Healthy_Control', 'Genetic_Unaffected','Prodromal' , 'SWEDD')) {
  for (i in list.files(paste0('PPMI_IP/Data_2/',DS,'/'))) {
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
meta_df <- read.csv('./PPMI_IP/meta_files/metadata_fianl.csv')

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

coldata <- data.frame( meta_df$HudAlphaSampleName, meta_df$PATNO, meta_df$Neutrophils,
                       meta_df$Lymphocytes, PD_status , meta_df$Genetic_Status , plate ,
                       visit, meta_df$Position , meta_df$PD_medication , meta_df$QC ,
                       sex , meta_df$Race , usableBases, age , row.names = names(count_mtx))

# Filter out non-expressed genes ------------------------------------------
keep <- filterByExpr(count_mtx, group = PD_status)
sum(keep) #18996 genes


dds = DESeqDataSetFromMatrix(countData = count_mtx[keep,], colData = coldata,  design = ~ 1)

# Estimate library size correction scaling factors
dds <- estimateSizeFactors(dds)
# identify genes that pass expression cutoff
isexpr <- rowSums(fpm(dds)>0.5) >= 0.5 * ncol(dds)
sum(isexpr)

quantLog <- log2( fpm( dds )[isexpr,] + 1)

colnames(coldata) <- c('Sample','Individual' , 'Neutrophils' , 'Lymphocytes',
                       'Disease' , 'Genetics', 'Plate' , 'Visit' , 'Position',
                       'PD_medication' , 'QC' , 'Sex' , 'Race' ,'Usable_bases' 
                       , 'Age')
coldata$Position <- as.factor(coldata$Position)
coldata$Genetics <- as.factor(coldata$Genetics)
coldata$QC <- as.factor(coldata$QC)

form <- ~   (1|Individual) + Neutrophils + Lymphocytes + 
  (1|Disease) + (1|Genetics) + (1|Plate) + (1|Visit) + 
  (1|Position) +  (1|QC)  + (1|Sex) + 
  (1|Race) + Usable_bases + (1|Age)
# Run variancePartition analysis
varPart <- fitExtractVarPartModel( quantLog, form, coldata)

vp <- sortCols( varPart )

plotPercentBars( vp[1:10,] )
#
# Figure 1b
# violin plot of contribution of each variable to total variance
plotVarPart( vp )
