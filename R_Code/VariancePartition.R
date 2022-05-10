library(DESeq2)
library(variancePartition)

count_mtx_PD <- NA
for (DS in c('Idiopathic_PD','Genetic_PD')) {
  for (i in c('BL' , 'V02' , 'V04' , 'V06' , 'V08')) {
    if (length(colnames(count_mtx_PD)) == 0) {
      count_mtx_PD <- read.csv(paste0('~/PPMI_IP/Data/',DS,'/',i,'/count_matrix.csv') , check.names = FALSE)
    }
    else {
      count_mtx_PD <- cbind(count_mtx_PD , read.csv(paste0('~/PPMI_IP/Data/',DS,'/',i,'/count_matrix.csv') , check.names = FALSE)) 
    }
  }
}

count_mtx_HC <- NA
for (HC in c('Healthy_Control' , 'Genetic_Unaffected')) {
  for (i in c('BL' , 'V02' , 'V04' , 'V06' , 'V08')) {
    if (length(colnames(count_mtx_PD)) == 0) {
      count_mtx_HC <- read.csv(paste0('~/PPMI_IP/Data/',HC,'/',i,'/count_matrix.csv') , check.names = FALSE)
    }
    else {
      count_mtx_HC <- cbind(count_mtx_HC , read.csv(paste0('~/PPMI_IP/Data/',HC,'/',i,'/count_matrix.csv') , check.names = FALSE)) 
    }
  }
}
count_mtx <- cbind(count_mtx_PD , count_mtx_HC)

gene_ids = count_mtx$Geneid
rownames(count_mtx) <- gene_ids
count_mtx = count_mtx[,-1] 

meta_df <- read.csv('./PPMI_IP/Data/meta_data.11192021.csv')
coldata <- meta_df[which(meta_df$HudAlphaSampleName %in% colnames(count_mtx)) , ]
coldata <- coldata[,c('HudAlphaSampleName' , 'PATNO','Neutrophils..GI.L.' , 'Lymphocytes..GI.L.', 'Disease.Status', 'Genetic.Status','Plate','Clinical.Event','Position','Hours.since.PD.Med','Sex' , 'Race','QC','Age.at.Consent')]
coldata$PATNO <- gsub("PP-", "", coldata$PATNO, fixed = TRUE)
rownames(coldata) <- coldata$HudAlphaSampleName

count_mtx <- count_mtx[, rownames(coldata)]

dds = DESeqDataSetFromMatrix(countData = count_mtx, colData = coldata,  design = ~ 1)

# Estimate library size correction scaling factors
dds <- estimateSizeFactors(dds)
# identify genes that pass expression cutoff
isexpr <- rowSums(fpm(dds)>0.5) >= 0.5 * ncol(dds)

quantLog <- log2( fpm( dds )[isexpr,] + 1)

colnames(coldata) <- c('Sample','Individual' , 'Neutrophils' , 'Lymphocytes', 'Disease' , 'Genetics', 'Plate' , 'Visit' , 'Position', 'PD_medication' , 'Sex' , 'Race' ,'Usable_bases' , 'Age')
coldata$Position <- as.factor(coldata$Position)

form <- ~   (1|Individual) + Neutrophils + Lymphocytes + (1|Disease) + (1|Genetics) + (1|Plate) + (1|Visit) + (1|Position) + PD_medication  + (1|Sex) + (1|Race) + (1|Usable_bases) + Age
# Run variancePartition analysis
varPart <- fitExtractVarPartModel( quantLog, form, coldata)

vp <- sortCols( varPart )

plotPercentBars( vp[1:10,] )
#
# Figure 1b
# violin plot of contribution of each variable to total variance
plotVarPart( vp )
