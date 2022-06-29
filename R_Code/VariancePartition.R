library(DESeq2)
library(variancePartition)
library(edgeR) #expression normalization

# Pull in Pre-processed Data --------------------------------------------------
load('./../Data/preprocessedData/preprocessed_data.RData')

# Create design matrix from meta data -------------------------------------
meta_df <- read.csv('./../Data/meta_files/metadata_final.csv')

i <- match(meta_df$HudAlphaSampleName, colnames(datExpr), nomatch = 0)
datExpr <- datExpr[,i]
meta_df <- meta_df[which(meta_df$HudAlphaSampleName %in% colnames(datExpr)),]

identical(colnames(datExpr) , meta_df$HudAlphaSampleName) #doube check that ids and genes match

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
                       sex , meta_df$Race , usableBases, age , row.names = names(datExpr))


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
form <- ~   (1|Individual) + Usable_bases + (1|PD_status)

# Run variancePartition analysis
varPart <- fitExtractVarPartModel( assay(dds), form, coldata)

vp <- sortCols( varPart )

plotPercentBars( vp[1:10,] )
#
# Figure 1b
# violin plot of contribution of each variable to total variance
plotVarPart( vp )
