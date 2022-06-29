library(DESeq2)
library(WGCNA)
library(dplyr)
library(ggplot2)
library(edgeR)
library(reshape2)
library(kableExtra)
library(plotly)
library(vsn)
library(tibble)

setwd('~/PPMI_IP/R_Code/')

# Pull in Count Matrices --------------------------------------------------
print('Pull in Raw Feature Counts')
count_mtx <- NA
for (DS in c('Idiopathic_PD' ,'Healthy_Control')) {
  print(DS)
  for (i in list.files(paste0('./../Data/Count_Data/Raw_Counts_TS/',DS,'/'))[c(5)]) {
    print(i)
    count_mtx_tmp <- read.csv(paste0('~/PPMI_IP/Data/Count_Data/Raw_Counts_TS/',DS,'/',i,'/count_matrix.csv') , check.names = FALSE)
    count_mtx <- cbind(count_mtx_tmp[sample(colnames(count_mtx_tmp[,2:length(count_mtx_tmp)]) , 150)] , count_mtx) 
    #count_mtx <- cbind(count_mtx_tmp[,2:length(count_mtx_tmp)] , count_mtx) 
  }
}
rownames(count_mtx) <- count_mtx_tmp$Geneid #add gene ids to rownames
rm(count_mtx_tmp)

# genes_removed <- read.csv('../Data/meta_files/genes_removed.csv' , header = FALSE)$V1 #remove genes as per paper
# genes_subset <- setdiff(rownames(count_mtx) , genes_removed)
# count_mtx <- count_mtx[genes_subset , ]
count_mtx <- count_mtx[,1:length(colnames(count_mtx))-1]

# Create coldata and condition table --------------------------------------
print('Create Coldata and Conditions')
meta_df <- read.csv('./../Data/meta_files/metadata_final.csv')

i <- match(meta_df$HudAlphaSampleName, names(count_mtx), nomatch = 0)
count_mtx <- count_mtx[,i]
meta_df <- meta_df[which(meta_df$HudAlphaSampleName %in% names(count_mtx)),]

identical(names(count_mtx) , meta_df$HudAlphaSampleName) #doube check that ids and genes match


# Create Coldata Matrix ---------------------------------------------------
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

coldata <- data.frame( meta_df$HudAlphaSampleName , meta_df$PATNO , meta_df$Genetic_Status , meta_df$Clinical_Event, sex , age ,plate , usableBases, PD_status ,  row.names = names(count_mtx))

rm(meta_df)

# Count Distribution ------------------------------------------------------
counts = count_mtx %>% melt

count_distr = data.frame('Statistic' = c('Min', '1st Quartile', 'Median', 'Mean', '3rd Quartile', 'Max'),
                         'Values' = c(min(counts$value), quantile(counts$value, probs = c(.25, .5)) %>% unname,
                                      mean(counts$value), quantile(counts$value, probs = c(.75)) %>% unname,
                                       max(counts$value)))

count_distr %>% kable(digits = 2, format.args = list(scientific = FALSE)) %>% kable_styling(full_width = F)

rm(counts, count_distr)


# Gene Annotation ---------------------------------------------------------
#Skipping for now


# Remove Genes with low level of expression -------------------------------
to_keep = rowSums(count_mtx) > 0 #removed 1157 genes
length(to_keep) - sum(to_keep)

count_mtx <- count_mtx[to_keep,]

datExpr <- count_mtx
datMeta <- coldata
datMeta_original = datMeta
datExpr_original = datExpr

# thresholds = round(100*(88-c(seq(72,7,-5),5,3,2,0))/88)
# 
# for(threshold in thresholds){
#   
#   datMeta = datMeta_original
#   datExpr = datExpr_original
#   
#   to_keep = apply(datExpr, 1, function(x) 100*mean(x>0)) >= threshold
#   datExpr = datExpr[to_keep,]
#   
#   # Filter outlier samples
#   absadj = datExpr %>% bicor %>% abs
#   netsummary = fundamentalNetworkConcepts(absadj)
#   ku = netsummary$Connectivity
#   z.ku = (ku-mean(ku))/sqrt(var(ku))
#   
#   to_keep = z.ku > -2
#   datMeta = datMeta[to_keep,]
#   datExpr = datExpr[,to_keep]
#   
#   rm(absadj, netsummary, ku, z.ku, to_keep)
#   
#   
#   # Create a DeseqDataSet object, estimate the library size correction and save the normalized counts matrix
#   counts = datExpr %>% as.matrix
# 
#   dds =  DESeqDataSetFromMatrix(countData = datExpr, colData = datMeta , design = ~PD_status)
#   
#   # Perform vst
#   vsd = vst(dds)
#   
#   datExpr_vst = assay(vsd)
#   datMeta_vst = colData(vsd)
#   datGenes_vst = rowRanges(vsd)
#   
#   rm(counts, rowRanges, se, vsd)
#   
#   # Save summary results in dataframe
#   if(threshold == thresholds[1]){
#     mean_vs_sd_data = data.frame('threshold'=threshold, 'ID'=rownames(datExpr_vst),
#                                  'Mean'=rowMeans(datExpr_vst), 'SD'=apply(datExpr_vst,1,sd))
#   } else {
#     new_entries = data.frame('threshold'=threshold, 'ID'=rownames(datExpr_vst),
#                              'Mean'=rowMeans(datExpr_vst), 'SD'=apply(datExpr_vst,1,sd))
#     mean_vs_sd_data = rbind(mean_vs_sd_data, new_entries)
#   }
# }
# 
# # Visualise the effects of different thresholds
# to_keep_1 = mean_vs_sd_data$ID[mean_vs_sd_data$threshold==thresholds[1] & mean_vs_sd_data$Mean<7] %>%
#   as.character
# to_keep_2 = mean_vs_sd_data$ID[mean_vs_sd_data$threshold==thresholds[1] & mean_vs_sd_data$Mean>=7]
# to_keep_2 = to_keep_2 %>% sample(round(length(to_keep_2)/10)) %>% as.character
# 
# plot_data = mean_vs_sd_data[mean_vs_sd_data$ID %in% c(to_keep_1, to_keep_2),]
# 
# ggplotly(plot_data %>% ggplot(aes(Mean, SD)) + 
#            geom_point(color='#0099cc', alpha=0.2, aes(id=ID, frame=threshold)) + 
#            scale_x_log10() + scale_y_log10() + theme_minimal())
# 
# # Plot remaining genes
# plot_data = mean_vs_sd_data %>% group_by(threshold) %>% tally
# 
# ggplotly(plot_data %>% ggplot(aes(threshold, n)) + geom_point() + geom_line() +
#            theme_minimal() + ggtitle('Remaining genes for each filtering threshold'))
# 
# rm(to_keep_1, to_keep_2, plot_data, dds, thresholds , mean_vs_sd_data , new_entries)
# 
# # Return to original variables
# datExpr = datExpr_original
# datMeta = datMeta_original
# 
# rm(datExpr_vst, datMeta_vst)

# Minimum percentage of non-zero entries allowed per gene
#threshold = 99.7

#to_keep = apply(datExpr, 1, function(x) 100*mean(x>0)) >= threshold
to_keep = filterByExpr(datExpr , group = PD_status)
print('keeping genes')
print(sum(to_keep))
#length(to_keep) - sum(to_keep)
count_mtx = count_mtx[to_keep,]

rm(datExpr_original, datMeta_original, threshold, to_keep , datExpr , datMeta)


# Remove Outliers ---------------------------------------------------------
print('removing outliers')
absadj = count_mtx %>% bicor %>% abs
netsummary = fundamentalNetworkConcepts(absadj)
ku = netsummary$Connectivity
z.ku = (ku-mean(ku))/sqrt(var(ku))

plot_data = data.frame('sample'=1:length(z.ku), 'distance'=z.ku, 'Sample_ID'=coldata$meta_df.HudAlphaSampleName,
                       'Subject_ID'=coldata$meta_df.PATNO, 'Plate'=coldata$plate,
                       'Usable Bases'=coldata$usableBases, 'Sex'=coldata$sex, 'Age'=coldata$age,
                       'Diagnosis'=coldata$PD_status)

ggplot(plot_data) + geom_point(aes(distance , Subject_ID , color = Diagnosis )  ) +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(),axis.ticks.y=element_blank())

to_keep = z.ku > -2
sum(to_keep)
length(to_keep) - sum(to_keep)

count_mtx <- count_mtx[,to_keep] #removed 95
coldata <- coldata[to_keep,]

rm(absadj, netsummary, ku, z.ku, plot_data , to_keep)

save(count_mtx, coldata, file='./../Data/preprocessedData/filtered_raw_data.RData')
load('./../Data/preprocessedData/filtered_raw_data.RData')


# Normalisation Using DESeq -----------------------------------------------
plot_data = data.frame('ID'=rownames(count_mtx), 'Mean'=rowMeans(count_mtx), 'SD'=apply(count_mtx,1,sd))

plot_data %>% ggplot(aes(Mean, SD)) + geom_point(color='#0099cc', alpha=0.1) + geom_abline(color='gray') +
  scale_x_log10() + scale_y_log10() + theme_minimal()

rm(plot_data)

dds = DESeqDataSetFromMatrix(countData = count_mtx, colData = coldata , design = ~ usableBases + sex + age + plate  + PD_status)

print('performing DESeq')
dds = DESeq(dds)

# DEA Plots ---------------------------------------------------------------
DE_info = results(dds)
DESeq2::plotMA(DE_info, main= 'Original LFC values')

#DE_info_shrunken = lfcShrink(dds, coef='PD_statusAffected', res = DE_info, lfcThreshold=0)
#plotMA(DE_info_shrunken, main = 'Shrunken LFC values')

rm(DE_info , DE_info_shrunken)
# VST Transformation of Data ----------------------------------------------
vsd = vst(dds)

datExpr_vst = assay(vsd)
datMeta_vst = colData(vsd)

rm(vsd)

meanSdPlot(datExpr_vst, plot=FALSE)$gg + theme_minimal() + ylim(c(0,2))

plot_data = data.frame('ID'=rownames(datExpr_vst), 'Mean'=rowMeans(datExpr_vst), 'SD'=apply(datExpr_vst,1,sd))

plot_data %>% ggplot(aes(Mean, SD)) + geom_point(color='#0099cc', alpha=0.2) + geom_smooth(color = 'gray') +
  scale_x_log10() + scale_y_log10() + theme_minimal()

rm(plot_data)


# Save Expression & Meta data ---------------------------------------------
datExpr = datExpr_vst
datMeta = datMeta_vst %>% data.frame

rm(datExpr_vst, datMeta_vst , count_mtx , coldata)

genes_info = DE_info %>% data.frame  %>% 
  mutate(significant=padj<0.05 & !is.na(padj) )
#%>%
#   cbind(DE_info_shrunken %>% data.frame %>% dplyr::select(log2FoldChange, lfcSE) %>% 
#           dplyr::rename('shrunken_log2FoldChange' = log2FoldChange, 'shrunken_lfcSE' = lfcSE)) %>%
#   add_column('ID'=rownames(DE_info), .before='baseMean')

counts = datExpr %>% data.frame %>% melt

count_distr = data.frame('Statistic' = c('Min', '1st Quartile', 'Median', 'Mean', '3rd Quartile', 'Max'),
                          'Values' = c(min(counts$value), quantile(counts$value, probs = c(.25, .5)) %>% unname,
                                       mean(counts$value), quantile(counts$value, probs = c(.75)) %>% unname,
                                       max(counts$value)))

count_distr %>% kable(digits = 2, format.args = list(scientific = FALSE)) %>% kable_styling(full_width = F)

save(datExpr, datMeta, dds, genes_info, file='./../Data/preprocessedData/preprocessed_gPD_BL.RData')
