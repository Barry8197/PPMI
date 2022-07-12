setwd('~/PPMI_IP/R_Code/')

# Pull in Count Matrices --------------------------------------------------
TS <- 'V04'
print('Pull in Raw Feature Counts')
count_mtx <- NA
for (DS in c('Idiopathic_PD')) {
  print(DS)
  for (TS in list.files(paste0('./../Data/Count_Data/Raw_Counts_TS/',DS,'/'))[c(3)]) {
    print(TS)
    count_mtx_tmp <- read.csv(paste0('~/PPMI_IP/Data/Count_Data/Raw_Counts_TS/',DS,'/',TS,'/count_matrix.csv') , check.names = FALSE)
    #count_mtx <- cbind(count_mtx_tmp[sample(colnames(count_mtx_tmp[,2:length(count_mtx_tmp)]) , 190)] , count_mtx) 
    count_mtx <- cbind(count_mtx_tmp[,2:length(count_mtx_tmp)] , count_mtx) 
  }
}
rownames(count_mtx) <- count_mtx_tmp$Geneid #add gene ids to rownames
rm(count_mtx_tmp)
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

my_kmeans_clust <- kmeans(t(count_mtx) , 3)$cluster
(groups <- table(my_kmeans_clust)/length(my_kmeans_clust))
g <- c()
for (i in c(1,2,3)) {
  g <- c(g, names(sample(my_kmeans_clust[my_kmeans_clust == i] , groups[[i]]*172)))
}
count_mtx_iPD <- count_mtx[,g]

