library(dplyr)
library(DESeq2)

dir_GPD <- './PPMI_IP/Data/Genetic_PD/BL/'
filenames_genetic_PD <- file.path(dir_GPD , list.files(dir_GPD))
dir_HC <- './PPMI_IP/Data/Healthy_Control/BL/'
filenames_HC <- file.path(dir_HC , list.files(dir_HC))
filenames_all <- c(filenames_genetic_PD , filenames_HC)

# Function to read featureCounts output and extract the raw counts from it:
extract_counts = function(filename){
  
  # Load sample info
  sample_info = read.delim(filename, comment.char='#')
  
  # Rename column with counts with the sample ID 
  # I don't know which part of the whole name is the ID, so I just kept all the weird stuff from the name
  sample_id = unlist(strsplit(filename , '[.]'))[6]
  colnames(sample_info)[ncol(sample_info)] = sample_id
  
  # Extract the raw counts. Each row corresponds to a gene identified by its Ensembl ID (Geneid)
  sample_counts = sample_info %>% dplyr::select(Geneid, all_of(sample_id))
  
  return(sample_counts)
}

################################################################################################################

# To create the count matrix:

# If all the samples contain exactly the same number of genes in the same order, you can just use cbind
# but in case they don't, I'll use join

# 1. Initiating the matrix with the first sample keeping the Geneids as the first column
sample_counts = extract_counts(filenames_all[1])
count_matrix = sample_counts

# 2. Adding columns iteratively:
for(i in 2:length(filenames_all)){
  (print(filenames_all[i]))
  sample_counts = extract_counts(filenames_all[i])
  count_matrix = count_matrix %>% full_join(sample_counts, by = 'Geneid')
}

# Convert the dataframe to matrix, saving the gene IDs in another vector
gene_ids = count_matrix$Geneid
count_matrix = count_matrix[,-1] %>% as.matrix
################################################################################################################

# Your output is now a regular R data.frame, so you can do plots, filter rows or columns or anything you want with it

# You could do your filtering, cleaning, batch effects modelling, initial visualisations, etc here

################################################################################################################

meta_df <- read.csv('./PPMI_IP/Data/meta_data.11192021.csv')
coldata <- meta_df[which(meta_df$HudAlphaSampleName %in% colnames(count_matrix)) , ]
coldata <- coldata[,c('HudAlphaSampleName','Disease.Status', 'Race')]
rownames(coldata) <- coldata$HudAlphaSampleName
colnames(coldata)[2] <- 'condition'
coldata$condition <- gsub(" ", "_", coldata$condition, fixed = TRUE)
coldata$condition <- factor(coldata$condition)
count_matrix <- count_matrix[, rownames(coldata)]
# To create a DESeq data set you would input it into the function like this:
# The only thing missing is colData, which would be the information of each sample, including condition as a column
# If you do batch effect modelling, you should include the batch effects into the design of the model
dds = DESeqDataSetFromMatrix(countData = count_matrix, colData = coldata,  design = ~ condition)

# To do the differential expression analysis (DEA) with this dataset you only do:
dds = DESeq(dds)

# And to extract the results from the DEA:

# If you want regular log-fold change values
DE_info = results(dds, lfcThreshold = 0, altHypothesis = 'greaterAbs')

# If you want shrunken log-fold change values
DE_info_shrunken = lfcShrink(dds, coef = 'Diagnosis_ASD_vs_CTL', lfcThreshold = 0, res = DE_info, quiet = TRUE)

# Notes:
# - I chose a threshold of 0, which gives you a lot of DE genes, if you  want you can increase the threshold to be more strict about this
# - Regular vs shrunken LFC: There's usually a relation between LFC and the mean expression of genes, so shrunken LFC corrects for this

library(ggplot2)
de <- as.data.frame(DE_info)
de$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
de$diffexpressed[de$log2FoldChange > 0.6 & de$pvalue < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
de$diffexpressed[de$log2FoldChange < -0.6 & de$pvalue < 0.05] <- "DOWN"
p <- ggplot(data=de, aes(x=log2FoldChange, y=-log10(pvalue))) + geom_point() + theme_minimal()
p2 <- p + geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")
mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")
p3 <- p2 + scale_colour_manual(values = mycolors)
p3

ggplot(data=de, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed)) +
  geom_point() + 
  theme_minimal() +
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")
