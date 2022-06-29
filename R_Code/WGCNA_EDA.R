library(WGCNA)
library(dplyr)
library(DESeq2)
library(polycor)
library(plotly)
library(doParallel)
library(gridExtra)

# Pull in Data ------------------------------------------------------------
clusters <- read.csv('./../Data/preprocessedData/clusters.csv')
load('./../Data/preprocessedData/preprocessed_gPD_BL.RData')
DE_info <- results(dds, lfcThreshold=0, altHypothesis='greaterAbs')
de_genes <- (DE_info %>% data.frame %>% mutate(de = (padj < 0.05) & (abs(log2FoldChange) > 0.1)))$de
#datExpr <- datExpr[which(de_genes == TRUE),]

genes_info$ID <- rownames(genes_info)
genes_info = genes_info %>% 
  left_join(clusters %>% dplyr::rename(Module = Modules_dh), by='ID') %>%
  dplyr::select(ID, baseMean, log2FoldChange, Module , significant) %>%
  mutate(module_number = Module %>% as.factor %>% as.numeric)
genes_info <- genes_info[de_genes ,]

datTraits <- datMeta %>% dplyr::select(PD_status ,sex , age ,  plate , usableBases) %>% mutate(across(where(is.factor), as.numeric)) 

# Module-Trait Associations -----------------------------------------------
# Define numbers of genes and samples
nGenes = nrow(datExpr);
nSamples = ncol(datExpr);

# Recalculate MEs with color labels
ME_object = datExpr %>% t %>% moduleEigengenes(colors = genes_info$Module)
MEs = orderMEs(ME_object$eigengenes)

# Calculate correlation between eigengenes and the traits and their p-values
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

# In case there are any NAs
if(sum(!complete.cases(moduleTraitCor))>0){
  print(paste0(sum(is.na(moduleTraitCor)),' correlation(s) could not be calculated')) 
}

rm(ME_object)

# Sort moduleTraitCor by Diagnosis
moduleTraitPvalue <- moduleTraitPvalue[order(moduleTraitCor[,1], decreasing=TRUE),]
moduleTraitCor <- moduleTraitCor[order(moduleTraitCor[,1], decreasing=TRUE),]

# Create text matrix for the Heatmap
textMatrix = paste0(signif(moduleTraitCor, 2), ' (', signif(moduleTraitPvalue, 1), ')')
dim(textMatrix) = dim(moduleTraitCor)

# Change name of gray cluster to reflect it's not a cluster
yLabels = paste0('Cluster ', rownames(moduleTraitCor) %>% as.factor %>% as.numeric)
yLabels[yLabels == paste0('Cluster ', ncol(MEs))] = 'No cluster'

labeledHeatmap(Matrix = moduleTraitCor, xLabels = names(datTraits), yLabels = yLabels, 
               yColorWidth = 0, colors = greenWhiteRed(50),
               bg.lab.y = gsub('ME', '', rownames(moduleTraitCor)), xLabelsPosition = 'top', xLabelsAngle = 0,
               textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 0.8, cex.lab.y = 0.75, zlim = c(-1,1))

pd_status_cor = data.frame('Module' = gsub('ME','',rownames(moduleTraitCor)),
                           'MTcor' = moduleTraitCor[,'PD_status'],
                           'MTpval' = moduleTraitPvalue[,'PD_status'])

genes_info = genes_info %>% left_join(pd_status_cor, by='Module')

rm(moduleTraitCor, moduleTraitPvalue, textMatrix, pd_status_cor, yLabels)


# Correlation between DE genes and PD_Status ------------------------------
plot_data = genes_info %>% group_by(Module, MTcor) %>% summarise(p = 100*mean(significant))

ggplotly(plot_data %>% ggplot(aes(MTcor, p)) + geom_smooth(color='gray', se=FALSE) +
           geom_point(color=plot_data$Module, alpha=0.7, aes(id=Module)) +
           geom_hline(yintercept=mean(plot_data$p),color='gray',linetype='dotted') +
           theme_minimal() + 
           xlab('Cluster-diagnosis correlation') + ylab('Percentage of DE genes'))

rm(plot_data)


# Gene-Level Metrics ------------------------------------------------------

############# 1. Calculate Gene Significance
GS_info = data.frame('ID'=rownames(datExpr),
                     'GS'=datExpr %>% apply(1, function(x) hetcor(x,datMeta$PD_status)$correlations[1,2])) %>%
  mutate('GSpval'=corPvalueStudent(GS, ncol(datExpr)))

#############  2. Calculate Module Membership

#setup parallel backend to use many processors
cores = detectCores()
cl = makeCluster(cores-1)
registerDoParallel(cl)

# Create matrix with MM by gene
MM = foreach(i=1:nrow(datExpr), .combine=rbind) %dopar% {
  library(polycor)
  tempMatrix = apply(MEs, 2, function(x) hetcor(as.numeric(datExpr[i,]), x)$correlations[1,2])
  tempMatrix
}

# Stop clusters
stopCluster(cl)

rownames(MM) = rownames(datExpr)
colnames(MM) = paste0('MM', gsub('ME','',colnames(MEs)))

# Calculate p-values
MMpval = MM %>% corPvalueStudent(ncol(datExpr)) %>% as.data.frame
colnames(MMpval) = paste0('MMpval', gsub('ME','',colnames(MEs)))

MM = MM %>% as.data.frame %>% mutate(ID = rownames(.))
MMpval = MMpval %>% as.data.frame %>% mutate(ID = rownames(.))

# Join and save results
dataset = genes_info %>% dplyr::select(ID, Module, module_number, MTcor, MTpval) %>%
  left_join(GS_info, by='ID') %>%
  left_join(MM, by='ID') %>%
  left_join(MMpval, by='ID')

write.csv(dataset, file = './../Data/preprocessedData/WGCNA_metrics.csv', row.names = FALSE)

rm(cores, cl) 


# Association between gene significance and cluster-Pd_status -------------
plot_data = dataset %>% dplyr::select(ID, MTcor, GS) %>% mutate(mean_expr = rowMeans(datExpr)) %>%
  left_join(genes_info %>% dplyr::select(ID,), by='ID') %>%
  left_join(genes_info %>% dplyr::select(ID, log2FoldChange, significant, 
                                         Module), by='ID') %>%
  left_join(data.frame(MTcor=unique(dataset$MTcor)) %>% arrange(by=MTcor) %>% 
              mutate(order=1:length(unique(dataset$MTcor))), by='MTcor')

plot_data %>% ggplot(aes(MTcor, GS)) + geom_point(color = plot_data$Module, alpha = 0.1) + 
  geom_smooth(color='gray', se = FALSE) + theme_minimal() + 
  xlab('Cluster-PD status correlation') + ylab('Gene Significance') + 
  ggtitle(paste0('R^2=',round(cor(plot_data$MTcor, plot_data$GS)^2,2)))

p1 = plot_data %>% ggplot(aes(MTcor, log2FoldChange)) + geom_point(color=plot_data$Module, alpha=0.1) + 
  geom_smooth(color='gray', se=FALSE) + xlab('Cluster-PD_status correlation') + ylab('LFC') +
  theme_minimal() + ggtitle(paste0('R^2 = ', round(cor(plot_data$log2FoldChange, plot_data$MTcor)[1]**2,2)))

p2 = plot_data %>% ggplot(aes(GS, log2FoldChange)) + geom_point(color=plot_data$Module, alpha=0.1) + 
  geom_smooth(color='gray', se=FALSE) + xlab('Gene Significance') + ylab('LFC') +
  theme_minimal() + ggtitle(paste0('R^2 = ', round(cor(plot_data$log2FoldChange, plot_data$GS)[1]**2,2)))

grid.arrange(p1, p2, nrow=1)

p1 = plot_data %>% ggplot(aes(MTcor, mean_expr)) + 
  geom_point(color=plot_data$Module, alpha=0.1) + geom_smooth(color='gray', se = TRUE, alpha = 0.2) + 
  xlab('Cluster-pd_status correlation') + ylab('Mean level of expression') + theme_minimal() + 
  ggtitle(paste0('R^2 = ', round(cor(plot_data$log2FoldChange, plot_data$MTcor)[1]**2,2)))

p2 = plot_data %>% ggplot(aes(GS, mean_expr)) + 
  geom_point(color=plot_data$Module, alpha=0.1) + geom_smooth(color='gray', se = TRUE, alpha = 0.2) + 
  xlab('Gene Significance') + ylab('Mean level of expression') + theme_minimal() + 
  ggtitle(paste0('R^2 = ', round(cor(plot_data$log2FoldChange, plot_data$GS)[1]**2,2)))

grid.arrange(p1, p2, nrow=1)

rm(p1,p2)

plot_ly(plot_data, x = ~MTcor, y = ~GS, z = ~mean_expr, color= ~Module, colors = sort(unique(plot_data$Module)), 
        alpha = 0.5, size = 0.5)


plot_data <- cluster_3 %>% mutate(mean_expr = rowMeans(cluster_3))
