library(WGCNA)
library(DESeq2)
library(dplyr)
library(expss)
library(gridExtra)
library(RColorBrewer) 
library(dendextend)
library(ggplot2)

# Get colors from the ggplot palette
gg_colour_hue = function(n) {
  hues = seq(15, 375, length = n+1)
  pal = hcl(h = hues, l = 65, c = 100)[1:n]
}

# Assign an HCL rainbow colour to each module
get_mod_colours = function(mods){
  
  n = length(unique(mods))-1
  set.seed(123) ; rand_order = sample(1:n)
  mod_colors = c('white', gg_colour_hue(n)[rand_order])
  names(mod_colors) = mods %>% table %>% names
  
  return(mod_colors)
}

setwd('~/PPMI_IP/R_Code/')

load('./../Data/preprocessedData/preprocessed_gPD_BL.RData')

# Obtain Scale Free Topology ----------------------------------------------
# Choose a set of soft-thresholding powers
powers = c(seq(from=1 , to=5 , by = 1) , seq(from = 5.5 , to =7.5 , by = 0.2) , seq(from =7 , to = 10 , by = 1) )
# Call the network topology analysis function
sft = pickSoftThreshold(t(datExpr), powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.80,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


best_power = datExpr %>% t %>% pickSoftThreshold(powerVector = powers, RsquaredCut=0.8)
S_sft = datExpr %>% t %>% adjacency(type='signed hybrid', power=6.1, corFnc='bicor')

# Create Distance Matrix --------------------------------------------------
dissTOM = S_sft %>% TOMdist

rownames(dissTOM) = rownames(S_sft)
colnames(dissTOM) = colnames(S_sft)

rm(S_sft)


# Create Hierarchical Cluster Tree ----------------------------------------
dend = dissTOM %>% as.dist %>% hclust(method='average')
plot(dend, hang=0, labels=FALSE, main=NULL, sub='', xlab='')

dend_complete = dissTOM %>% as.dist %>% hclust
plot(dend_complete, hang=0, labels=FALSE, main=NULL, sub='', xlab='')


# Cluster Extraction ------------------------------------------------------
modules_dt = cutreeDynamic(dend, method = 'tree', minClusterSize = 10)

table(modules_dt) %>% data.frame %>% filter(modules_dt != 0) %>%
  ggplot(aes(x=modules_dt, y=Freq)) + geom_bar(stat='identity', color = 'white', fill = '#009999') + 
  xlab('Modules') + ylab('Number of genes') +
  theme_minimal() + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())

sum(modules_dt == 0)
#1846 genes without a Cluster

pca = datExpr %>% prcomp

plot_data = data.frame('ID'=rownames(datExpr), 'PC1' = pca$x[,1], 'PC2' = pca$x[,2], 
                       'hasCluster' = modules_dt!=0) %>%
  left_join(de , by = 'ID') %>% apply_labels(significant = 'padj < 0.05' , hasCluster = 'Belongs to a Module')


p1 = plot_data %>% ggplot(aes(PC1, PC2, color=hasCluster)) + geom_point(alpha=0.2) + 
  theme_minimal() + ggtitle('Genes are assigned to a cluster') + theme(legend.position='bottom')

p1
p2 = plot_data %>% ggplot(aes(PC1, PC2, color=significant)) + geom_point(alpha=0.2) + 
  theme_minimal() + ggtitle('Genes found to be DE') + theme(legend.position='bottom')

grid.arrange(p1, p2, nrow=1, top = 'Exploring unassigned genes using the Dynamic Tree algorithm')

cro(plot_data$significant, list(plot_data$hasCluster, total()))

modules_dh = cutreeDynamic(dend, minClusterSize = 10, distM = dissTOM)

table(modules_dh) %>% data.frame %>% filter(modules_dh != 0) %>%
  ggplot(aes(x=modules_dh, y=Freq)) + geom_bar(stat='identity', color = 'white', fill = '#009999') + 
  xlab('Modules') + ylab('Number of genes') + 
  theme_minimal() + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())


# Merge Similar Clusters --------------------------------------------------
# Calculate eigengenes
MEList = datExpr %>% t %>% moduleEigengenes(colors = modules_dt)
MEs = MEList$eigengenes

# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs)

# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = 'average')

METree %>% as.dendrogram %>% plot(leaflab = 'none')
abline(h=0.4, col='#009999')

merge_top_dt = datExpr %>% t %>% mergeCloseModules(modules_dt, cutHeight = 1)

merge_similar_dt = datExpr %>% t %>% mergeCloseModules(modules_dt, cutHeight = 0.4)

table(merge_similar_dt$colors) %>% data.frame %>% filter(Var1 != 0) %>%
  ggplot(aes(x=reorder(Var1, -Freq), y=Freq)) + geom_bar(stat='identity', color = 'white', fill = '#009999') + 
  xlab('Modules') + ylab('Number of genes') + 
  theme_minimal() + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())

rm(MEList, MEs, MEDiss, METree)

MEList = datExpr %>% t %>% moduleEigengenes(colors = modules_dh)
MEs = MEList$eigengenes

# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs)

# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = 'average')

METree %>% as.dendrogram %>% plot(leaflab = 'none')

merge_top_dh = datExpr %>% t %>% mergeCloseModules(modules_dh, cutHeight = 1)

rm(MEList, MEs, MEDiss, METree)

dend_colors_all = data.frame('ID' = rownames(datExpr),
                             'cl_dh' = get_mod_colours(modules_dh)[as.character(modules_dh)],
                             'top_dh' = get_mod_colours(merge_top_dh$colors)[as.character(merge_top_dh$colors)],
                             'dh' = rep('white', nrow(datExpr)),
                             'cl_dt' = get_mod_colours(merge_similar_dt$colors)[as.character(merge_similar_dt$colors)],
                             'top_dt' = get_mod_colours(merge_top_dt$colors)[as.character(merge_top_dt$colors)],
                             'dt' = rep('white', nrow(datExpr)))

dend %>% as.dendrogram(hang=0) %>% plot(ylim=c(min(dend$height),1), leaflab='none')
colored_bars(colors=dend_colors_all[dend$order,-1], cex.rowLabels = 0.8, text_shift = 0,
             rowLabels = rev(c('Dynamic Tree:','Top clusters', 'Clusters', 'Dynamic Hybrid:', 'Top clusters', 
                               'Clusters')))

# Save Clusters -----------------------------------------------------------
dt_info = data.frame('ID' = rownames(datExpr),
                     'Modules_dt' = get_mod_colours(modules_dt)[as.character(modules_dt)],
                     'Merged_dt' = get_mod_colours(merge_similar_dt$colors)[as.character(merge_similar_dt$colors)],
                     'Top_dt' = get_mod_colours(merge_top_dt$colors)[as.character(merge_top_dt$colors)])

dh_info = data.frame('ID' = rownames(datExpr),
                     'Modules_dh' = get_mod_colours(modules_dh)[as.character(modules_dh)],
                     'Top_dh' = get_mod_colours(merge_top_dh$colors)[as.character(merge_top_dh$colors)])

clusters = dt_info %>% left_join(dh_info, by='ID') %>% lapply(function(x) gsub('white','gray',x))

write.csv(clusters, file='./../Data/preprocessedData/clusters.csv', row.names=FALSE)

