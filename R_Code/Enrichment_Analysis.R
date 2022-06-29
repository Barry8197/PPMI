library(org.Hs.eg.db)
library(clusterProfiler)
library(WGCNA)
library(dplyr)
library(biomaRt)
library(DOSE)
library(ReactomePA)

# Pull in Data and get Entrez ID's ----------------------------------------
top_modules <- c("#E08A00" , "#00BECA")
clusters <- read.csv('./../Data/preprocessedData/clusters.csv')
load('./../Data/preprocessedData/preprocessed_data.RData')

# Get entrez ID of genes
genes_info$ID <- rownames(genes_info)
genes_info <- genes_info %>% left_join(clusters %>% dplyr::rename(Module = Modules_dh), by='ID') %>%
  dplyr::select(ID, baseMean, log2FoldChange, significant, Module , padj) %>%
  mutate(module_number = Module %>% as.factor %>% as.numeric)
genes_info$ID <- gsub("\\..*","", genes_info$ID)
gene_names = genes_info %>% dplyr::rename('ensembl_gene_id' = ID) %>% filter(Module!='gray')

# Prepare dataset for Enrichment Analysis

EA_dataset = genes_info %>% dplyr::rename('ensembl_gene_id' = ID) %>% filter(Module!='gray')

# ClusterProfile works with Entrez Gene Ids, o we have to assign one to each gene
getinfo = c('ensembl_gene_id','entrezgene','hgnc_symbol')
mart=useMart(biomart='ENSEMBL_MART_ENSEMBL',dataset='hsapiens_gene_ensembl',host='feb2014.archive.ensembl.org')
biomart_output = getBM(attributes=getinfo, filters=c('ensembl_gene_id'), 
                       values=EA_dataset$ensembl_gene_id, mart=mart)

EA_dataset = biomart_output %>% left_join(EA_dataset, by='ensembl_gene_id') %>% dplyr::rename('ID'=ensembl_gene_id)

rm(getinfo, mart, biomart_output , gene_names)


# Enrichment Analysis -----------------------------------------------------
universe = EA_dataset$entrezgene %>% as.character
genes_in_module = EA_dataset %>% filter(Module == top_modules[1]) %>% pull(entrezgene)

ORA_GO = enrichGO(gene = genes_in_module, universe = universe, OrgDb = org.Hs.eg.db, ont = 'BP', 
                  pAdjustMethod = 'bonferroni')

goplot(ORA_GO)


ORA_KEGG = enrichKEGG(gene = genes_in_module, universe = universe, pAdjustMethod = 'bonferroni')

library(pathview)
logFC <- DE_info$log2FoldChange
names(logFC) <- DE_info$entrezgene
pathview(gene.data = logFC, 
         pathway.id = "hsa05417", 
         species = "hsa", 
         limit = list(gene=100, cpd=10))

DE_info <- EA_dataset %>% data.frame  %>% 
  mutate(DEG = (significant== TRUE & abs(log2FoldChange) > 0.1 ) )

library(fgsea)
gseaDat <- filter(DE_info, !is.na(entrezgene))

rankData <- gseaDat$log2FoldChange
names(rankData) <- gseaDat$entrezgene
head(rankData)

geneList <- gseaDat[which(gseaDat$module_number == 37),]$entrezgene
names(geneList) = geneList %>% pull(entrezgene) %>% as.character
geneList = sort(geneList, decreasing = TRUE)


GSEA_GO = gseGO(geneList, OrgDb = org.Hs.eg.db, pAdjustMethod = 'bonferroni', pvalueCutoff = 0.1, 
                nPerm = nPerm, verbose = FALSE, seed = TRUE)

fgseaRes <- fgsea(pathways = geneList, 
                  stats = rankData, 
                  minSize = 15, 
                  maxSize = 500)

fgseaRes %>% 
  arrange(desc(abs(NES))) %>% 
  top_n(10, -padj)

EA_dataset = genes_info %>% dplyr::rename('ensembl_gene_id' = ID) %>% filter(Module!='gray')

# ClusterProfile works with Entrez Gene Ids, o we have to assign one to each gene
getinfo = c('ensembl_gene_id','entrezgene')
mart=useMart(biomart='ENSEMBL_MART_ENSEMBL',dataset='hsapiens_gene_ensembl',host='feb2014.archive.ensembl.org')
biomart_output = getBM(attributes=getinfo, filters=c('ensembl_gene_id'), 
                       values=EA_dataset$ensembl_gene_id, mart=mart)

EA_dataset = biomart_output %>% left_join(EA_dataset, by='ensembl_gene_id') %>% dplyr::rename('ID'=ensembl_gene_id)

rm(getinfo, mart, biomart_output)


plotEnrichment(DE_info$entrezgene, rankData)

plot_data <- EA_dataset %>% filter(Module == top_modules)

ggplot(plot_data) + geom_boxplot(aes(x = as.factor(module_number) , y = baseMean))
# ORA enrichment
file_name = './../Data/preprocessedData/ORA_results.RData'
if(file.exists(file_name)){
  load(file_name)
} else {
  # Prepare input
  universe = EA_dataset$entrezgene %>% as.character
  
  # Perform Enrichment
  ORA_enrichment = list()
  
  for(module in top_modules){
    
    genes_in_module = EA_dataset %>% filter(Module == module) %>% pull(entrezgene)
    
    ORA_GO = enrichGO(gene = genes_in_module, universe = universe, OrgDb = org.Hs.eg.db, ont = 'All', 
                      pAdjustMethod = 'bonferroni', pvalueCutoff = 0.1, qvalueCutoff = 1)
    
    ORA_DO = enrichDO(gene = genes_in_module, universe = universe, qvalueCutoff = 1,
                      pAdjustMethod = 'bonferroni', pvalueCutoff = 0.1)
    
    ORA_DGN = enrichDGN(gene = genes_in_module, universe = universe, qvalueCutoff = 1,
                        pAdjustMethod = 'bonferroni', pvalueCutoff = 0.1)
    
    ORA_KEGG = enrichKEGG(gene = genes_in_module, universe = universe, qvalueCutoff = 1,
                          pAdjustMethod = 'bonferroni', pvalueCutoff = 0.1) 
    
    ORA_Reactome = enrichPathway(gene = genes_in_module, universe = universe, qvalueCutoff = 1,
                                 pAdjustMethod = 'bonferroni', pvalueCutoff = 0.1)
    
    ORA_enrichment[[module]] = list('GO' = ORA_GO, 'DO' = ORA_DO, 'DGN' = ORA_DGN, 'KEGG' = ORA_KEGG, 
                                    'Reactome' = ORA_Reactome)
    
    # Save after each iteration
    save(ORA_enrichment, file = file_name)
  }
  
  rm(universe, genes_in_module, module, ORA_GO, ORA_DGN, ORA_DO, ORA_KEGG, ORA_Reactome)
  
}
      
# Get shared enrichment for each module
top_modules_enrichment = list()

for(module in top_modules){
  
  ORA_enrichment_for_module = ORA_enrichment[[module]]
  
  for(dataset in c('KEGG', 'Reactome', 'GO', 'DO', 'DGN')){
    
    print(dataset)
    ORA_enrichment_dataset = ORA_enrichment_for_module[[dataset]] %>% data.frame %>%
      dplyr::rename('pvalue_ORA' = pvalue, 'p.adjust_ORA' = p.adjust, 'qvalue_ORA' = qvalue)
  
  }
  
  top_modules_enrichment[[module]] = ORA_enrichment_dataset  
}

save(top_modules_enrichment, file = './../Data/preprocessedData/top_modules_enrichment.RData')

rm(module, ORA_enrichment_for_module, dataset, 
   ORA_enrichment_dataset)
}

