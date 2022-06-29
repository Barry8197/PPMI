library(DESeq2)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(GenomicSuperSignature)
library(scran)
library(biomaRt)

setwd('~/PPMI_IP/R_Code/')

# Load Normalised Dataset -------------------------------------------------
load('./../Data/preprocessedData/preprocessed_data.RData')

RAVmodel <- getModel("PLIERpriors", load=TRUE)

# GenomicSS works with Entrez Gene Ids, so we have to assign one to each gene
rownames(datExpr) <- gsub("\\..*","", rownames(datExpr))
getinfo = c('ensembl_gene_id','hgnc_symbol')
mart=useMart(biomart='ENSEMBL_MART_ENSEMBL',dataset='hsapiens_gene_ensembl',host='feb2014.archive.ensembl.org')
biomart_output = getBM(attributes=getinfo, filters=c('ensembl_gene_id'), 
                       values=rownames(datExpr), mart=mart)

datExpr <- datExpr[biomart_output$ensembl_gene_id,]
rownames(datExpr) <- biomart_output$hgnc_symbol

val_all <- validate(datExpr, RAVmodel)

heatmapTable(val_all, RAVmodel, num.out = 5, swCutoff = 0)

plotValidate(val_all, interactive = FALSE)

validated_ind <- validatedSignatures(val_all, RAVmodel, num.out = 3, 
                                     swCutoff = 0, indexOnly = TRUE)
validated_ind
drawWordcloud(RAVmodel, validated_ind[1])
drawWordcloud(RAVmodel, validated_ind[2])
drawWordcloud(RAVmodel, validated_ind[3])

RAVnum <- validated_ind[2]  # RAV1139
res <- gsea(RAVmodel)[[RAVnum]]   
head(res)

findSignature(RAVmodel, "Parkinson")
findSignature(RAVmodel, "Parkinson", k = 1)

findKeywordInRAV(RAVmodel, "Parkinson", ind = 131)

subsetEnrichedPathways(RAVmodel, ind = RAVnum, n = , both = TRUE)
subsetEnrichedPathways(RAVmodel, ind = 91, n = 3, both = TRUE)
subsetEnrichedPathways(RAVmodel, ind = 1376, n = 3, both = TRUE)

findStudiesInCluster(RAVmodel, validated_ind[2])

findStudiesInCluster(RAVmodel, validated_ind[2], studyTitle = TRUE)

RAVnum <- 131  # RAV131
res <- gsea(RAVmodel)[[RAVnum]]   
head(res)

val_PD <- val_all[findSignature(RAVmodel, "Parkinson", k = 1),]
heatmapTable(val_PD, RAVmodel, num.out = 5, swCutoff = 0)
rav_131 <- subsetEnrichedPathways(RAVmodel, ind = 131, n = 3, both = TRUE)
