library(tidyverse) #42; base for data/table manipulation
library(tidyr) #for pivot tables
library(edgeR) #expression normalization
library(limma) #more expression normalization
library(corrplot) #correlation plot matrix

# Load Expression Data --------------------------------------------------
load('./../Data/preprocessedData/filtered_raw_data.RData')


PD_status <- coldata$PD_status
sex <- coldata$sex
plate <- coldata$plate
age <- coldata$age
usableBases <- coldata$usableBases

design <- model.matrix(~ 0 + PD_status + sex + plate + age + usableBases)
colnames(design) <- gsub("PD_status", "", colnames(design))
colnames(design)

# Differential Expression using Limma Fit ---------------------------------
dge <- DGEList(count_mtx)
dge  <- calcNormFactors(dge)

logCPM <- cpm(dge, log = TRUE, prior.count = 3)


#make table of mean expression for each group
PDGroup <- meta_df$HudAlphaSampleName[meta_df$Disease.Status %in% c("Idiopathic PD", "Genetic PD")]
conGroup <- meta_df$HudAlphaSampleName[meta_df$Disease.Status %in% c("Healthy Control", "Genetic Unaffected")]


PDCPM <- logCPM[,PDGroup]
conCPM <- logCPM[,conGroup]

meanCPM <- data.frame("gene_id" = row.names(logCPM),
                      "AvgExpr_allPD" = rowMeans(PDCPM),
                      "AvgExpr_allControls" = rowMeans(conCPM))

v <- voom(dge, design, plot=TRUE) #check voom plot for curve to see if we need to do more filtering
vfit <- lmFit(v, design)

contr.matrix <- makeContrasts(
  PD = Affected - Unaffected,
  levels = colnames(design)
)

vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit)
summary(decideTests(efit))


# Identify differentially expressed genes and visualize -------------------
dt <- decideTests(efit)
contrName <- colnames(contr.matrix)[1]
print(contrName)
df <- rownames_to_column(topTable(efit, coef = 1, n = Inf), var = "gene_id")

genes_info = df %>% data.frame  %>% 
  mutate(significant=adj.P.Val<0.05 & !is.na(adj.P.Val) )

df$diffexpressed <- "NO"
# if log2Foldchange > 0.1 and adjusted pvalue < 0.05, set as "UP" 
df$diffexpressed[df$logFC > 0.1 & df$adj.P.Val < 0.05] <- "UP"
# if log2Foldchange < -0.1 and adjusted pvalue < 0.05, set as "DOWN"
df$diffexpressed[df$logFC < -0.1 & df$adj.P.Val < 0.05] <- "DOWN"
df$delabel <- NA
df$delabel[df$diffexpressed != "NO"] <- df$gene_id[df$diffexpressed != "NO"]

ggplot(data=df, aes(x=logFC, y=-log10(adj.P.Val), col=diffexpressed , label=delabel)) +
  geom_point() + 
  theme_minimal() +
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=c(-0.1, 0.1), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red") + 
  ggtitle('Case-Control Volcano Plot Limma-Fit') + theme(plot.title = element_text(hjust = 0.5)) + 
  ylim(c(0,15))

write.csv(df , file = './PPMI_IP/output/case_control_limma.csv')

