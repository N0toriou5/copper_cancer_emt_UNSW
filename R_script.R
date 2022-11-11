### Analysis with DESeq2 MB231 ----------------------------------------------------
## This code part allows to produce the images from the Copper EMT paper
## First part is focused on the analysis of breast cancer 
# CODE BY Dr. Daniele Mercatelli, University of Bologna

setwd("C:/set_your_workdir")
dir.create("plots/")
dir.create("results/")

### load all required libraries
# If the packages are not available, we will install them using Bioconductor.

if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}
packages <- c("corto", "DESeq2", "EnhancedVolcano", "enrichR", "factoextra", "fgsea", "ggplot2", "magrittr",
            "matrixStats", "msigdbr", "pheatmap", "stringr", "tidyverse", "xlsx"
)
for(p in packages){
  if (!p %in% rownames(installed.packages())){
    BiocManager::install(p)
  }
  library(p, character.only = TRUE)
}

options(java.parameters = "-Xmx8G")
require(corto)
require(DESeq2)
require(EnhancedVolcano)
require(enrichR)
require(factoextra)
require(fgsea)
require(ggplot2)
require(magrittr)
require(matrixStats)
require(msigdbr)
require(pheatmap)
require(stringr)
require(tidyverse)
require(xlsx)
source("textplot3.R")
source("heatmaps.R")

####### DOWNLOAD RAW DATA FROM GSE
## Breast cancer experiments: GSE185760
## DIPG Experiments: GSE215220
## SH-SY5Y Experiments: GSE155031

##### Process FASTQ files in bash! Take a look at the RNAseq_bash_pipeline.sh code

####### Build rawcounts data matrix from read counts --------------

# Import & pre-process ----------------------------------------------------
rawcounts <- read.table("featureCounts.txt", header=TRUE, row.names=1)
# Remove first five columns (chr, start, end, strand, length)
rawcounts <- rawcounts[ ,6:ncol(rawcounts)]
colnames(rawcounts)
# Remove unnecessary characters in colnames
colnames(rawcounts) <- gsub("[^S]*$", "", colnames(rawcounts))
colnames(rawcounts) <- gsub("_S", "", colnames(rawcounts))

# #### Label Samples and Treatments
annotation<-matrix(nrow = ncol(rawcounts),ncol = 2, 
                   dimnames = list(colnames(rawcounts), c("Cell","Treat")))
annotation[,"Cell"] <- ifelse(str_detect(colnames(rawcounts), "DIP") == TRUE,"DIPG","MB231")
annotation[,"Treat"] <- c(rep("CTRL", 3),rep("TEPA", 3),rep("CTRL", 3),rep("TEPA", 3),
                          rep("CTRL", 3),rep("TEPA24", 3),rep("TEPA8", 3))
annotation <- as.data.frame(annotation, stringsAsFactors=FALSE)
save(rawcounts, annotation, file = "results/000_rawcounts.rda")

# Variance Stabilizing Transformation
rawcounts <- as.matrix(rawcounts)
expmat <- vst(rawcounts, blind = TRUE, nsub = 1000, fitType = "parametric")
save(expmat, file = "results/000_expmat.rda")

# PCA
color <- c(rep("green", 6),rep("blue", 6),rep("red", 9))
shapes <- c(rep(1, 3),rep(2, 3),rep(1, 3),rep(2, 3),rep(1, 3),rep(2, 6))
png("plots/000_PCA_Allsamples.png", w = 4000, h = 4000, res = 300)
pca <- prcomp(t(expmat))
totvar <- sum(pca$sdev^2)
pcavar <- ((pca$sdev^2) / totvar)*100
x <- setNames(pca$x[,1], colnames(expmat))
y <- setNames(pca$x[,2], colnames(expmat))
plot(x, y, pch = shapes, main = paste0("TEPA Dataset"),
     xlim = c(min(x)*1.5, max(x)*1.5), cex = 2, col = color,
     xlab = paste0("PC", 1, " (",signif(pcavar[1], 3), "%)"),
     ylab = paste0("PC", 2, " (",signif(pcavar[2], 3), "%)")
)
textplot3(x, y, colnames(expmat), new = FALSE)
grid()
dev.off()

# subset breast
color <- c(rep("green", 3), rep("red", 3), rep("tomato", 3))
breastmat <- expmat[,str_detect(colnames(expmat), "MB")]
png("plots/000_PCA_MB231.png", w = 4000, h = 4000, res = 450)
pca <- prcomp(t(breastmat))
totvar <- sum(pca$sdev^2)
pcavar <- ((pca$sdev^2) / totvar)*100
x <- setNames(pca$x[,1], colnames(breastmat))
y <- setNames(pca$x[,2], colnames(breastmat))
plot(x, y, pch = 20, main = paste0("MB231 Dataset"),
     xlim = c(min(x)*1.5, max(x)*1.5), cex = 2, col = color,
     xlab = paste0("PC", 1, " (", signif(pcavar[1],3), "%)"),
     ylab = paste0("PC", 2, " (", signif(pcavar[2],3), "%)")
)
textplot3(x, y, colnames(breastmat), new = FALSE)
grid()
dev.off()

# DIPG007
dipg007 <- expmat[,str_detect(colnames(expmat), "DIPG007")]
color <- c(rep("green", 3), rep("red", 3))
png("plots/000_PCA_DIPG007.png", w = 4000, h = 4000, res = 450)
pca <- prcomp(t(dipg007))
totvar <- sum(pca$sdev^2)
pcavar <- ((pca$sdev^2) / totvar)*100
x <- setNames(pca$x[,1], colnames(dipg007))
y <- setNames(pca$x[,2], colnames(dipg007))
plot(x, y, pch = 20, main = paste0("DIPG007 Dataset"),
     xlim = c(min(x)*1.5, max(x)*1.5), cex = 2, col = color,
     xlab = paste0("PC", 1, " (", signif(pcavar[1], 3), "%)"),
     ylab=paste0("PC", 2, " (", signif(pcavar[2], 3), "%)")
)
textplot3(x, y, colnames(dipg007), new = FALSE)
grid()
dev.off()

#DIPG010
dipg010 <- expmat[,str_detect(colnames(expmat), "DIPG010")]
color <- c(rep("green", 3), rep("red", 3))
png("plots/000_PCA_DIPG010.png", w = 4000, h = 4000, res = 450)
pca <- prcomp(t(dipg010))
totvar <- sum(pca$sdev^2)
pcavar <- ((pca$sdev^2) / totvar)*100
x <- setNames(pca$x[,1], colnames(dipg010))
y <- setNames(pca$x[,2], colnames(dipg010))
plot(x, y, pch = 20, main = paste0("DIPG010 Dataset"),
     xlim = c(min(x)*1.5, max(x)*1.5), cex = 2,col = color,
     xlab = paste0("PC", 1, " (", signif(pcavar[1], 3), "%)"),
     ylab = paste0("PC", 2, " (", signif(pcavar[2], 3), "%)")
)
textplot3(x, y, colnames(dipg010), new = FALSE)
grid()
dev.off()

rawcounts <- as.matrix(rawcounts)# if not already a matrix!

# Start with BREAST Analysis!
subraw <- rawcounts[,str_detect(colnames(rawcounts), "MB")]
subannot <- annotation[colnames(subraw),]

# Variance Stabilizing Transformation
expmat <- vst(subraw, blind = TRUE, nsub = 1000, fitType = "parametric")
save(expmat, file = "Ema/results/000_expmat.rda")

# Check the complex design to control for time and treat https://support.bioconductor.org/p/62684/
expmat <- expmat[,c(1:3,7:9,4:6)]

# PCA
# subset breast

color <- c(rep("green", 3), rep("red", 3), rep("tomato", 3))
png("plots/000_PCA_MB231.png", w = 3000, h = 3000, res = 450)
pca <- prcomp(t(expmat))
totvar <- sum(pca$sdev^2)
pcavar <- ((pca$sdev^2) / totvar)*100
x <- setNames(pca$x[,1], colnames(expmat))
y <- setNames(pca$x[,2], colnames(expmat))
plot(x, y, pch = 20, main = paste0("MB231 Dataset PCA"),
     xlim = c(min(x)*1.5, max(x)*1.5), cex = 2, col = color,
     xlab = paste0("PC", 1, " (", signif(pcavar[1], 3), "%)"),
     ylab = paste0("PC", 2, " (",signif(pcavar[2], 3), "%)")
)
textplot3(x, y, colnames(expmat), new = FALSE)
grid()
dev.off()

# DESeq2 block (filter out poorly expressed genes)
dds <- DESeqDataSetFromMatrix(countData = subraw, colData = subannot, design = ~ Treat)
dds <- dds[rowSums(counts(dds))>=5,]
dds$Treat <- relevel(dds$Treat, ref = "CTRL")
dea <- DESeq(dds, parallel = TRUE)
resultsNames(dea) #"Treat_TEPA24_vs_CTRL" "Treat_TEPA8_vs_CTRL"
reslist <- list()
reslist[["TEPA8"]] <- as.data.frame(results(dea, name = "Treat_TEPA8_vs_CTRL"))
reslist[["TEPA24"]] <- as.data.frame(results(dea, name = "Treat_TEPA24_vs_CTRL"))
save(reslist, file = "results/001_DiffeExpMB.rda")
write.xlsx(reslist$TEPA8, file = "results/DEA.xlsx", sheetName = "TEPA8h")
write.xlsx(reslist$TEPA24, file = "results/DEA.xlsx", sheetName = "TEPA24h", append = TRUE)

#### DEGs plots -------------
# Volcano Plots

results <- reslist$TEPA24
results <- na.omit(results)
topgenes <- na.omit(results)
topgenes <- topgenes[order(topgenes$padj),]
topgenesup <- rownames(topgenes[topgenes$log2FoldChange>1,])[1:25]
topgenesdn <- rownames(topgenes[topgenes$log2FoldChange<(-1),])[1:25]
top <- c(topgenesup,topgenesdn)

png("plots/001_MB231_TEPA24_vs_Ctrl_volcano.png", w = 2500, h = 2500, res = 300)
EnhancedVolcano(results, subtitle = "",
                lab = rownames(results),
                x = 'log2FoldChange',
                y = 'padj',
                selectLab = top,
                title = 'TEPA 24h vs. CTRL',
                pCutoff = 0.05, #0.05 cutoff
                FCcutoff = 0.5, # 40% fold change
                labFace = "bold",
                labSize = 3,
                col = c('lightgrey', 'pink', 'lightblue', 'salmon'),
                colAlpha = 4/5,
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.3, colConnectors = 'gray51', maxoverlapsConnectors = Inf,
                caption = paste0('Upregulated = ', nrow(results[results$log2FoldChange>0.5&results$padj<=0.05,]), ' genes',"\n",'Downregulated = ',
                                 nrow(results[results$log2FoldChange< -0.5&results$padj<=0.05,]), ' genes'))+
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

# GSEA ------------------
## Load the msigdb database 7.4.1
msigdbr_species()
mdf <- msigdbr(species = "Homo sapiens") # Retrieve all human gene sets
length(mdf$gs_name) # Number of associations 4263110
length(length(unique(mdf$gs_name))) # Number of pathways
head(mdf)
mlist <- mdf %>% split(x=.$gene_symbol, f=.$gs_name)
plength<-sapply(mlist, length)
max(plength) # 2496, the biggest pathway N.B. Numbers may change according to msigdb version!

# Calculate enrichment for TEPA vs CTRL
signature <- setNames(results$stat,rownames(results))
set.seed(1)
gseas <- fgseaMultilevel(pathways = mlist, stats = signature, eps = 0,
                         minSize = 15, maxSize = Inf, nproc = 7, nPermSimple = 10000)
gseas <- gseas[order(gseas$pval),]
write.xlsx2(gseas, file = "results/001_GSEA.xlsx", sheetName = "TEPA24_vs_Ctrl", row.names = FALSE, append = FALSE)
save(gseas, file = "results/001_gsea_MB231_TEPA24.rda")

# keep only sig pathways
gsig <- gseas[gseas$padj<=0.05,]
met <- gsig[grep(("EMT|METAS|MESENCH"), gsig$pathway),] # HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION
tgf <- gsig[grep(("TRANSFORMING|TGF"), gsig$pathway),] # VERRECCHIA_EARLY_RESPONSE_TO_TGFB1
lipid <- gsig[grep(("LIPID"), gsig$pathway),] # REACTOME_METABOLISM_OF_LIPIDS

# GSEA plots with corto -----------------
library(corto) 
paths <- c("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
           "VERRECCHIA_EARLY_RESPONSE_TO_TGFB1",
           "HALLMARK_CHOLESTEROL_HOMEOSTASIS")


# EMT
for (path in paths){
  set1 <- mlist[[path]]
  set.seed(1)
  obj <- gsea(signature, set1, method = 'pareto', np = 1000)
  png(paste0("plots/001_", path, "_24h.png"), w = 1500, h = 1000, res = 300)
  plot_gsea(obj,ext_nes = signif(gseas$NES[gseas$pathway==path], 3),
            ext_pvalue = signif(gseas$padj[gseas$pathway==path],2),
            title=paste0(gsub("_"," ",path), " (24h)"), colBarcode = "#1e90ff",
            twoColors = c("firebrick4", "#1e90ff"))
  dev.off()
}

# 8 hours analysis 
results <- reslist$TEPA8
results <- na.omit(results)
topgenes <- na.omit(results)
topgenes <- topgenes[order(topgenes$padj),]
topgenesup <- rownames(topgenes[topgenes$log2FoldChange>1,])[1:25]
topgenesdn <- rownames(topgenes[topgenes$log2FoldChange<(-1),])[1:25]
top <- c(topgenesup,topgenesdn)

png("plots/001_MB231_TEPA8_vs_Ctrl_volcano.png", w = 2500, h = 2500, res = 300)
EnhancedVolcano(results, subtitle = "",
                lab = rownames(results),
                x = 'log2FoldChange',
                y = 'padj',
                selectLab = top,
                xlim = c(-5, 5),
                title = 'TEPA 8h vs. CTRL',
                pCutoff = 0.05, #0.05 cutoff
                FCcutoff = 0.5, # 40% fold change
                labFace = "bold",
                labSize = 3,
                col = c('lightgrey', 'pink', 'lightblue', 'salmon'),
                colAlpha = 4/5,
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.3, colConnectors = 'gray51', maxoverlapsConnectors = Inf,
                caption = paste0('Upregulated = ', nrow(results[results$log2FoldChange>0.5&results$padj<=0.05,]), ' genes',"\n",'Downregulated = ',
                                 nrow(results[results$log2FoldChange< -0.5&results$padj<=0.05,]), ' genes'))+
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

# Calculate enrichment for TEPA vs CTRL
signature <- setNames(results$stat,rownames(results))
set.seed(1)
gseas <- fgseaMultilevel(pathways = mlist, stats = signature, eps = 0,
                         minSize = 15, maxSize = Inf, nproc = 7, nPermSimple = 10000)
gseas <- gseas[order(gseas$pval),]
write.xlsx2(gseas,file="results/001_GSEA.xlsx", sheetName = "TEPA8_vs_Ctrl", row.names = FALSE, append = TRUE)
save(gseas, file = "results/001_gsea_MB231_TEPA8.rda")

# keep only sig pathways

for (path in paths){
  set1 <- mlist[[path]]
  set.seed(1)
  obj <- gsea(signature, set1, method = 'pareto', np = 1000)
  png(paste0("plots/001_", path, "_8h.png"), w = 1500, h = 1000, res = 300)
  plot_gsea(obj, ext_nes = signif(gseas$NES[gseas$pathway==path], 3),
            ext_pvalue = signif(gseas$padj[gseas$pathway==path], 2),
            title = paste0(gsub("_"," ",path), " (8h)"),
            colBarcode = "#1e90ff", twoColors = c("firebrick4", "#1e90ff"))
  dev.off()
  ledge <- unlist(gseas$leadingEdge[gseas$pathway==path])
  ledge <- unique(ledge)
}

# design heatmaps
load("results/001_gsea_MB231_TEPA8.rda")
gse8 <- gseas
load("results/001_gsea_MB231_TEPA24.rda")
gse24 <- gseas
paletteLength <- 50
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)

#### Plot all 8 hours leading edges
for (path in paths){
  ledge1 <- unlist(gse8$leadingEdge[gse8$pathway==path])
  plotmat <- t(scale(t(expmat[ledge1,c(1:6)])))
  myBreaks <- c(seq(min(plotmat), 0, length.out = ceiling(paletteLength/2) + 1), 
                seq(max(plotmat)/paletteLength, max(plotmat), length.out = floor(paletteLength/2)))
  png(paste0("plots/001_8h_", path, ".png"), w = 1500, h = 3000, res = 300) #h = 3500 if the plot shows bad
  pheatmap(plotmat, cluster_cols = F , cluster_rows = T, color = myColor, breaks = myBreaks) #,border_color="white"
  dev.off()
}

###### Plot all 24 hours leading edge
for (path in paths){
  ledge2 <- unlist(gse24$leadingEdge[gse24$pathway==path])
  plotmat <- t(scale(t(expmat[ledge2,c(1:3,7:9)])))
  myBreaks <- c(seq(min(plotmat), 0, length.out = ceiling(paletteLength/2) + 1), 
                seq(max(plotmat)/paletteLength, max(plotmat), length.out = floor(paletteLength/2)))
  png(paste0("plots/001_24h_", path, ".png"), w = 1500, h = 2200, res = 300) #h = 3500 if the plot shows bad
  pheatmap(plotmat, cluster_cols = F , cluster_rows = T, color = myColor, breaks = myBreaks) #,border_color="white"
  dev.off()
}

### Combined leading edge
for (path in paths){
  ledge1 <- unlist(gse8$leadingEdge[gse8$pathway==path])
  ledge2 <- unlist(gse24$leadingEdge[gse24$pathway==path])
  ledge <- intersect(ledge1,ledge2)
  plotmat <- t(scale(t(expmat[ledge2,])))
  myBreaks <- c(seq(min(plotmat), 0, length.out = ceiling(paletteLength/2) + 1), 
                seq(max(plotmat)/paletteLength, max(plotmat), length.out = floor(paletteLength/2)))
  png(paste0("plots/001_combined_", path, ".png"), w = 1500, h = 2000, res = 300) #h = 3500 if the plot shows bad
  pheatmap(plotmat, cluster_cols = F , cluster_rows = T, color = myColor, breaks = myBreaks) #,border_color="white"
  dev.off()
}

##### Go with the clustering
tepa24 <- reslist$TEPA24
tepa24 <- na.omit(tepa24)
sig_genes24 <- rownames(tepa24[tepa24$padj<=0.05&abs(tepa24$log2FoldChange)>=0.5,])
tepa8 <- reslist$TEPA8
tepa8 <- na.omit(tepa8)
sig_genes8 <- rownames(tepa8[tepa8$padj<=0.05&abs(tepa8$log2FoldChange)>=0.5,])
sig_genes <- unique(c(sig_genes24,sig_genes8))
load("results/000_expmat.rda")
df <- expmat
mat <- matrix(ncol = 3, nrow=nrow(expmat))
rownames(mat) <- rownames(df)
mat[,1] < -rowMeans(df[,1:3]) 
mat[,2] < -rowMeans(df[,7:9]) 
mat[,3] <- rowMeans(df[,4:6]) 
mat <- mat - rowMeans(mat)
mat <- as.data.frame(mat)
mat <- cbind(rownames(mat), mat)
colnames(mat) <- c("geneName", "T0", "T1", "T2")
str(mat)
plotmat <- mat[sig_genes,]
# lineplot https://www.biostars.org/p/343055/
## prepare data for cluster 
for_clust  <- plotmat %>% 
  select(-1) ## remove first column which is gene id 
# Elbow method for number of cluster selection
png("plots/001_elbowplot.png", w = 2000, h = 1500, res = 450)
fviz_nbclust(for_clust, kmeans, method = "wss", k.max = 10) +
  geom_vline(xintercept = 4, linetype = 2)+
  labs(subtitle = "Elbow method")
dev.off()

### kmeans
max_itr <-  10
n_clust  <-  4  ## number of cluster 
set.seed(123) ## reproduce the cluster 
kmeans_out  <- kmeans(for_clust, n_clust, iter.max = max_itr)
## add cluster info to orig matrix 
data_with_cust_info <- plotmat %>% 
  mutate(clust = paste("cluster", kmeans_out$cluster, sep = ""))
## visualise  each cluster
png("plots/001_lineplot.png", w = 2500, h = 2500, res = 450)
data_with_cust_info %>% 
  gather(key = "Time" , value = "StdExp", -c(1,5)) %>%  ### 1 is the index of column 'geneName' and 5 is the index of column 'clust'
  group_by(Time) %>%  
  mutate(row_num =  1:n()) %>% 
  ggplot(aes(x =  Time , y = StdExp , group = row_num,))+
  geom_point() +  
  geom_line(alpha = 0.1 , aes(col = as.character(clust))) + 
  theme_bw() +  
  theme(legend.position = "none" , axis.text.x = element_text(angle = 90 , vjust = 0.4)) +
  facet_wrap(~clust)
dev.off()

clust1 <- rownames(data_with_cust_info[data_with_cust_info$clust=="cluster1",])
clust2 <- rownames(data_with_cust_info[data_with_cust_info$clust=="cluster2",])
clust3 <- rownames(data_with_cust_info[data_with_cust_info$clust=="cluster3",])
clust4 <- rownames(data_with_cust_info[data_with_cust_info$clust=="cluster4",])

## gene enrichment analysis
dbs <- listEnrichrDbs()
dbs <- c("WikiPathways_2019_Human","KEGG_2019_Human","MSigDB_Hallmark_2020","BioPlanet_2019","GO_Biological_Process_2018")
# Cluster 1
go <- enrichr(clust1,dbs)
# plot cluster 1 results
# MSigDB
msigdb <- go$MSigDB_Hallmark_2020
msigdb <- msigdb[1:10,]
msigdb$Term <- paste0(msigdb$Term,ifelse(msigdb$P.value<=0.05, "*", ""))
msigdb <- msigdb[order(msigdb$Combined.Score, decreasing = FALSE),]
msigdb$Term <- factor(msigdb$Term, levels = msigdb$Term)

png("plots/001_MB_cluster1_msigdb.png", w = 3500, h = 2000, res = 450)
par(mar = c(4, 1, 3, 1))
ggplot(msigdb, aes(x = Term, y = Combined.Score, label = Combined.Score)) +
  geom_bar(stat = 'identity', aes(fill = Combined.Score), width = .5, position = 'dodge') +
  labs(subtitle = "Combined scores from MSigDB Hallmark 2020",
       title = "Enriched in Cluster 1")+
  coord_flip() + ylab("Enricher Combined Score") + xlab("") + 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  # get rid of panel grids
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
dev.off()

#### The SY5Y: load processed data from SY5Y data anlysis (a DESEQ2 res object!)
### These data have been previusly analyzed and published here: https://doi.org/10.1158/0008-5472.CAN-20-0471
load("results/003_DESeq2_Chelators_sy5y.rda")
res <- as.data.frame(results(dea, name = "Treatment_TEPA_vs_C2"))
res < -na.omit(res)
sig <- rownames(res[res$log2FoldChange<= 0&res$padj<=0.05,])
gob <- enrichr(sig,dbs)
msigdb <- gob$MSigDB_Hallmark_2020
msigdb <- msigdb[1:11,]
msigdb$Term <- paste0(msigdb$Term, ifelse(msigdb$P.value<=0.05, "*", ""))
msigdb <- msigdb[order(msigdb$Combined.Score, decreasing = FALSE),]
msigdb$Term <- factor(msigdb$Term, levels = msigdb$Term)

png("plots/001_SY5Y_msigdb.png", w = 3500, h = 2000, res = 450)
par(mar = c(4, 1, 3, 1))
ggplot(msigdb, aes(x = Term, y = Combined.Score,label = Combined.Score)) +
  geom_bar(stat = 'identity', aes(fill = Combined.Score), width = .5, position = 'dodge') +
  labs(subtitle = "Combined scores from MSigDB Hallmark 2020",
       title = "Enriched in SY5Y") +
  coord_flip() + ylab("Enricher Combined Score") + xlab("") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  # get rid of panel grids
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
dev.off()
coregenes <- intersect(clust1, sig)
x1 <- res[coregenes,]
x2 <- reslist$TEPA8[coregenes,]
x3 <- reslist$TEPA24[coregenes,]
write.xlsx(x1, file = "results/coregenes.xlsx", sheetName = "SY5Y")
write.xlsx(x2, file = "results/coregenes.xlsx", sheetName = "MB231_8h", append = TRUE)
write.xlsx(x3, file = "results/coregenes.xlsx", sheetName = "MB231_24h", append = TRUE)
##### all downregulated genes
sigdn8 <- reslist$TEPA8
sigdn8 <- na.omit(sigdn8)
sigdn8 <- rownames(sigdn8[sigdn8$log2FoldChange<= -0.5&sigdn8$padj<=0.05,])
sigdn24 <- reslist$TEPA24
sigdn24 <- na.omit(sigdn24)
sigdn24 <- rownames(sigdn24[sigdn24$log2FoldChange<= -0.5&sigdn24$padj<=0.05,])
common <- intersect(sigdn24, sigdn8)
coregenes <- intersect(common, sig)
write.xlsx(coregenes, file = "results/coregenes.xlsx", sheetName = "Globally_down", append = TRUE,
           row.names = FALSE, col.names = FALSE)


### Perform Master Regulator Analysis with corto: DOI 10.1093/bioinformatics/btaa223
# create BRCA TF network: download data from the TCGA consortium!!!
# Tumor tissue counts
tcga_counts <- read.delim("your_TCGA_BRCA_data", sep = "\t", header = TRUE, as.is = TRUE)
traw <- sapply(tcga_counts[,3:ncol(tcga_counts)], as.integer)
samplenames <- gsub("\\.", "-", colnames(traw))
samplenames <- substr(samplenames, 1, 16)
colnames(traw) <- samplenames
rownames(traw) <- tcga_counts$Hugo_Symbol
traw <- traw[, unique(colnames(traw))]
# What is normal, what is cancer
tcgacodes <- substr(colnames(traw), 14, 16)
table(tcgacodes) # 01A 954, 01B 14
tumor <- traw[,tcgacodes=="01A"]
ncol(tumor) #954
rawcounts <- tumor
# VST
## Use DESeq2 function
expmat <- vst(rawcounts, blind = TRUE, nsub = 1000, fitType = "parametric")
save(expmat, file = "results/000_BRCA-expmat.rda")
if (!file.exists("results/TCGA_BRCA_regulon2020.rda")){
  # Create BRCA network
  load("F:/Datasets/Datasets/Liste/tfs_2020.rda")
  # the list can be downloaded from the corto paper referenced above
  # create a corto regulon with GTEx Hippocampus dataset
  regulon <- corto(expmat, centroids = centroids, nbootstraps = 1000, p = 1e-10, nthreads = 7)
  save(regulon, file = "results/TCGA_BRCA_regulon2020.rda")
} else {load("results/TCGA_BRCA_regulon2020.rda")}

# load Orazio TEPA expmat
load("MB231_EXPMAT!")
load("results/000_rawcounts.rda")
submat <- expmat[,str_detect(colnames(expmat), "MB")]
subannot <- annotation[str_detect(colnames(expmat), "MB"),]
h8 <- submat[,c(1:3,7:9)]
h8 <- h8[rowVars(h8)>0,]
#### TEPA 8h vs. Ctrl different contrasts
trt <- h8[,4:6]
ctr <- h8[,1:3]
mra <- corto::mra(trt, ctr, regulon = regulon, minsize = 15, nperm = 1000, nthreads = 7)
save(mra, file = "results/002_MRA_TEPA8_MB231.rda")
# plot 8h
dn <- names(sort(mra$nes))[1:10]
up <- names(sort(mra$nes, decreasing = TRUE))[1:10]
png("plots/001_MRA_DN.png", w = 3000, h = 4500, res = 300)
mraplot(mra, mrs = dn, pthr = 0.05, title = "TEPA 8h vs. Ctrl")
dev.off()
png("plots/001_MRA_UP.png",w = 3000, h = 4500, res = 300)
mraplot(mra,mrs = up, pthr = 0.05, title = "Tepa 8h vs. Ctrl")
dev.off()

#### TEPA 24h vs. Ctrl different contrasts
h24 <- submat[,1:6]
h24 <- h24[rowVars(h24)>0,]
### prepare the contrasts
trt <- h24[,4:6]
ctr <- h24[,1:3]
mra <- corto::mra(trt, ctr, regulon = regulon, minsize = 15, nperm = 1000, nthreads = 7)
save(mra, file = "results/002_MRA_TEPA24_MB231.rda")
# plot 24h
dn <- names(sort(mra$nes))[1:10]
up <- names(sort(mra$nes, decreasing = TRUE))[1:10]
png("plots/001_MRA24_DN.png", w = 3000, h = 4500, res = 300)
mraplot(mra, mrs = dn, pthr = 0.05, title = "TEPA 24h vs. Ctrl")
dev.off()
png("plots/001_MRA24_UP.png", w = 3000, h = 4500, res = 300)
mraplot(mra, mrs = up, pthr = 0.05, title = "Tepa 24h vs. Ctrl")
dev.off()

load("results/002_MRA_TEPA8_MB231.rda")
png("plots/001_MRA8_SNAI2.png", w = 2500, h = 1000, res = 300)
mraplot(mra, mrs = "SNAI2", pthr = 0.05, title = "Tepa 8h vs. Ctrl")
dev.off()

load("results/002_MRA_TEPA24_MB231.rda")
png("plots/001_MRA24_SNAI2.png", w = 2500, h = 1000, res = 300)
mraplot(mra, mrs = "SNAI2", pthr = 0.05, title = "Tepa 24h vs. Ctrl")
dev.off()

###### This code part is for the analysis of DIPG
### Start with DIPG010
load("data/000_rawcounts.rda")
rawcounts <- as.matrix(rawcounts)# if not already a matrix!
subraw <- rawcounts[,str_detect(colnames(rawcounts), "DIPG010")]
subannot <- annotation[colnames(subraw),]

# DESeq2 block (filter out poorly expressed genes)
dds <- DESeqDataSetFromMatrix(countData = subraw, colData = subannot, design = ~ Treat)
dds <- dds[rowSums(counts(dds))>=5,]
dds$Treat <- relevel(dds$Treat, ref = "CTRL")
dea <- DESeq(dds, parallel = TRUE)
resultsNames(dea) #"Treat_TEPA_vs_CTRL" 
res <- as.data.frame(results(dea, name = "Treat_TEPA_vs_CTRL"))
results <- res[order(res$padj, decreasing = F),]
results <- na.omit(results)
write.xlsx2(results, file = "results/000_DiffExp.xlsx", sheetName = "DIPG010_TEPA24")
save(dea, file = "results/000_DE_DIPG010.rda")

#### DEGs plots

# Volcano Plots
# Volcano Plots
topgenes <- na.omit(results)
topgenes <- topgenes[order(topgenes$padj),]
topgenesup <- rownames(topgenes[topgenes$log2FoldChange>1,])[1:25]
topgenesdn <- rownames(topgenes[topgenes$log2FoldChange<(-1),])[1:25]
top <- c(topgenesup,topgenesdn)

#results<-na.omit(results)
png("plots/000_DIPG010_TEPA24_vs_Ctrl_volcano.png", w = 2500, h = 2500, res = 300)
EnhancedVolcano(results, subtitle = "",
                lab = rownames(results),
                x = 'log2FoldChange',
                y = 'padj',
                selectLab = top,
                title = 'TEPA 24h vs. CTRL',
                pCutoff = 0.05, #0.05 cutoff
                FCcutoff = 1, # 2-fold change
                labFace = "bold",
                labSize = 3,
                col = c('lightgrey', 'pink', 'lightblue', 'salmon'),
                colAlpha = 4/5,
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.3, colConnectors = 'gray51', maxoverlapsConnectors = Inf,
                caption = paste0('Upregulated = ', nrow(results[results$log2FoldChange>1&results$padj<=0.05,]), ' genes',"\n",'Downregulated = ',
                                 nrow(results[results$log2FoldChange< -1&results$padj<=0.05,]), ' genes'))+
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

# # Gene Ontology
dbs <- listEnrichrDbs()
dbs <- c("WikiPathways_2019_Human", "MSigDB_Hallmark_2020")

de <- results
up <- rownames(de[de$log2FoldChange>1&de$padj<=0.05,])
dn <- rownames(de[de$log2FoldChange< -1&de$padj<=0.05,])

eup <- enrichr(up, dbs)
edn <- enrichr(dn, dbs)

up <- eup$WikiPathways_2019_Human
down <- edn$WikiPathways_2019_Human

up$type <- "up"
down$type <- "down"

up <- up[c(1:10),]
up <- up[order(up$Combined.Score),]
down <- down[c(1:10),]
down$Combined.Score <- (-1)*down$Combined.Score
down <- down[order(down$Combined.Score),]
gos <- rbind(down,up)
gos$Term <- gsub("WP(.*)", "", gos$Term)
gos$Term <- paste0(gos$Term, ifelse(gos$P.value<=0.05, "*", ""))
gos$Term <- factor(gos$Term, levels = gos$Term)

#Diverging Barchart
png("plots/000_GO_bar_TEPA24_Wiki_DIPG007.png", w = 2500, h = 1500, res = 300)
ggplot(gos, aes(x = Term, y = Combined.Score, label = Combined.Score)) +
  geom_bar(stat = 'identity', aes(fill = type), width = .5, position = 'dodge') +
  scale_fill_manual(name = "Expression",
                    labels = c("Down regulated", "Up regulated"),
                    values = c("down" = "lightblue","up" = "#f8766d")) +
  labs(subtitle = "Combined scores from Wiki Pathways",
       title = "Enriched in 24h TEPA-treated DIPG007 cells") +
  coord_flip() +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

# # Hallmark
up <- eup$MSigDB_Hallmark_2020
down <- edn$MSigDB_Hallmark_2020

up$type <- "up"
down$type <- "down"

up <- up[c(1:5),]
up <- up[order(up$Combined.Score),]
down <- down[c(1:5),]
down$Combined.Score <- (-1)*down$Combined.Score
down <- down[order(down$Combined.Score),]
gos <- rbind(down,up)
#gos$Term<-gsub("WP(.*)","",gos$Term)
gos$Term <- paste0(gos$Term, ifelse(gos$P.value<=0.05, "*", ""))
gos$Term <- factor(gos$Term, levels = gos$Term)

#Diverging Barchart
png("plots/000_GO_bar_TEPA24_MSigDB_DIPG010.png", w = 2500, h = 1500, res = 300)
ggplot(gos, aes(x = Term, y = Combined.Score, label = Combined.Score)) +
  geom_bar(stat = 'identity', aes(fill = type), width = .5, position = 'dodge') +
  scale_fill_manual(name = "Expression",
                    labels = c("Down regulated", "Up regulated"),
                    values = c("down" = "lightblue", "up" = "#f8766d"))+
  labs(subtitle = "Combined scores from MSigDB Pathways",
       title = "Enriched in 24h TEPA-treated DIPG007 cells") +
  coord_flip() +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

# GSEA ------------------
msigdbr_species()
mdf <- msigdbr(species = "Homo sapiens") # Retrieve all human gene sets
length(mdf$gs_name) # Number of associations 4263110
length(length(unique(mdf$gs_name))) # Number of pathways
head(mdf)
mlist <- mdf %>% split(x=.$gene_symbol,f=.$gs_name)
#grep("(STAT1|STAT2|STAT3)",names(mlist),value="TRUE")
plength <- sapply(mlist,length)
max(plength) # 2496, the biggest pathway

# Calculate enrichment for TEPA vs CTRL
signature <- setNames(results$stat,rownames(results))
set.seed(1)
gseas <- fgseaMultilevel(pathways = mlist,stats = signature, eps = 0, minSize = 15, maxSize = Inf,
                         nproc=14, nPermSimple = 10000)
gseas <- gseas[order(gseas$pval),]
write.xlsx2(gseas, file = "results/000_GSEA.xlsx", sheetName = "TEPA24_vs_CtrlDIPG010", row.names = FALSE,
            append = FALSE)
save(gseas, file = "results/001_gsea_DIPG010_TEPA24.rda")

# keep only sig pathways
gsig <- gseas[gseas$padj<=0.05,]
met <- gsig[grep(("EMT|METAS|MESENCH"),gsig$pathway),] # HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION
tgf <- gsig[grep(("TRANSFORMING|TGF"),gsig$pathway),] # VERRECCHIA_EARLY_RESPONSE_TO_TGFB1
#lipid<-gsig[grep(("LIPID"),gsig$pathway),] #REACTOME_METABOLISM_OF_LIPIDS

# GSEA plots with corto -----------------
paths <- c("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "VERRECCHIA_EARLY_RESPONSE_TO_TGFB1")

# EMT
for (path in paths){
  set1 <- mlist[[path]]
  set.seed(1)
  obj <- gsea(signature, set1, method = 'pareto', np = 1000)
  png(paste0("plots/001_DIPG010_",path,"_24h.png"), w = 1500, h = 1000, res = 300)
  plot_gsea(obj, ext_nes = signif(gseas$NES[gseas$pathway==path], 3),
            ext_pvalue = signif(gseas$padj[gseas$pathway==path], 2),
            title = paste0(gsub("_", " ", path), " (24h)"),
            colBarcode = "#1e90ff", twoColors = c("firebrick4", "#1e90ff"))
  dev.off()
}

# heatmap
paletteLength <- 50
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
###### Plot all 24 hours leading edge
for (path in paths){
  ledge2 <- unlist(gseas$leadingEdge[gseas$pathway==path])
  plotmat <- t(scale(t(expmat[ledge2, c(7:12)])))
  myBreaks <- c(seq(min(plotmat), 0, length.out = ceiling(paletteLength/2) + 1), 
                seq(max(plotmat)/paletteLength, max(plotmat), length.out = floor(paletteLength/2)))
  png(paste0("plots/001_24h_DIPG010_heatmaps_",path,".png"), w = 1500, h = 3500, res = 300) #3500 h if preferred
  pheatmap(plotmat, cluster_cols = F ,cluster_rows = T, color = myColor, breaks = myBreaks) #,border_color="white"
  dev.off()
}
# Diverging barplots
# All pathways
gsea <- as.data.frame(gseas)
gsea <- gsea[order(-abs(gsea$NES)),]

### Table of top pathways ----
top <- gsea[gsea$NES<0,][1:15,]
top <- rbind(top,gsea[gsea$NES>0,][15:1,])
top <- top[order(top$NES),]
toplot <- setNames(top$NES, top$pathway)
# Format
#names(toplot)<-gsub("GO_","",names(toplot))
names(toplot) <- gsub("_", " ", names(toplot))
names(toplot) <- str_to_title(names(toplot))
png("plots/004_gsea_DIPG007_TEPA24_top.png", w = 6000, h = 3000, res = 500)
par(mar = c(4, 1, 3, 1))
bp <- barplot(toplot, horiz = TRUE, xlab = "Normalized Enrichment Score",
            xlim = 1.3*c(-max(abs(toplot)), max(abs(toplot))),
            main = "TEPA 24 vs. CTRL, top Pathways",
            col = rep(c("cornflowerblue", "salmon"), each = 15),
            yaxt = "n", cex.main = 2
)
text(0,bp[1:15,1], names(toplot)[1:15], pos = 4)
text(0,bp[16:30,1], names(toplot)[16:30], pos = 2)
abline(v = c(-p2z(0.05), p2z(0.05)), lty = 2)
dev.off()

# plot Sarrio leading edge

ledge <- gsea[gsea$pathway=="SARRIO_EPITHELIAL_MESENCHYMAL_TRANSITION_UP"|
              gsea$pathway=="SARRIO_EPITHELIAL_MESENCHYMAL_TRANSITION_DN",]
genelist <- unlist(ledge$leadingEdge)
upg <- unlist(ledge[1,8])
dpg <- unlist(ledge[2,8])
colside <- c(rep("red", length(upg)), rep("blue", length(dpg)))
names(colside) <- c(upg, dpg)
mat_zscore <- t(scale(t(expmat[,rownames(subannot)])))
mat_zscore <- mat_zscore[genelist,]
png("plots/004_heatmap_sarrio.png", w = 2500, h = 6000, res = 300)
heatmap.3(mat_zscore, KeyValueName = "VST-normalized exprs", RowSideColors = colside)
dev.off()

# heatmap
cols <- as.data.frame(colside)
colnames(cols) <- "EMT"
cols[cols$EMT=="red",] <- "up"
cols[cols$EMT=="blue",] <- "dn"
cols$EMT <- as.factor(cols$EMT)

# Specify colors
ann_colors = list(
  EMT = c(up = "red", dn ="blue"))

paletteLength <- 50
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
mat_zscore<-t(scale(t(expmat[,rownames(subannot)])))
mat_zscore<-mat_zscore[genelist,]
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(min(mat_zscore), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(mat_zscore)/paletteLength, max(mat_zscore), length.out=floor(paletteLength/2)))
png("plots/001_EMT_heatmap.png", w = 1500, h = 5500, res = 300)
pheatmap(mat_zscore, cluster_cols = , cluster_rows = T, color = myColor, breaks = myBreaks,
         annotation_row = cols, 
         annotation_colors = ann_colors)
dev.off()

# sarrio 2 way

# the SARRIO 2way GSEA
set1 <- mlist[["SARRIO_EPITHELIAL_MESENCHYMAL_TRANSITION_UP"]]
set2 <- mlist[["SARRIO_EPITHELIAL_MESENCHYMAL_TRANSITION_DN"]]
set.seed(1)
obj <- gsea2(signature, set2, set1, method = 'pareto', np = 1000)

png("plots/004_007_2WGSEA_SARRIO.png", w = 2500, h = 2500, res = 450)
plot_gsea2(obj, twoColors = c("red3", "dodgerblue3"), bottomTitle = "SARRIO_EPITHELIAL_MESENCHYMAL_TRANSITION_UP",
           title = "SARRIO_EPITHELIAL_MESENCHYMAL_TRANSITION_DN")
dev.off()

# Epigenetic
ledge <- gseas[gseas$pathway=="REACTOME_EPIGENETIC_REGULATION_OF_GENE_EXPRESSION",]
genelist <- unlist(ledge$leadingEdge)
exp <- expmat
mat_zscore <- t(scale(t(expmat[,rownames(subannot)])))
mat_zscore <- mat_zscore[genelist,]
png("004_heatmap_Epigenetic.png", w = 2500, h = 4000, res = 450)
heatmap.3(mat_zscore, KeyValueName = "VST-norm exprs")
dev.off()

##### Analyze DIPG 007
subraw<-rawcounts[,str_detect(colnames(rawcounts),"DIPG007")]
subannot<-annotation[colnames(subraw),]

# DESeq2 block (filter out poorly expressed genes)
dds <- DESeqDataSetFromMatrix(countData = subraw, colData = subannot, design = ~ Treat)
dds <- dds[rowSums(counts(dds))>=5,]
dds$Treat <- relevel(dds$Treat, ref = "CTRL")
dea <- DESeq(dds, parallel = TRUE)
resultsNames(dea) #"Treat_TEPA_vs_CTRL" 
res <- as.data.frame(results(dea, name = "Treat_TEPA_vs_CTRL"))
results <- res[order(res$padj,decreasing = F),]
results <- na.omit(results)
write.xlsx2(results, file = "results/000_DiffExp.xlsx", sheetName = "DIPG007_TEPA24")
save(dea,file = "results/000_DE_DIPG007.rda")

#### DEGs plots
# Volcano Plots

topgenes <- na.omit(results)
topgenes <- topgenes[order(topgenes$padj),]
topgenesup <- rownames(topgenes[topgenes$log2FoldChange>1,])[1:25]
topgenesdn <- rownames(topgenes[topgenes$log2FoldChange<(-1),])[1:25]
top <- c(topgenesup, topgenesdn)

png("plots/000_DIPG007_TEPA24_vs_Ctrl_volcano.png", w = 2500, h = 2500, res = 300)
EnhancedVolcano(results, subtitle = "",
                lab = rownames(results),
                x = 'log2FoldChange',
                y = 'padj',
                selectLab = top,
                title = 'TEPA 24h vs. CTRL',
                pCutoff = 0.05, #0.05 cutoff
                FCcutoff = 1, # 2-fold change
                labFace = "bold",
                labSize = 3,
                col = c('lightgrey', 'pink', 'lightblue', 'salmon'),
                colAlpha = 4/5,
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.3,colConnectors = 'gray51',maxoverlapsConnectors = Inf,
                caption = paste0('Upregulated = ', nrow(results[results$log2FoldChange>1&results$padj<=0.05,]), ' genes',"\n",'Downregulated = ',
                                 nrow(results[results$log2FoldChange< -1&results$padj<=0.05,]), ' genes'))+
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

# # Gene Ontology
dbs <- listEnrichrDbs()
dbs <- c("WikiPathways_2019_Human", "MSigDB_Hallmark_2020")

de <- results
up <- rownames(de[de$log2FoldChange>1&de$padj<=0.05,])
dn <- rownames(de[de$log2FoldChange< -1&de$padj<=0.05,])

eup <- enrichr(up,dbs)
edn <- enrichr(dn,dbs)

up <- eup$WikiPathways_2019_Human
down <- edn$WikiPathways_2019_Human

up$type <- "up"
down$type <- "down"

up <- up[c(1:10),]
up <- up[order(up$Combined.Score),]
down <- down[c(1:10),]
down$Combined.Score <- (-1)*down$Combined.Score
down <- down[order(down$Combined.Score),]
gos <- rbind(down,up)
gos$Term <- gsub("WP(.*)", "", gos$Term)
gos$Term <- paste0(gos$Term, ifelse(gos$P.value<=0.05, "*", ""))
gos$Term <- factor(gos$Term, levels = gos$Term)

#Diverging Barchart
png("plots/000_GO_bar_TEPA24_Wiki_DIPG007.png", w = 2500, h = 1500, res = 300)
ggplot(gos,aes(x = Term, y = Combined.Score, label = Combined.Score)) +
  geom_bar(stat = 'identity', aes(fill = type), width = .5, position = 'dodge') +
  scale_fill_manual(name = "Expression",
                    labels = c("Down regulated", "Up regulated"),
                    values = c("down" = "lightblue", "up" = "#f8766d")) +
  labs(subtitle = "Combined scores from Wiki Pathways",
       title = "Enriched in 24h TEPA-treated DIPG007 cells") +
  coord_flip() +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

# GSEA ------------------

# Calculate enrichment for TEPA vs CTRL
signature <- setNames(results$stat, rownames(results))
set.seed(1)
gseas <- fgseaMultilevel(pathways = mlist, stats = signature, eps = 0,
                         minSize = 15, maxSize = Inf, nproc = 14, nPermSimple = 10000)
gseas <- gseas[order(gseas$pval),]
write.xlsx2(gseas, file = "results/000_GSEA.xlsx", sheetName = "TEPA24_vs_Ctrl", row.names = FALSE, append = FALSE)
save(gseas, file = "results/001_gsea_MB231_TEPA24.rda")

# keep only sig pathways
gsig <- gseas[gseas$padj<=0.05,]
met <- gsig[grep(("EMT|METAS|MESENCH"), gsig$pathway),] # SARRIO_EPITHELIAL_MESENCHYMAL_TRANSITION_UP
tgf <- gsig[grep(("TRANSFORMING|TGF"), gsig$pathway),] # KARLSSON_TGFB1_TARGETS_UP

# Diverging barplots
# All pathways
gsea <- as.data.frame(gseas)
gsea <- gsea[order(-abs(gsea$NES)),]

### Table of top pathways ----
top <- gsea[gsea$NES<0,][1:15,]
top <- rbind(top,gsea[gsea$NES>0,][15:1,])
top <- top[order(top$NES),]
toplot <- setNames(top$NES, top$pathway)
# Format
names(toplot) <- gsub("_", " ", names(toplot))
names(toplot) <- str_to_title(names(toplot))
png("plots/004_gsea_DIPG007_TEPA24_top.png", w = 6000, h = 3000, res = 500)
par(mar = c(4, 1, 3, 1))
bp <- barplot(toplot, horiz = TRUE, xlab = "Normalized Enrichment Score",
            xlim = 1.3*c(-max(abs(toplot)), max(abs(toplot))),
            main = "TEPA 24 vs. CTRL, top Pathways",
            col = rep(c("cornflowerblue", "salmon"), each = 15),
            yaxt = "n", cex.main = 2
)
text(0, bp[1:15,1], names(toplot)[1:15], pos=4)
text(0, bp[16:30,1], names(toplot)[16:30], pos=2)
abline(v = c(-p2z(0.05), p2z(0.05)), lty=2)
dev.off()

# plot Sarrio leading edge
ledge <- gsea[gsea$pathway=="SARRIO_EPITHELIAL_MESENCHYMAL_TRANSITION_UP"|
              gsea$pathway=="SARRIO_EPITHELIAL_MESENCHYMAL_TRANSITION_DN",]
genelist <- unlist(ledge$leadingEdge)
upg <- unlist(ledge[1,8])
dpg <- unlist(ledge[2,8])
colside <- c(rep("red", length(upg)), rep("blue", length(dpg)))
names(colside) <- c(upg, dpg)
mat_zscore <- t(scale(t(expmat[,rownames(subannot)])))
mat_zscore <- mat_zscore[genelist,]
png("plots/004_heatmap_sarrio.png", w = 2500, h = 6000, res = 300)
heatmap.3(mat_zscore, KeyValueName = "VST-normalized exprs", RowSideColors = colside)
dev.off()

# heatmap
cols <- as.data.frame(colside)
colnames(cols) <- "EMT"
cols[cols$EMT=="red",] <- "up"
cols[cols$EMT=="blue",] <- "dn"
cols$EMT <- as.factor(cols$EMT)

# Specify colors
ann_colors = list(
  EMT = c(up = "red", dn ="blue"))

paletteLength <- 50
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
#load("MB231/results/000_expmat.rda")
mat_zscore<-t(scale(t(expmat[,rownames(subannot)])))
#ledge<-c(ledge,ledge8)
#ledge<-unique(ledge)
mat_zscore<-mat_zscore[genelist,]
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(min(mat_zscore), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(mat_zscore)/paletteLength, max(mat_zscore), length.out=floor(paletteLength/2)))
png("plots/001_EMT_heatmap.png", w = 1500, h = 5500, res = 300)
pheatmap(mat_zscore,cluster_cols = ,cluster_rows = T, color = myColor, breaks = myBreaks,
         annotation_row = cols, 
         annotation_colors = ann_colors)
dev.off()

# sarrio 2 way
# the SARRIO 2way GSEA
set1 <- mlist[["SARRIO_EPITHELIAL_MESENCHYMAL_TRANSITION_UP"]]
set2 <- mlist[["SARRIO_EPITHELIAL_MESENCHYMAL_TRANSITION_DN"]]
set.seed(1)
obj <- gsea2(signature, set2, set1, method = 'pareto', np = 1000)

png("plots/004_007_2WGSEA_SARRIO.png", w = 2500, h = 2500, res = 450)
plot_gsea2(obj, twoColors = c("red3", "dodgerblue3"),
           bottomTitle = "SARRIO_EPITHELIAL_MESENCHYMAL_TRANSITION_UP",
           title = "SARRIO_EPITHELIAL_MESENCHYMAL_TRANSITION_DN")
dev.off()

set <- mlist[["REACTOME_EPIGENETIC_REGULATION_OF_GENE_EXPRESSION"]]
load("results/004_DE_DIPG007.rda")
resultsNames(dea)
results <- as.data.frame(results(dea, name = "Treat_TEPA_vs_CTRL"))
signature <- setNames(results$stat, rownames(results))
set.seed(1)
obj <- gsea(signature, set, method = 'pareto', np = 1000)
png("DIPG/plots/004_EPIGENETIC.png", w = 2000, h = 2000, res = 300)
plot_gsea(obj, twoColors = c("red", "deepskyblue3"),
          colBarcode = "deepskyblue3", title = "REACTOME EPIGENETIC REGULATION OF GENE EXPRESSION")
dev.off()