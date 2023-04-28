#-------------------------------------------------------------------------------
# Thermofly: RNA-seq library QC analysis  (of D. mel. acute thermal shock)
# Comparing polyA enhancement to rRNA depletion for SALMON transcript data
# Sara Filler & Nivea Patel
# 7/25/2022
#
# based off of DESeq2_salmon.Rmd
#-------------------------------------------------------------------------------

# Load libraries ----
library(tximport)
library(GenomicFeatures)
library(readr)
library(gprofiler2)
library(DESeq2)
library(tidyverse)
library(hciR)
library(pheatmap)
library(vidger)
library(ggVennDiagram)
library(UpSetR)
library(ggdendro)
library(geneplotter)
library(kableExtra)
library(openxlsx)

# RNA-seq QC all together ----

txdb <- makeTxDbFromGFF("~/../joeboyd/indexes/DM6/GTF/dm6.ensGene.gtf")
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")
head(tx2gene)

dir = "/slipstream_old/home/sethandy/sethandy/thermofly/salmon"

sample <- list.files(dir) %>% str_replace_all( ".salmon_quant", "")
#group <- sample %>% str_replace_all( "_R.", "")
library = sapply(strsplit(sample, '_'), `[`, 1)
condition = sapply(strsplit(sample, '_'), `[`, 2)
condition<- factor(condition,levels=c("25","37","34","10","4"))
replicate = sapply(strsplit(sample, '_'), `[`, 3)
directory = list.files(dir)
#make df
annotation <- data.frame(sample = sample, directory= directory, library = library, condition = condition, replicate = replicate, batch = c("A"), stringsAsFactors = TRUE)

# mark rows for batch effect
batch_rows<- (annotation["sample"]=="polyA_10_2")|
  (annotation["sample"]=="polyA_34_2")|
  (annotation["sample"]=="polyA_34_3")|
  (annotation["sample"]=="polyA_37_1")|
  (annotation["sample"]=="polyA_37_2")
annotation$batch<-as.character(annotation$batch)
annotation[batch_rows,"batch"] = "B"
annotation$batch<-as.factor(annotation$batch)

# continue making df
files <- file.path(dir,annotation$directory, "quant.sf")
names(files) <- paste0(annotation$sample)
all(file.exists(files))
txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene)
salmon.df = as.data.frame(txi.salmon$counts) %>% dplyr::select(-rRNA_10_4)
salmon.mat = as.matrix(salmon.df)
# getting rid of outlier extra sample
anno = annotation[annotation$sample != "rRNA_10_4",]
rownames(anno) = colnames(salmon.mat)

#setup DESeq2 design and run
dds_all <- DESeqDataSetFromMatrix(round(salmon.mat), 
                                colData = anno, 
                                design = ~batch + condition)

dds_all <- DESeq(dds_all)

# Save as RDS object so dds object is easily loadable later
# wd = "/slipstream_old/home/sethandy/sethandy/thermofly"
# setwd(wd)
# saveRDS(dds_all,"all_salmon_dds.RDS")
# # readRDS("/slipstream_old/home/sethandy/sethandy/thermofly/all_salmon_dds.rds")

vst_all <- varianceStabilizingTransformation(dds_all, blind=TRUE)
vst_mat_all <- assay(vst_all)


## library sizes and gene count densities ----

#barplot of library sizes
library(RColorBrewer)
librarySizes <- colSums(salmon.df)/10^6
allSample.cols <- rep(brewer.pal(2, "Set1"), each=15)
### Change the margins so the sample names will fit
mymar <- par()$mar
mymar[2] <- 7
mymar[1] <- 4
par(mar = mymar)

# default
par(mar=c(5.1, 4.1, 4.1, 2.1))

par(mar=c(5.1, 3.1, 3.1, 2.1))

par(mar=c(7, 4.1, 4.1, 2.1))

###Use the barplot function
barplot(librarySizes, 
        col = allSample.cols, 
        xlab="Reads per million", 
        main="Reads mapped to genes", 
        horiz=TRUE,
        cex.names = .35,
        #legend.text=c("rRNA","polyA"),
        ylab="library",
        las=1)
legend("topright", inset = 0.0125,legend= c("rRNA","polyA")
       ,fill = c("steelblue","red"), horiz=TRUE)

legend(fill=c("blue","red"))


# original plot
barplot(librarySizes, 
        names=names(librarySizes), 
        las=2, 
        main="Reads mapped to genes")

# Normalized Library sizes ----
# typical expression level of a gene after it's been normalized
# (just want to be relatively even)

logcounts <- log2(vst_mat_all + 1)
# make a color vector
statusCol <- as.numeric(factor(annotation$condition)) + 1
# Check distributions of samples using boxplots
boxplot(logcounts, 
        xlab="", 
        ylab="Log2(counts)",
        las=2,
        col=statusCol)

#number of genes that have non–zero counts in all samples
GeneCounts <- counts(dds_all)
idx.nz <- apply(GeneCounts, 1, function(x) { all(x > 0)})
sum(idx.nz)

#estimate the size factors
DESeq2Table <- estimateSizeFactors(dds_all)
sizeFactors(DESeq2Table)
#densities of counts
multidensity( counts(DESeq2Table, normalized = T)[idx.nz ,],
              xlab="mean counts", xlim=c(0, 1000), legend = "center")

#empircal cumulative distribution functions (ECDFs)
# essentially integrals of the densities and 
# give the probability of observing a certain number of counts equal 
# to x or less given the data.
multiecdf( counts(DESeq2Table, normalized = T)[idx.nz ,],
           ylab="estimated probability(Fn(x))", xlab="mean counts", 
           xlim=c(0, 1000), 
           legend = ("center")) # drops the legend

## sample correlation ----

d <- cor((vst_mat_all), method="spearman")
hc <- hclust(dist(1-d))
ggdendrogram(hc, rotate = TRUE, theme_dendro = FALSE)

anno2 = anno
anno2$directory = NULL
anno2$replicate = NULL
anno2$batch = NULL
anno2$sample = NULL
vst_cor <- cor(vst_mat_all)
pheatmap(vst_cor, annotation_row = anno2, show_rownames = F, show_colnames = F)


## PCA ----

# Comparing libraries
# Principal components analysis
# examine how batch variability contributes to a PCA plot

# pre-batch removal, color batch
plotPCA(vst_all, intgroup="library")

#pre batch PCA all together
pcaData <- plotPCA(vst_all, intgroup=c("condition", "library"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition, shape=library)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

# #remove batch effect then replot PCA
# assay(vst_all) <- limma::removeBatchEffect(assay(vst_all), vst_all$batch)
# 
# # post-batch removal, color condition, shape replicate
# pcaData <- plotPCA(vst_all, intgroup=c("condition", "batch"), returnData=TRUE)
# percentVar <- round(100 * attr(pcaData, "percentVar"))
# ggplot(pcaData, aes(PC1, PC2, color=condition, shape=batch)) +
#   geom_point(size=3) +
#   xlab(paste0("PC1: ",percentVar[1],"% variance")) +
#   ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
#   coord_fixed()


# RNA-seq QC rRNA-depletion only ----

df.salmon_A = as.data.frame(txi.salmon$counts) 
df.A = df.salmon_A %>% dplyr::select(contains("rRNA")) %>% dplyr::select(-rRNA_10_4)
df.A <- cbind(rownames(df.A), data.frame(df.A, row.names=NULL))
names(df.A)[1] <- "GeneID"
mat.A <- as_tibble(df.A) %>% as_matrix()

#annotation
file = colnames(mat.A)
condition = sapply(strsplit(colnames(mat.A), '_'), `[`, 2) 
replicate = sapply(strsplit(colnames(mat.A), '_'), `[`, 3) 
#make df
annotation <- data.frame(file = file, condition = condition, replicate = replicate, stringsAsFactors = TRUE)
row.names(annotation$file) 
anno <- annotation[,-1]
rownames(anno) <- annotation[,1]

dds.A <- DESeqDataSetFromMatrix(countData  = round(mat.A), 
                                  colData = annotation, 
                                  design = ~condition)


dds <- DESeq(dds.A)

vst <- varianceStabilizingTransformation(dds, blind=TRUE)
vst_mat <- assay(vst)


## rRNA-depletion correlation ----

d <- cor((vst_mat), method="spearman")
hc <- hclust(dist(1-d))
ggdendrogram(hc, rotate = TRUE, theme_dendro = FALSE)

#PCA

pcaData <- plotPCA(vst, intgroup=c("condition", "replicate"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition, shape=replicate)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()



# RNA-seq QC polyA-enhancment only ----

df.salmon_A = as.data.frame(txi.salmon$counts) 
df.A = df.salmon_A %>% dplyr::select(contains("polyA"))
df.A <- cbind(rownames(df.A), data.frame(df.A, row.names=NULL))
names(df.A)[1] <- "GeneID"
mat.A <- as_tibble(df.A) %>% as_matrix()
# rename polyA_10_4 to polyA_10_3 for plot
colnames(mat.A) =  c("polyA_10_1", "polyA_10_2", "polyA_10_3", "polyA_25_1", "polyA_25_2", "polyA_25_3", "polyA_34_1", "polyA_34_2", "polyA_34_3", "polyA_37_1", "polyA_37_2", "polyA_37_3", "polyA_4_1", "polyA_4_2", "polyA_4_3")

#annotation
file = colnames(mat.A)
condition = sapply(strsplit(colnames(mat.A), '_'), `[`, 2) 
replicate = sapply(strsplit(colnames(mat.A), '_'), `[`, 3) 
#batch = c("A", "B", rep("A",5), rep("B",4), rep("A",4)) # not reliable

#make df
annotation <- data.frame(file = file, condition = condition, replicate = replicate, batch = c("A"), stringsAsFactors = TRUE)

# mark rows for batch effect
batch_rows<- (annotation["file"]=="polyA_10_2")|
  (annotation["file"]=="polyA_34_2")|
  (annotation["file"]=="polyA_34_3")|
  (annotation["file"]=="polyA_37_1")|
  (annotation["file"]=="polyA_37_2")
annotation$batch<-as.character(annotation$batch)
annotation[batch_rows,"batch"] = "B"
annotation$batch<-as.factor(annotation$batch)

# finish making df
row.names(annotation$file) 
anno <- annotation[,-1]
rownames(anno) <- annotation[,1]

dds.A <- DESeqDataSetFromMatrix(countData  = round(mat.A), 
                                  colData = annotation, 
                                  design = ~batch + condition)


dds <- DESeq(dds.A)

vst <- varianceStabilizingTransformation(dds, blind=TRUE)
vst_mat <- assay(vst)


## polyA-depletion correlation ----

d <- cor((vst_mat), method="spearman")
hc <- hclust(dist(1-d))
ggdendrogram(hc, rotate = TRUE, theme_dendro = FALSE)

#PCA
# pre-batch removal, color batch
pcaData <- plotPCA(vst, intgroup=c("condition", "batch"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition, shape=batch)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

# use limma to remove batch effect, then redo PCA (DOES THIS look right?)
assay(vst) <- limma::removeBatchEffect(assay(vst), vst$batch)

# post-batch removal, color condition, shape replicate
pcaData <- plotPCA(vst, intgroup=c("condition", "batch"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition, shape=batch)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()


## Use DEG_library_analysis.R to continue --------------------------------------
#polyA DEGS
## 25 vs 34
resA <- results(dds, contrast=c('condition','25','34'))
summary_deseq(resA)
resA <- resA[order(resA$padj), ]
resA_data <- merge(as.data.frame(resA), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resA_data)[1] <- "Gene"

DEG_25v34 <- resA_data %>% dplyr::filter(abs(log2FoldChange) > 1 & padj < 0.05)
tableDEGs <- DEG_25v34[order(DEG_25v34$log2FoldChange), ] %>% dplyr::select("Gene", "log2FoldChange", "padj")
rownames(tableDEGs) = NULL

knitr::kable(tableDEGs, caption = "25 vs 34 polyA DEGs", digits = 3) %>%
  kable_classic(full_width = F, html_font = "Cambria") %>%
  scroll_box(width = "500px", height = "200px")

vsVolcano(
    x = '25', y = '34', 
    data = dds, d.factor = 'condition', type = 'deseq', 
    padj = 0.05, x.lim = NULL, lfc = NULL, title = TRUE, 
    legend = TRUE, grid = TRUE, data.return = FALSE
)


## polyA DEGs accounting for batch - using LRT
# across all conditions not just two?
batch_dds<- DESeq(dds, test = "LRT",reduced = ~batch)
batch_res<-results(batch_dds)
summary_deseq(batch_res)
batch_res_data <- merge(as.data.frame(batch_res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(batch_res_data)[1] <- "Gene"

DEGs <- resA_data %>% dplyr::filter(abs(log2FoldChange) > 1 & padj < 0.05)
lrt_tableDEGs <- DEGs[order(DEGs$log2FoldChange), ] %>% dplyr::select("Gene", "log2FoldChange", "padj")
rownames(tableDEGs) = NULL
#total = 17807 variables

# DEGS for ONLY one comparison at a time- DOES THIS ACCOUNT FOR BATCH? Since we changed the original dds design for polyA

# DEG_25v4 = 36
# DEG_25v10 = 87
# DEG_25v34 = 83
# DEG_25v37 = 402
## total = 608

# what do these commands tell us?
resultsNames(dds)

mm <- model.matrix(~ batch + condition + batch:condition,
                   as.data.frame(colData(dds)))




## enhanced volcano? this works
```{r Enhanced Volcano plots, warning=FALSE, echo=FALSE, cache= T} 
# Repeat for each condition (must highlight and run each separately)
# Hot 37º vs Control 25º
res_37v25 <- lfcShrink(dds, contrast=c('condition','37','25'), type = "normal") 

data<- res_37v25
title<- "Hot 37º vs Control 25º"

# Hot 34º vs Control 25º
res_34v25 <- lfcShrink(dds, contrast=c('condition','34','25'), type = "normal") 

data<- res_34v25
title<- "Hot 34º vs Control 25º"

# Cold 10º vs Control 25º
res_10v25 <- lfcShrink(dds, contrast=c('condition','10','25'), type = "normal") 

data<- res_10v25
title<- "Cold 10º vs Control 25º"

# Cold 4º vs Control 25º
res_4v25 <- lfcShrink(dds, contrast=c('condition','4','25'), type = "normal") 

data<- res_4v25
title<- "Cold 4º vs Control 25º"

## Labelling the plot (legend)

FC <- 1
p <- 5e-2

keyvals <- rep('grey50', nrow(data))
names(keyvals) <- rep('NS', nrow(data))

# keyvals[which(abs(data$log2FoldChange) > FC & data$padj > p)] <- 'grey50'
# names(keyvals)[which(abs(data$log2FoldChange) > FC & data$padj > p)] <- 'log2FoldChange'
# 
# keyvals[which(abs(data$log2FoldChange) < FC & data$padj < p)] <- 'grey50'
# names(keyvals)[which(abs(data$log2FoldChange)  < FC & data$padj < p)] <- '-Log10Q'

# mark the lncRNAs from flybase
lncs<-read.table("/slipstream/home/sethandy/thermofly/lncRNAs_flybase.csv")
keyvals[which(lncs$V1 %in% rownames(data))] <- 'orange'
names(keyvals)[which(lncs$V1 %in% rownames(data))] <- 'lncRNA'

keyvals[which(data$log2FoldChange < -FC & data$padj < p)] <- 'darkcyan'
names(keyvals)[which(data$log2FoldChange  < -FC & data$padj < p)] <- 'Signif. down-regulated'

keyvals[which(data$log2FoldChange > FC & data$padj < p)] <- 'red2'
names(keyvals)[which(data$log2FoldChange > FC & data$padj < p)] <- 'Signif. up-regulated'


unique(keyvals)
unique(names(keyvals))


#plot
volcano <- EnhancedVolcano(data,
                           lab = "",
                          # lab = rownames(data),
                           x = 'log2FoldChange',
                           y = 'pvalue',
                           selectLab = NULL,
                           ylim = c(0,20),
                           xlim = c(-3,3),
                           xlab = bquote(~Log[2]~ 'fold change'),
                           ylab = bquote(~-Log[10] ~ italic(P)),
                           title = paste0(title),
                           pCutoff = 5e-2,
                           FCcutoff = 1,
                           pointSize = 1,
                           subtitle=NULL,
                           labSize = 3.0,
                           labCol = 'black',
                           labFace = 'bold',
                           boxedLabels = F,
                           #shape = c(6, 4, 2, 11, 15),
                           colCustom = keyvals,
                           colAlpha = 1,
                           legendPosition = 'right',
                           legendLabSize = 5,
                           legendIconSize = 5.0,
                           drawConnectors = FALSE,
                           widthConnectors = 0.5,
                           colConnectors = 'grey50',
                           gridlines.major = FALSE,
                           gridlines.minor = FALSE,
                           border = 'partial',
                           borderWidth = 0.5,
                           borderColour = 'black')

volcano

