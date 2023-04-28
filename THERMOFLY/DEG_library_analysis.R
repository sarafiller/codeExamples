#-------------------------------------------------------------------------------
# Thermofly: RNA-seq library DEG analysis (of D. mel. acute thermal shock)
# Comparing polyA enhancement to rRNA depletion for SALMON transcript data
# July 18, 2022
# Sara Filler
#-------------------------------------------------------------------------------

# Load packages ----
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

require(ggpattern)
require(ggrepel)

# Load files ---- 

# for STAR data
#setwd("/Users/sarafiller/Desktop/Thermofly/reu_workshop/saraNivea")
# # polyA
# polyA<-readRDS("polyA_dds.rds")
# # rRNA
# rRNA<-readRDS("reu_workshop_nivea_dds_rRNA.rds")

# for SALMON data
txdb <- makeTxDbFromGFF("~/../joeboyd/indexes/DM6/GTF/dm6.ensGene.gtf")
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")
head(tx2gene)

dir = "/slipstream_old/home/sethandy/sethandy/thermofly/salmon"

sample <- list.files(dir) %>% str_replace_all( ".salmon_quant", "")
library = sapply(strsplit(sample, '_'), `[`, 1)
condition = sapply(strsplit(sample, '_'), `[`, 2)
condition<- factor(condition,levels=c("25","37","34","10","4"))
replicate = sapply(strsplit(sample, '_'), `[`, 3)
directory = list.files(dir)
# Make df ----
annotation <- data.frame(sample = sample, directory= directory, library = library, condition = condition, replicate = replicate, batch=c("A"), stringsAsFactors = TRUE)

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
# getting rid of outlier extra sample
salmon.df = as.data.frame(txi.salmon$counts) %>% dplyr::select(-rRNA_10_4)
salmon.mat = as.matrix(salmon.df)
anno = annotation[annotation$sample != "rRNA_10_4",]
rownames(anno) = colnames(salmon.mat)

#setup DESeq2 design and run ----
dds_all <- DESeqDataSetFromMatrix(round(salmon.mat), 
                                  colData = anno, 
                                  design = ~batch + condition)

# keep only rows where row sums are >= 10 
keep <- rowSums(counts(dds_all)) >= 10
dds_all<- dds_all[keep,]

dds_all <- DESeq(dds_all)

# check that this is true
all(colnames(dds_all) == colnames(salmon.mat))

# only need for qc
# vst_all <- varianceStabilizingTransformation(dds_all, blind=TRUE)
# vst_mat_all <- assay(vst_all)

# make polyA dds object
# has a batch
polyA_salmon.df<-salmon.df %>% dplyr:: select(starts_with("polyA"))
polyA_salmon.mat<-as.matrix(polyA_salmon.df)
polyA_anno<-anno[anno$library=="polyA",]
rownames(polyA_anno) = colnames(polyA_salmon.mat)
dds_polyA<- DESeqDataSetFromMatrix(round(polyA_salmon.mat), 
                                  colData = polyA_anno, 
                                  design = ~batch + condition)

# keep only rows where row sums are >= 10 
keep <- rowSums(counts(dds_polyA)) >= 10
dds_polyA<- dds_polyA[keep,]

dds_polyA <- DESeq(dds_polyA)

# save as RDS object so dds object is easily loaded later
# wd = "/slipstream_old/home/sethandy/sethandy/thermofly"
# setwd(wd)
# saveRDS(dds_polyA,"polyA_salmon_dds.RDS")
# # readRDS("/slipstream_old/home/sethandy/sethandy/thermofly/polyA_salmon_dds.rds")


# make rRNA dds object
# no batch
rRNA_salmon.df<-salmon.df %>% dplyr:: select(starts_with("rRNA"))
rRNA_salmon.mat<-as.matrix(rRNA_salmon.df)
rRNA_anno<-anno[anno$library=="rRNA",]
rownames(rRNA_anno) = colnames(rRNA_salmon.mat)
dds_rRNA<- DESeqDataSetFromMatrix(round(rRNA_salmon.mat), 
                                   colData = rRNA_anno, 
                                   design = ~condition)

# keep only rows where row sums are >= 10 
keep <- rowSums(counts(dds_rRNA)) >= 10
dds_rRNA<- dds_rRNA[keep,]

dds_rRNA <- DESeq(dds_rRNA)

# Save as RDS object
# wd = "/slipstream_old/home/sethandy/sethandy/thermofly"
# setwd(wd)
# saveRDS(dds_rRNA,"rRNA_salmon_dds.RDS")
# # readRDS("/slipstream_old/home/sethandy/sethandy/thermofly/rRNA_salmon_dds.rds")


# Load functions ----
res_function<-function(dds){
  # Analyze data
  res_37<-results(dds,
                  contrast=c('condition','37','25'),
                  #name = "condition_Hot37_vs_Control",
                  alpha=0.05) # can do here or can set padj value later
  
  res_34<-results(dds,
                  contrast=c('condition','34','25'),
                  #name = "condition_Hot34_vs_Control",
                  alpha=0.05)
  res_10<-results(dds,
                  contrast=c('condition','10','25'),
                  #name = "condition_Cold10_vs_Control",
                  alpha=0.05)
  res_4<-results(dds,
                 contrast=c('condition','4','25'),
                 #name = "condition_Cold4_vs_Control",
                 alpha=0.05)
  
  # p adjusted value - after making a comparison for every gene, chances of genes 
  # being different is much more likely to be due to chance
  # so p adj makes it more stringent, so what you find is actually different
  
  # Take a peak at summary values
  # summary(res_37)
  # summary(res_34)
  # summary(res_10)
  # summary(res_4)
  
  # convert to Tibbles (like dfs but better)
  res_37 <- res_37 %>%
    as_tibble(rownames="FBgn") #genes row names are lost, unless rownames are specified
  res_34 <- res_34 %>%
    as_tibble(rownames="FBgn")
  res_10 <- res_10 %>%
    as_tibble(rownames="FBgn")
  res_4 <- res_4 %>%
    as_tibble(rownames="FBgn")
  
  # Join all these tibbles together
  res <- bind_rows(list("res_37"=res_37,
                        "res_34"=res_34,
                        "res_10"=res_10,
                        "res_4"=res_4),
                   .id="temp") # a column is created to link each row to its original df

  return(res)
}

# Make DF ----

# DEG_25v34 <- resA_data %>% dplyr::filter(abs(log2FoldChange) > 1 & padj < 0.05)
# tableDEGs <- DEG_25v34[order(DEG_25v34$log2FoldChange), ] %>% dplyr::select("Gene", "log2FoldChange", "padj")
# rownames(tableDEGs) = NULL

# Making a results df for rRNA and polyA, then join
# splitting them up is necessary because we are doing intra-library comparisons
# polyA has a batch effect
# rRNA does not
polyA_res<-res_function(dds_polyA)
rRNA_res<-res_function(dds_rRNA)

# bind and get rid of rows with missing values
total_res <- bind_rows(list("polyA"=polyA_res,
                            "rRNA"=rRNA_res),
                       .id="library")

# Converting to gene symbol from FBgn
#install.packages("BiocManager")
# BiocManager::install("AnnotationDbi")
# BiocManager::install("org.Dm.eg.db")

library(BiocManager)
library(AnnotationDbi) # annotation DB
library(org.Dm.eg.db) # specifically D. melanogaster DB
total_res$symbol <-mapIds(org.Dm.eg.db,
                          keys = total_res$FBgn,
                          column="SYMBOL",
                          keytype="FLYBASE",
                          multiVals="first")


# get rid of na's later (during plotting) via filter(padj<0.05)
# total_res_com<- total_res[complete.cases(total_res),] # was in Sophie pipeline
# dim(total_res_com)
# total_res_noNa<-total_res%>%filter(!is.na(padj))
# dim(total_res_noNa)

# DEG counts ----
# find all p values less than .05
FC <- 1

filtered_deg_counts<-total_res %>%
              filter(padj<0.05) %>%
              filter(abs(log2FoldChange)>FC) %>% # 1 - more stringent vs 0
              group_by(library, temp, log2FoldChange > 1) %>% # false=downreg/true=upreg
              #group_by(library)%>%
              tally() %>%
              dplyr::rename(direction = `log2FoldChange > 1`)%>%
              mutate(direction = ifelse(direction, "up", "down"))

# UNfiltered degs
# polyA 49160
# rRNA  46152

# filtered degs, not unique
# polyA 562
# rRNA  855

# filtered degs, UNIQUE
# polyA 448
# rRNA  590

# lncRNA counts ----

## Use lncRNAs_flybase.csv to go through
lncs <-read.table("/slipstream_old/home/sethandy/sethandy/thermofly/lncRNAs_flybase.csv")

lncs_total_res<-total_res[which(total_res$FBgn %in% lncs$V1),]
# do instead of filtering total_res with
 #filter(str_detect(symbol,"lnc")) %>%

# filter more to plot
lncs_plot<-lncs_total_res%>%
  filter(padj<0.05) %>%
  filter(abs(log2FoldChange)>FC) %>% # 1 - more stringent vs 0
  group_by(temp, library, log2FoldChange > 1) %>% # false=downreg/true=upreg
  #group_by(library)%>% # to get numbers per lib
  tally() %>%
  dplyr::rename(direction = `log2FoldChange > 1`)%>%
  mutate(direction = ifelse(direction, "up", "down"))

# add these rows to make plot nicer
new_row1<-list("res_4","rRNA","down",as.integer(0))
new_row2<-list("res_4","rRNA","up",as.integer(0))
lncs_plot[nrow(lncs_plot) + 1,] <- new_row1
lncs_plot[nrow(lncs_plot) + 1,] <- new_row2

# filtered, NOT unique lnc DEGs
# polyA 27
# rRNA  23

# UNfiltered lncRNAs
# polyA 2908
# rRNA  1988

# filtered, unique lncRNAs
# polyA 20
# rRNA  19

# IN COMMON ----
# common lncsDEGs ----

# not filtered  (if not filtering total DEGs)
#lncs<- lncs_total_res # unfiltered by anything

# filtered
lncs<-  lncs_total_res %>%
        filter(padj<0.05) %>%
        filter(abs(log2FoldChange)>FC)
polyA_lncs<-lncs%>%
  filter(library=="polyA")

rRNA_lncs<-lncs%>%
  filter(library=="rRNA")

#Can do this to get unique numbers but ggVennDiagram does automatically
# unique_lncs_polyA<-unique(polyA_lncs$FBgn)
# unique_lncs_rRNA<-unique(rRNA_lncs$FBgn)

x<-list(`polyA enhancement`=polyA_lncs$FBgn,`rRNA depletion`=rRNA_lncs$FBgn)
ggVennDiagram(x)

# common DEGs ----
degs<-total_res%>%
  filter(padj<0.05) %>%
  filter(abs(log2FoldChange)>FC) 

polyA_degs<-  degs%>%
              filter(library=="polyA")

rRNA_degs<- degs%>%
            filter(library=="rRNA") 

unique_degs_polyA<-unique(polyA_degs$FBgn)
unique_degs_rRNA<-unique(rRNA_degs$FBgn)



# venn diagram
x<-list(`polyA enhancement`= polyA_degs$FBgn,`rRNA depletion`=rRNA_degs$FBgn)
ggVennDiagram(x)

same_degs<-polyA_degs[which(polyA_degs$FBgn %in% rRNA_degs$FBgn),]
#unique(same_degs$FBgn)

# find and plot where in common degs are found
category_same_degs<-same_degs %>%
  group_by(library,temp, log2FoldChange > 1) %>% # false=downreg/true=upreg
  tally() %>%
  dplyr::rename(direction = `log2FoldChange > 1`)%>%
  mutate(direction = ifelse(direction, "up", "down"))

# plot
category_same_degs%>%
  mutate(library = ifelse(library == "polyA", 
                          "polyA enhancement", 
                          "rRNA depletion")) %>%
  mutate(temp = case_when(temp=="res_4"~ "4°C",
                          temp=="res_10"~ "10°C",
                          temp=="res_34"~ "34°C",
                          temp=="res_37"~ "37°C"))%>%
  mutate(temp=factor(temp,levels=c("4°C",
                                   "10°C",
                                   "34°C",
                                   "37°C")))%>%
  dplyr::rename(`Number of DEGs found in common by both libraries`=n)%>%
  ggplot(aes(y = temp,
             x =`Number of DEGs found in common by both libraries`,
             fill = direction)) +
  geom_bar(stat = "identity",
           position = "stack") +
  geom_vline(xintercept = 0, color = "grey20") +
  scale_fill_manual(values = c("#457b9d", "#e63946")) +
  theme_classic(base_size = 12) +
  scale_x_continuous(breaks = seq(from = 0, to = 200, by = 20),
                     label = abs(seq(from = 0, to = 200, by = 20)),
                     name = "Number of differentially expressed genes",
                     position = "top",
                     limits=c(0,200)) +
  ylab("Acute shock temperature relative to control 25°C") +
  scale_y_discrete(drop=FALSE) +
  theme(legend.position = "right",
        axis.line.y = element_blank(),
        axis.line.x = element_line(color = "grey50"),
        axis.ticks = element_line(color = "grey50"),
        panel.grid.major.x = element_line(color = "grey70", 
                                          size = 0.5),
        panel.grid.minor.x = element_line(color = "grey90", 
                                          size = 0.5))

# Stacked bar plots -----

# for all DEGs

filtered_deg_counts %>%
  mutate(library = ifelse(library == "polyA", 
                        "polyA enhancement", 
                        "rRNA depletion")) %>%
  mutate(temp = case_when(temp=="res_4"~ "4°C",
                          temp=="res_10"~ "10°C",
                          temp=="res_34"~ "34°C",
                          temp=="res_37"~ "37°C"))%>%
  mutate(temp=factor(temp,levels=c("4°C",
                                   "10°C",
                                   "34°C",
                                   "37°C")))%>%
  mutate(DEGs = ifelse(direction == "up", n, -n)) %>%
  ggplot(aes(y = temp,
             x = DEGs,
             fill = direction,
             pattern = library)) +
  annotate("rect",
           fill = "grey95",
           xmin = -Inf, xmax = Inf,
           ymin = 0.6, 
           ymax = 1.4,
           alpha = 0.5) +
  annotate("rect",
           fill = "grey95",
           xmin = -Inf, xmax = Inf,
           ymin = 1.6, 
           ymax = 2.4,
           alpha = 0.5) +
  annotate("rect",
           fill = "grey95",
           xmin = -Inf, xmax = Inf,
           ymin = 2.6, 
           ymax = 3.4,
           alpha = 0.5) +
  annotate("rect",
           fill = "grey95",
           xmin = -Inf, xmax = Inf,
           ymin = 3.6, 
           ymax = 4.4,
           alpha = 0.5) +
  ggpattern::geom_bar_pattern(aes(group = library), pattern_color = NA,
                              color = "grey50",
                              width = 0.65,
                              stat = "identity",
                              position = position_dodge(width = 0.7),
                              pattern_fill = "grey80",
                              pattern_angle = 45,
                              pattern_density = 0.3,
                              pattern_spacing = 0.0125,
                              pattern_key_scale_factor = .5) +
  geom_vline(xintercept = 0, color = "grey20") +
  scale_fill_manual(values = c("#457b9d", "#e63946")) +
  theme_classic(base_size = 12) +
  scale_x_continuous(breaks = seq(from = -400, to = 400, by = 100),
                     label = abs(seq(from = -400, to = 400, by = 100)),
                     name = "Number of differentially expressed genes",
                     position = "top",
                     limits=c(-400,400)) +
  ggpattern::scale_pattern_manual(values = c(`polyA enhancement` = "none", 
                                             `rRNA depletion` = "stripe")) +
  ylab("Acute shock temperature relative to control 25°C") +
  scale_y_discrete(drop = FALSE) +
  theme(legend.position = "right",
        axis.line.y = element_blank(),
        axis.line.x = element_line(color = "grey50"),
        axis.ticks = element_line(color = "grey50"),
        panel.grid.major.x = element_line(color = "grey70", 
                                          size = 0.5),
        panel.grid.minor.x = element_line(color = "grey90", 
                                          size = 0.5)) 

# for only lncDEGs

lncs_plot %>%
  mutate(library = ifelse(library == "polyA", 
                          "polyA enhancement", 
                          "rRNA depletion")) %>%
  mutate(temp = case_when(temp=="res_4"~ "4°C",
                          temp=="res_10"~ "10°C",
                          temp=="res_34"~ "34°C",
                          temp=="res_37"~ "37°C"))%>%
  mutate(temp=factor(temp,levels=c("4°C",
                                   "10°C",
                                   "34°C",
                                   "37°C")))%>%
  mutate(DEGs = ifelse(direction == "up", n, -n)) %>%
  ggplot(aes(y = temp,
             x = DEGs,
             fill = direction,
             pattern = library)) +
  annotate("rect",
           fill = "grey95",
           xmin = -Inf, xmax = Inf,
           ymin = 0.6, 
           ymax = 1.4,
           alpha = 0.5) +
  annotate("rect",
           fill = "grey95",
           xmin = -Inf, xmax = Inf,
           ymin = 1.6, 
           ymax = 2.4,
           alpha = 0.5) +
  annotate("rect",
           fill = "grey95",
           xmin = -Inf, xmax = Inf,
           ymin = 2.6, 
           ymax = 3.4,
           alpha = 0.5) +
  annotate("rect",
           fill = "grey95",
           xmin = -Inf, xmax = Inf,
           ymin = 3.6, 
           ymax = 4.4,
           alpha = 0.5) +
  ggpattern::geom_bar_pattern(aes(group = library), pattern_color = NA,
                              color = "grey50",
                              width = 0.65,
                              stat = "identity",
                              position = position_dodge(width = 0.7),
                              pattern_fill = "grey80",
                              pattern_angle = 45,
                              pattern_density = 0.3,
                              pattern_spacing = 0.0125,
                              pattern_key_scale_factor = .5) +
  geom_vline(xintercept = 0, color = "grey20") +
  scale_fill_manual(values = c("#457b9d", "#e63946")) +
  theme_classic(base_size = 12) +
  scale_x_continuous(breaks = seq(from = -400, to = 400, by = 2),
                     label = abs(seq(from = -400, to = 400, by = 2)),
                     name = "Number of differentially expressed genes encoding for lncRNA",
                     position = "top", limits=c(-10,10)
                     ) +
  ggpattern::scale_pattern_manual(values = c(`polyA enhancement` = "none", 
                                             `rRNA depletion` = "stripe")) +
  ylab("Acute shock temperature relative to control 25°C") +
  scale_y_discrete(drop = FALSE) +
  theme(legend.position = "right",
        axis.line.y = element_blank(),
        axis.line.x = element_line(color = "grey50"),
        axis.ticks = element_line(color = "grey50"),
        panel.grid.major.x = element_line(color = "grey70", 
                                          size = 0.5),
        panel.grid.minor.x = element_line(color = "grey90", 
                                          size = 0.5)) 

# Volcano Plots ----
total_res%>%
  #filter(temp=="res_37") %>%
  #filter(temp=="res_10" | temp=="res_4") %>%
  #filter(library=="polyA")%>%
  #filter(str_detect(symbol,"lnc"))%>%   # filters to just lncs
  mutate(temp = case_when(temp=="res_4"~ "4°C vs 25°C",
                          temp=="res_10"~ "10°C vs 25°C",
                          temp=="res_34"~ "34°C vs 25°C",
                          temp=="res_37"~ "37°C vs 25°C"))%>%
  mutate(temp=factor(temp,levels=c("4°C vs 25°C",
                                   "10°C vs 25°C",
                                   "34°C vs 25°C",
                                   "37°C vs 25°C")))%>%
  filter(!is.na(padj)) %>% # na is genes that aren't expressed at high enough level
  ggplot()+
  geom_point(aes(x=log2FoldChange,y=-log10(padj),
                 color= padj <0.05),
             size= 0.5) + 
  scale_color_manual(values = c("grey70",
                                "firebrick"),
                     name = "Differentially expressed genes",
                     label = c("Not significant", "Significant"))+
  geom_vline(xintercept=0,
             color="grey70",
             linetype= 2)+
  geom_hline(yintercept=-log10(0.05),
             color= "grey 70",
             linetype=2)+
  #xlim(c(-14,14))+
  #ylim(c(0,150))+
  theme_classic()+
  theme(legend.position="top")+
  facet_wrap(~temp~library, nrow=2, ncol=4)

# some genes have p values super close to 0, so they are reaching infinity on the plot
total_res$padj[order(total_res$padj)]

# Fold changes 
# log2FoldChange > 1 - things that are increasing by double or more
# log2FoldChange < -1 - things that are decreasing by half or more
# enhanced volcano plots these lines automatically

# use ggrepel to make nice labels that don't overlap
#focus on 25v37 volcano and facet by library
total_res%>%
  filter(temp=="res_37") %>%
  filter(!is.na(padj)) %>% # na is genes that aren't expressed at high enough level
  ggplot(aes(x=log2FoldChange,y=-log10(padj)))+
  geom_point(aes(color= padj <0.05)) + 
  geom_vline(xintercept=0,
             color="grey70",
             linetype= 2)+
  geom_vline(xintercept=-1,
             color="grey70",
             linetype= 2)+
  geom_vline(xintercept=1,
             color="grey70",
             linetype= 2)+
  geom_hline(yintercept=-log10(0.05),
             color= "grey 70",
             linetype=2)+
  geom_label_repel(data=total_res %>%
                     filter(temp=="res_37"& library=="polyA") %>%
                     arrange(padj)%>%
                     head(10), # get top 10 lowest p values
                   aes(label=symbol),
                   min.segment.length = unit(0, 'lines'),
                   max.iter = 30000,
                   force=10)+ # add lines
  geom_label_repel(data=total_res %>%
                     filter(temp=="res_37"& library=="rRNA") %>%
                     arrange(padj)%>%
                     head(10), # get top 10 lowest p values
                   aes(label=symbol),
                   min.segment.length = unit(0, 'lines'),
                   max.iter=30000,
                   force = 10)+ # add lines
  scale_color_manual(values = c("grey70",
                                "firebrick"),
                     name = "Differentially expressed genes",
                     label = c("Not significant", "Significant"))+
  xlim(c(-12,12))+
  theme_classic()+
  theme(legend.position="top")+
  facet_wrap(~library)

# ind DEG analysis ----
# Example of individual comparison (Seth/Sophie way) from DESeq2_salmon.Rmd script
# 25 vs 34
resA <- results(dds, contrast=c('condition','25','34'))
summary_deseq(resA)
resA <- resA[order(resA$padj), ]
resA_data <- merge(as.data.frame(resA), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resA_data)[1] <- "Gene"

DEG_25v34 <- resA_data %>% 
             dplyr::filter(abs(log2FoldChange) > 1 & padj < 0.05)
tableDEGs <- DEG_25v34[order(DEG_25v34$log2FoldChange), ] %>% 
             dplyr::select("Gene", "log2FoldChange", "padj")
rownames(tableDEGs) = NULL

# Volcano for polyA Hot34vControl25
vsVolcano(
  x = '25', y = '34', 
  data = dds, d.factor = 'condition', type = 'deseq', 
  padj = 0.05, x.lim = NULL, lfc = NULL, title = TRUE, 
  legend = TRUE, grid = TRUE, data.return = FALSE
)
# ----
