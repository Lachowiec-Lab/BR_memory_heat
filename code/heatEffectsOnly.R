
### Libraries ####

library(dplyr)
library(DESeq2)
library(RColorBrewer)
library(viridis)
library(pheatmap)
library(goseq)
library(stringr)
library(VennDiagram)
library(simplifyEnrichment)
library(biomaRt)

### Functions ####
GOenrichTable <- function(gene_df) {
  # Create df of gene ID, DE and length - keep unique
  de_len <- gene_df %>% 
    dplyr::select(gene_id, is_de, gen_len) %>%
    unique()
  # Gene vector
  de_vector <- de_len$is_de
  names(de_vector) <- de_len$gene_id
  #GO analysis
  pwf <- nullp(de_vector, bias.data = de_len[,c("gen_len")])
  go_tab <- goseq(pwf=pwf, gene2cat = gene_df[,c("gene_id", "go_id")]) %>%
    dplyr::mutate(padj_over = p.adjust(over_represented_pvalue, method="BH")) %>%
    dplyr::filter(over_represented_pvalue < 0.05)
  return(go_tab)
}

# function for the color scale on heatmaps
quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

# Calculating means for each treatment
makeMeans <- function(countsAll) {
  C0 <- c(grep("CN0", colnames(countsAll), value = TRUE), grep("CY0", colnames(countsAll), value = TRUE))
  mock_preHeat_time0 <- rowMeans(countsAll[,C0])
  S0 <- c(grep("SN0", colnames(countsAll), value = TRUE), grep("SY0", colnames(countsAll), value = TRUE))
  seed_preHeat_time0 <- rowMeans(countsAll[,S0])
  
  CN8 <- c(grep("CN8", colnames(countsAll), value = TRUE))
  mock_noHeat_time8 <- rowMeans(countsAll[,CN8])
  CY8 <- c(grep("CY8", colnames(countsAll), value = TRUE))
  mock_heat_time8 <- rowMeans(countsAll[,CY8])
  SN8 <- c(grep("SN8", colnames(countsAll), value = TRUE))
  seed_noHeat_time8 <- rowMeans(countsAll[,SN8])
  SY8 <- c(grep("SY8", colnames(countsAll), value = TRUE))
  seed_heat_time8 <- rowMeans(countsAll[,SY8])
  
  CN16 <- c(grep("CN16", colnames(countsAll), value = TRUE))
  mock_noHeat_time16 <- rowMeans(countsAll[,CN16])
  CY16 <- c(grep("CY16", colnames(countsAll), value = TRUE))
  mock_heat_time16 <- rowMeans(countsAll[,CY16])
  SN16 <- c(grep("SN16", colnames(countsAll), value = TRUE))
  seed_noHeat_time16 <- rowMeans(countsAll[,SN16])
  SY16 <- c(grep("SY16", colnames(countsAll), value = TRUE))
  seed_heat_time16 <- rowMeans(countsAll[,SY16])
  
  CN24 <- c(grep("CN24", colnames(countsAll), value = TRUE))
  mock_noHeat_time24 <- rowMeans(countsAll[,CN24])
  CY24 <- c(grep("CY24", colnames(countsAll), value = TRUE))
  mock_heat_time24 <- rowMeans(countsAll[,CY24])
  SN24 <- c(grep("SN24", colnames(countsAll), value = TRUE))
  seed_noHeat_time24 <- rowMeans(countsAll[,SN24])
  SY24 <- c(grep("SY24", colnames(countsAll), value = TRUE))
  seed_heat_time24 <- rowMeans(countsAll[,SY24])
  
  countsMean <- data.frame(mock_preHeat_time0, seed_preHeat_time0, mock_noHeat_time8, seed_noHeat_time8, 
                           mock_noHeat_time16, seed_noHeat_time16, mock_noHeat_time24, seed_noHeat_time24,
                           mock_heat_time8, seed_heat_time8, 
                           mock_heat_time16, seed_heat_time16, 
                           mock_heat_time24, seed_heat_time24)
  return(countsMean)
}


### Read in Data ####
# info for GO term enrichment
features <- read.csv("GOterms.csv")

countsRaw <- read.table("RNA-seq-counts.txt", header=T, row.name=1)

##Only using high confidence genes; remove rows containing LC
selected_rows <- grep("LC", rownames(countsRaw), value = TRUE)
countsRaw <- countsRaw[!rownames(countsRaw) %in% selected_rows,]

expdesignRaw <- read.table("expDesign7.txt", row.names=1, header=T)
expdesignRaw$trt <- factor(expdesignRaw$trt)
expdesignRaw$trt <- relevel(expdesignRaw$trt, "ck")
expdesignRaw$stress <- factor(expdesignRaw$stress, levels = c("no_heat","heat"))
expdesignRaw$stress <- relevel(expdesignRaw$stress, "no_heat")
expdesignRaw$time <- factor(expdesignRaw$time, levels = c("0", "8", "16", "24"))
expdesignRaw$time <- relevel(expdesignRaw$time, "0")

levels(expdesignRaw$trt)
levels(expdesignRaw$stress)
levels(expdesignRaw$time)

### Preparing for plotting genes####
# Normalizing read counts
ddsAll <- DESeqDataSetFromMatrix(countData=countsRaw, colData=expdesignRaw, design= ~ time)
smallestGroupSize <- 3
keep <- rowSums(counts(ddsAll) >= 10) >= smallestGroupSize
ddsAll  <- ddsAll[keep,]
ddsAll <- estimateSizeFactors(ddsAll)
countsAll <- counts(ddsAll, normalized=TRUE)
countsMean <- makeMeans(countsAll)

### Main effect of heat ####
#####All timepoints used

## Data clean up--we need to rename the timepoint zero as all being no heat 
# because there is no heat applied at time 0
expdesign <- expdesignRaw %>%
  as.data.frame() %>%
  mutate(trueStress = ifelse(stress == "heat" & time != 0, "heat", "no_heat"))

expdesign$trueStress <- factor(expdesign$trueStress, levels = c("no_heat","heat"))
expdesign$trueStress <- relevel(expdesign$trueStress, "no_heat")

counts <- countsRaw
ncol(counts) == nrow(expdesign)

##### DESeq2 differential expression ####
dds <- DESeqDataSetFromMatrix(countData=counts, colData=expdesign, design = ~ trueStress + trt + time)

smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds  <- dds[keep,]
dds <- DESeq(dds)
resultsNames(dds)
resHeat <- results(dds, name="trueStress_heat_vs_no_heat", alpha = 0.05)
resHeat
summary(resHeat)

# GO term enrichment on significant
# up-regulated and down-regulated genes

# Identifying relevant DEGs
res_de_tab_upreg <- resHeat%>% 
  as.data.frame() %>%
  filter(log2FoldChange > 0) %>%
  mutate(is_de=ifelse(padj >= 0.05|is.na(padj), 0, 1),
         gene_id=rownames(.)) 

res_de_tab_dnreg <- resHeat %>% 
  as.data.frame() %>%
  filter(log2FoldChange < 0) %>%
  mutate(is_de=ifelse(padj>=0.05|is.na(padj), 0, 1),
         gene_id=rownames(.)) 

# Merging GO terms with DEGs
heat_up_df <- res_de_tab_upreg %>%
  left_join(features, by=c("gene_id"="ensembl_gene_id")) %>%
  mutate(gen_len=end_position-start_position)

heat_down_df <- res_de_tab_dnreg %>%
  left_join(features, by=c("gene_id"="ensembl_gene_id")) %>%
  mutate(gen_len=end_position-start_position)

##### GO enrichment ####
heatDownGOtable <- GOenrichTable(heat_down_df)
heatUpGOtable <- GOenrichTable(heat_up_df)

heatDownGOtable %>%
  filter(padj_over < 0.05 & ontology == "BP")
heatUpGOtable %>%
  filter(padj_over < 0.05 & ontology == "BP")

##### Heat shock gene expression ####
## Filter for HSP annotated genes differentially expressed
hsps <- features %>%
  filter(str_detect(description, "Heat shock protein|Hsp|HSP")|str_detect(name_1006, "Heat shock protein|Hsp"))
resHeat_hsps <- resHeat[rownames(resHeat) %in% hsps$ensembl_gene_id,]

nrow(resHeat_hsps[resHeat_hsps$padj < 0.05,]) # 34 show differential expression
nrow(resHeat_hsps[resHeat_hsps$padj < 0.05 & resHeat_hsps$log2FoldChange > 0,]) # 30 show up-regulation
diffExpHsps <- resHeat_hsps[resHeat_hsps$padj < 0.05,] 
hspsCounts <- countsMean[rownames(countsMean) %in% rownames(diffExpHsps),]

## Preparing for plotting annotations
hspsGroups <- unique(features[features$ensembl_gene_id %in% rownames(hspsCounts), c("ensembl_gene_id","description")])
rownames(hspsGroups) <- hspsGroups$ensembl_gene_id
hspsGroups <- hspsGroups[-1]
hspsGroups$description <- gsub("\\[[^][]*]", "", hspsGroups$description)
hspsGroups$description <- sub("^$", "Hsp90 or Hsp70 binding protein", hspsGroups$description)

ann_colors = list(
  description = c(`Hsp90 or Hsp70 binding protein` = "navy",
                  `Co-chaperone protein p23-1 ` = "blue", `Heat shock protein 101 ` = "darkorange", 
                  `Heat shock protein 90 ` = "darkgreen", `Heat shock protein 90-6, mitochondrial ` = "green",
                  `HSP70 ` = "yellow", `RAR1 ` = "lightblue",
                  `Similar to HSP protein (Fragment) ` = "violet", `Small heat shock protein HSP17.8 ` = "gray"))


mat_breaks <- quantile_breaks(as.matrix(hspsCounts), n = 10)

######Figure 3B #####
pheatmap(hspsCounts, color = brewer.pal(9, "Reds"),
         #breaks = mat_breaks, 
         show_rownames = F,
         annotation_row = hspsGroups, 
         annotation_colors = ann_colors,
         cluster_cols = T,
         border_color = "NA",
         cutree_rows = 2, 
         cutree_cols = 2,
         scale = "row"
)

#### Histone gene expression ####
###Histones --this may not be a part of this section...
histones <- features %>%
  filter(go_id == "GO:0000786") %>%
  filter(str_detect(description, "Histone|histone|H2A.5|H2A.7")) # should be thorough

length(unique(histones$ensembl_gene_id)) # no gene id's are repeated
resHeat_histones <- resHeat[rownames(resHeat) %in% histones$ensembl_gene_id,]
nrow(resHeat_histones) #433 expressed of the 456 histones found

nrow(resHeat_histones[resHeat_histones$padj < 0.05,]) # 185 show differential expression
nrow(resHeat_histones[resHeat_histones$padj < 0.05 & resHeat_histones$log2FoldChange < 0,]) # 183 show down-regulation

### Heat effect over time ####
#### Timepoint 0 v 8 #####
# subset data
###### 0 v 8 no heat ####
countsTemp <- countsRaw[, which(names(countsRaw) %in% c("CN01", "CN02", "CN03", "CN81", "CN82", "CN83"))]
expdesignTemp <- expdesignRaw[which(names(countsRaw) %in% c("CN01", "CN02", "CN03", "CN81", "CN82", "CN83")),]
expdesignTemp

dds <- DESeqDataSetFromMatrix(countData = countsTemp, colData = expdesignTemp, design = ~ time)
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds  <- dds[keep,]
dds <- DESeq(dds)
resultsNames(dds)

res0_8nHeat <- results(dds, alpha = 0.05)

summary(res0_8nHeat)
sum(res0_8nHeat$padj < 0.05, na.rm=TRUE) #4028 DEGs due to time only

###### 0 v 8 heat ####
countsTemp <- countsRaw[, which(names(countsRaw) %in% c("CN01", "CN02", "CN03", "CY81", "CY82", "CY83"))]
expdesignTemp <- expdesignRaw[which(names(countsRaw) %in% c("CN01", "CN02", "CN03", "CY81", "CY82", "CY83")),]

dds <- DESeqDataSetFromMatrix(countData=countsTemp, colData=expdesignTemp, design= ~ time)
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds  <- dds[keep,]
dds <- DESeq(dds)
resultsNames(dds)

res0_8Heat <- results(dds, alpha = 0.05)
summary(res0_8Heat)
sum(res0_8Heat$padj < 0.05, na.rm=TRUE) #11736 DEGs due to time and heat

###### Intersection no heat versus heat ####
DEGUpNoHeat <- as.data.frame(res0_8nHeat) %>%
  filter(padj < 0.05, log2FoldChange > 0) 
nrow(DEGUpNoHeat) #1235 DEGs

DEGUpHeat <- as.data.frame(res0_8Heat) %>%
  filter(padj < 0.05, log2FoldChange > 0) 
nrow(DEGUpHeat) #4882 DEGs

overlap1 <- calculate.overlap(x = list(rownames(DEGUpNoHeat), rownames(DEGUpHeat)))
str(overlap1)

##The genes that are only upregulated due to heat (not time) are of interest
# This is captured by not including the intersection genes
overlap1[3] #list of intersecting genes
uniqueHeat0_8upDEGs <- unlist(overlap1[2])[!unlist(overlap1[2]) %in% unlist(overlap1[3])]

DEGDownNoHeat <- as.data.frame(res0_8nHeat) %>%
  filter(padj < 0.05, log2FoldChange < 0) #2793 DEGs
nrow(DEGDownNoHeat)

DEGDownHeat <- as.data.frame(res0_8Heat) %>%
  filter(padj < 0.05, log2FoldChange < 0) #6854 DEGs
nrow(DEGDownHeat)

overlap2 <- calculate.overlap(x = list(rownames(DEGDownNoHeat), rownames(DEGDownHeat)))
str(overlap2)
##The genes that are only downregulated due to heat (not time) are of interest
# This is captured by not including the intersection genes
overlap2[3] #list of intersecting genes
uniqueHeat0_8downDEGs <- unlist(overlap2[2])[!unlist(overlap2[2]) %in% unlist(overlap2[3])]

#####GO terms 0 v 8 unique to heat ####

######Down-regulated GO terms ####
## These are taken from the heat DEGs but only those that ARE NOT also DEGs without heat

# First we need to make a table that indicates which genes have differential expression of interest
# In this case, of interest are genes that are uniquely differentially expressed due to heat not time
# The "base" set are those that are expressed,

# This creates adds a column to the table indicating the genes of interest with a 1 
# and those not with a 0
uniq_0_8Heat_dnreg <- res0_8Heat %>% 
  as.data.frame() %>%
  mutate(is_de=ifelse(rownames(res0_8Heat) %in% uniqueHeat0_8downDEGs, 1, 0),
         gene_id=rownames(.)) 

# Now we merge this dataframe with the features dataframe that includes gene lengths and GO terms
# Genes can have multiple GO terms associated with them
heat_down_df <- uniq_0_8Heat_dnreg %>%
  left_join(features, by=c("gene_id"="ensembl_gene_id")) %>%
  mutate(gen_len=end_position-start_position)

# Complete GO enrichment and output a table
# Uses goseq approach
heatDownGOtable8 <- GOenrichTable(heat_down_df)

dim(heatDownGOtable8 %>% 
      mutate(term=forcats::fct_reorder(.f = term,.x = numInCat, .fun = sort)) %>% 
      filter(ontology=="BP", padj_over < 0.05))
# 101 terms total enriched

######Up-regulated GO terms ####
## These are taken from the heat DEGs but only those that ARE NOT also DEGs without heat
uniq_0_8Heat_upreg <- res0_8Heat %>% 
  as.data.frame() %>%
  #filter(rownames(res0_8Heat) %in% uniqueHeat0_8downDEGs) %>%
  mutate(is_de=ifelse(rownames(res0_8Heat) %in% uniqueHeat0_8upDEGs, 1, 0),
         gene_id=rownames(.)) 

heat_up_df <- uniq_0_8Heat_upreg %>%
  left_join(features, by=c("gene_id"="ensembl_gene_id")) %>%
  mutate(gen_len=end_position-start_position)

heatUpGOtable8 <- GOenrichTable(heat_up_df)
nrow(heatUpGOtable8 %>% 
      mutate(term=forcats::fct_reorder(.f = term,.x = numInCat, .fun = sort)) %>% 
      filter(ontology=="BP", padj_over < 0.05))
# 15 terms total enriched

#### Timepoint 0 v 16 #####
# subset data
###### 0 v 16 no heat ####
countsTemp <- countsRaw[, which(names(countsRaw) %in% c("CN01", "CN02", "CN03", "CN161", "CN162", "CN163"))]
dim(countsTemp)
expdesignTemp <- expdesignRaw[which(names(countsRaw) %in% c("CN01", "CN02", "CN03", "CN161", "CN162", "CN163")),]

dds <- DESeqDataSetFromMatrix(countData=countsTemp, colData=expdesignTemp, design= ~ time)
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds  <- dds[keep,]
dds <- DESeq(dds)
resultsNames(dds)

res0_16nHeat <- results(dds, alpha = 0.05)
summary(res0_16nHeat)
sum(res0_16nHeat$padj < 0.05, na.rm=TRUE) #11125

###### 0 v 16 heat ####
countsTemp <- countsRaw[, which(names(countsRaw) %in% c("CN01", "CN02", "CN03", "CY161", "CY162", "CY163"))]
dim(countsTemp)
expdesignTemp <- expdesignRaw[which(names(countsRaw) %in% c("CN01", "CN02", "CN03", "CY161", "CY162", "CY163")),]

dds <- DESeqDataSetFromMatrix(countData=countsTemp, colData=expdesignTemp, design= ~ time)
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds  <- dds[keep,]
dds <- DESeq(dds)
resultsNames(dds)

res0_16Heat <- results(dds, alpha = 0.05)
summary(res0_16Heat)
sum(res0_16Heat$padj < 0.05, na.rm=TRUE) #7767 DEGs

###### Intersection no heat versus heat ####

DEGUpNoHeat <- as.data.frame(res0_16nHeat) %>%
  filter(padj < 0.05, log2FoldChange > 0) 
nrow(DEGUpNoHeat) #4477 DEGs

DEGUpHeat <- as.data.frame(res0_16Heat) %>%
  filter(padj < 0.05, log2FoldChange > 0) 
nrow(DEGUpHeat) #3234 DEGs

overlap1 <- calculate.overlap(x = list(rownames(DEGUpNoHeat), rownames(DEGUpHeat)))
str(overlap1)
##The genes that are only upregulated due to heat (not time) are of interest
# This is captured by not including the intersection genes
# overlap1[3] #list of intersecting genes
uniqueHeat0_16upDEGs <- unlist(overlap1[2])[!unlist(overlap1[2]) %in% unlist(overlap1[3])]


DEGDownNoHeat <- as.data.frame(res0_16nHeat) %>%
  filter(padj < 0.05, log2FoldChange < 0) 
nrow(DEGDownNoHeat) #6648 DEGs

DEGDownHeat <- as.data.frame(res0_16Heat) %>%
  filter(padj < 0.05, log2FoldChange < 0) 
nrow(DEGDownHeat) #4533 DEGs

overlap2 <- calculate.overlap(x = list(rownames(DEGDownNoHeat), rownames(DEGDownHeat)))
str(overlap2)
##The genes that are only downregulated due to heat (not time) are of interest
# This is captured by not including the intersection genes
# overlap2[3] #list of intersecting genes
uniqueHeat0_16downDEGs <- unlist(overlap2[2])[!unlist(overlap2[2]) %in% unlist(overlap2[3])]

######Down-regulated GO terms ####
## These are taken from the heat DEGs but only those that ARE NOT also DEGs without heat

# First we need to make a table that indicates which genes have differential expression of interest
# In this case, of interest are genes that are uniquely differentially expressed due to heat not time
# The "base" set are those that are expressed,

# This creates adds a column to the table indicating the genes of interest with a 1 
# and those not with a 0
uniq_0_16Heat_dnreg <- res0_16Heat %>% 
  as.data.frame() %>%
  mutate(is_de=ifelse(rownames(res0_16Heat) %in% uniqueHeat0_16downDEGs, 1, 0),
         gene_id=rownames(.)) 

# Now we merge this dataframe with the features dataframe that includes gene lengths and GO terms
# Genes can have multiple GO terms associated with them
heat_down_df <- uniq_0_16Heat_dnreg %>%
  left_join(features, by=c("gene_id"="ensembl_gene_id")) %>%
  mutate(gen_len=end_position-start_position)

# Complete GO enrichment and output a table
# Uses goseq approach
heatDownGOtable16 <- GOenrichTable(heat_down_df)

dim(heatDownGOtable16 %>% 
      mutate(term=forcats::fct_reorder(.f = term,.x = numInCat, .fun = sort)) %>% 
      filter(ontology=="BP", padj_over < 0.05))
# 66 terms total enriched

######Up-regulated GO terms ####
## These are taken from the heat DEGs but only those that ARE NOT also DEGs without heat
uniq_0_16Heat_upreg <- res0_16Heat %>% 
  as.data.frame() %>%
  mutate(is_de=ifelse(rownames(res0_16Heat) %in% uniqueHeat0_16upDEGs, 1, 0),
         gene_id=rownames(.)) 

heat_up_df <- uniq_0_16Heat_upreg %>%
  left_join(features, by=c("gene_id"="ensembl_gene_id")) %>%
  mutate(gen_len=end_position-start_position)

heatUpGOtable16 <- GOenrichTable(heat_up_df)
dim(heatUpGOtable16 %>% 
      mutate(term=forcats::fct_reorder(.f = term,.x = numInCat, .fun = sort)) %>% 
      filter(ontology=="BP", padj_over < 0.05))
# 8 terms total enriched

#### Timepoint 0 v 24 #####

# subset data
###### 0 v 24 no heat ####
countsTemp <- countsRaw[, which(names(countsRaw) %in% c("CN01", "CN02", "CN03", "CN241", "CN242", "CN243"))]
dim(countsTemp)
expdesignTemp <- expdesignRaw[which(names(countsRaw) %in% c("CN01", "CN02", "CN03", "CN241", "CN242", "CN243")),]
expdesignTemp

dds <- DESeqDataSetFromMatrix(countData=countsTemp, colData=expdesignTemp, design= ~ time)
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds  <- dds[keep,]
dds <- DESeq(dds)
resultsNames(dds)

res0_24nHeat <- results(dds, alpha = 0.05)
summary(res0_24nHeat)
sum(res0_24nHeat$padj < 0.05, na.rm=TRUE) #7416 due to 24 hours time

###### 0 v 24 heat ####
countsTemp <- countsRaw[, which(names(countsRaw) %in% c("CN01", "CN02", "CN03", "CY241", "CY242", "CY243"))]
dim(countsTemp)
expdesignTemp <- expdesignRaw[which(names(countsRaw) %in% c("CN01", "CN02", "CN03", "CY241", "CY242", "CY243")),]
expdesignTemp

dds <- DESeqDataSetFromMatrix(countData=countsTemp, colData=expdesignTemp, design= ~ time)
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds  <- dds[keep,]
dds <- DESeq(dds)
resultsNames(dds)

res0_24Heat <- results(dds, alpha = 0.05)
summary(res0_24Heat)
sum(res0_24Heat$padj < 0.05, na.rm=TRUE) #13514 DEGs

###### Intersection no heat versus heat ####

DEGUpNoHeat <- as.data.frame(res0_24nHeat) %>%
  filter(padj < 0.05, log2FoldChange > 0) 
nrow(DEGUpNoHeat) #2644 DEGs

DEGUpHeat <- as.data.frame(res0_24Heat) %>%
  filter(padj < 0.05, log2FoldChange > 0) 
nrow(DEGUpHeat) #5900

overlap1 <- calculate.overlap(x = list(rownames(DEGUpNoHeat), rownames(DEGUpHeat)))
str(overlap1)
##The genes that are only upregulated due to heat (not time) are of interest
# This is captured by not including the intersection genes
# overlap1[3] #list of intersecting genes
uniqueHeat0_24upDEGs <- unlist(overlap1[2])[!unlist(overlap1[2]) %in% unlist(overlap1[3])]


DEGDownNoHeat <- as.data.frame(res0_24nHeat) %>%
  filter(padj < 0.05, log2FoldChange < 0) 
nrow(DEGDownNoHeat) #4772 DEGs

DEGDownHeat <- as.data.frame(res0_24Heat) %>%
  filter(padj < 0.05, log2FoldChange < 0) 
nrow(DEGDownHeat) #7614 DEGs

overlap2 <- calculate.overlap(x = list(rownames(DEGDownNoHeat), rownames(DEGDownHeat)))
str(overlap2)
##The genes that are only downregulated due to heat (not time) are of interest
# This is captured by not including the intersection genes
# overlap2[3] #list of intersecting genes
uniqueHeat0_24downDEGs <- unlist(overlap2[2])[!unlist(overlap2[2]) %in% unlist(overlap2[3])]

####GO terms 0 v 24 unique to heat #####

######Down-regulated GO terms ####
## These are taken from the heat DEGs but only those that ARE NOT also DEGs without heat

# First we need to make a table that indicates which genes have differential expression of interest
# In this case, of interest are genes that are uniquely differentially expressed due to heat not time
# The "base" set are those that are expressed,

# This creates adds a column to the table indicating the genes of interest with a 1 
# and those not with a 0
uniq_0_24Heat_dnreg <- res0_24Heat %>% 
  as.data.frame() %>%
  mutate(is_de=ifelse(rownames(res0_24Heat) %in% uniqueHeat0_24downDEGs, 1, 0),
         gene_id=rownames(.)) 

# Now we merge this dataframe with the features dataframe that includes gene lengths and GO terms
# Genes can have multiple GO terms associated with them
heat_down_df <- uniq_0_24Heat_dnreg %>%
  left_join(features, by=c("gene_id"="ensembl_gene_id")) %>%
  mutate(gen_len=end_position-start_position)

# Complete GO enrichment and output a table
# Uses goseq approach
heatDownGOtable24 <- GOenrichTable(heat_down_df)
dim(heatDownGOtable24 %>% 
      mutate(term=forcats::fct_reorder(.f = term,.x = numInCat, .fun = sort)) %>% 
      filter(ontology=="BP", padj_over < 0.05))
# 48 terms total enriched

######Up-regulated GO terms ####
## These are taken from the heat DEGs but only those that ARE NOT also DEGs without heat
uniq_0_24Heat_upreg <- res0_24Heat %>% 
  as.data.frame() %>%
  mutate(is_de=ifelse(rownames(res0_24Heat) %in% uniqueHeat0_24upDEGs, 1, 0),
         gene_id=rownames(.)) 

heat_up_df <- uniq_0_24Heat_upreg %>%
  left_join(features, by=c("gene_id"="ensembl_gene_id")) %>%
  mutate(gen_len=end_position-start_position)

heatUpGOtable24 <- GOenrichTable(heat_up_df)
dim(heatUpGOtable24 %>% 
      mutate(term=forcats::fct_reorder(.f = term,.x = numInCat, .fun = sort)) %>% 
      filter(ontology=="BP", padj_over < 0.05))
# 4 terms total enriched


heatUpGOtable24BP <- heatUpGOtable24 %>% 
  mutate(term=forcats::fct_reorder(.f = term,.x = numInCat, .fun = sort)) %>% 
  filter(ontology=="BP", padj_over < 0.05)

heatDownGOtable24BP <- heatDownGOtable24 %>% 
  mutate(term=forcats::fct_reorder(.f = term,.x = numInCat, .fun = sort)) %>% 
  filter(ontology=="BP", padj_over < 0.05)

heatUpGOtable16BP <- heatUpGOtable16 %>% 
  mutate(term=forcats::fct_reorder(.f = term,.x = numInCat, .fun = sort)) %>% 
  filter(ontology=="BP", padj_over < 0.05)

heatDownGOtable16BP <- heatDownGOtable16 %>% 
  mutate(term=forcats::fct_reorder(.f = term,.x = numInCat, .fun = sort)) %>% 
  filter(ontology=="BP", padj_over < 0.05)

heatUpGOtable8BP <- heatUpGOtable8 %>% 
  mutate(term=forcats::fct_reorder(.f = term,.x = numInCat, .fun = sort)) %>% 
  filter(ontology=="BP", padj_over < 0.05)

heatDownGOtable8BP <- heatDownGOtable8 %>% 
  mutate(term=forcats::fct_reorder(.f = term,.x = numInCat, .fun = sort)) %>% 
  filter(ontology=="BP", padj_over < 0.05)

#### Writing out GO enrich tables and DEGs #####
write.table(heatUpGOtable8BP, "hr0_8GOenrichBPUpDEGs.csv", quote = F, row.names = F, sep = "\t")
write.table(heatDownGOtable8BP, "hr0_8GOenrichBPDownDEGs.csv", quote = F, row.names = F, sep = "\t")

write.table(heatUpGOtable16BP, "hr0_16GOenrichBPUpDEGs.csv", quote = F, row.names = F, sep = "\t")
write.table(heatDownGOtable16BP, "hr0_16GOenrichBPDownDEGs.csv", quote = F, row.names = F, sep = "\t")

write.table(heatUpGOtable24BP, "hr0_24GOenrichBPUpDEGs.csv", quote = F, row.names = F, sep = "\t")
write.table(heatDownGOtable24BP, "hr0_24GOenrichBPDownDEGs.csv", quote = F, row.names = F, sep = "\t")


write.csv(as.data.frame(res0_8Heat)[rownames(as.data.frame(res0_8Heat)) %in% c(uniqueHeat0_8downDEGs, uniqueHeat0_8upDEGs),], "hr0_8DEGs.csv", quote = F)

write.csv(as.data.frame(res0_16Heat)[rownames(as.data.frame(res0_16Heat)) %in% c(uniqueHeat0_16downDEGs, uniqueHeat0_16upDEGs),], "hr0_16DEGs.csv", quote = F)

write.csv(as.data.frame(res0_24Heat)[rownames(as.data.frame(res0_24Heat)) %in% c(uniqueHeat0_24downDEGs, uniqueHeat0_24upDEGs),], "hr0_24DEGs.csv", quote = F)


####Visualizing GO term comparisons #####

heatDownGOtable8 <- heatDownGOtable8 %>% 
  filter(ontology=="BP", padj_over < 0.05)
heatDownGOtable16 <- heatDownGOtable16  %>% 
  filter(ontology=="BP", padj_over < 0.05)
heatDownGOtable24 <- heatDownGOtable24  %>% 
  filter(ontology=="BP", padj_over < 0.05)

lt <- list(hr8 = heatDownGOtable8, hr16 = heatDownGOtable16, hr24 = heatDownGOtable24)
BiocManager::install("org.At.tair.db")

pdf( file="simplifyMultiGODownRegkmeans05.pdf", width=10, height=8 ) # width and height are in inches
simplifyGOFromMultipleLists(lt, padj_cutoff = 0.05, db = "org.At.tair.db", exclude_words = c("ii"),
                            min_term = 5, 
                            padj_column = "padj_over",
                            show_barplot = TRUE,
                            method = "kmeans", max_words = 20,
                            word_cloud_grob_param = list(max_width = 80),
                            fontsize_range = c(7,15),
                            #heatmap_param = list(col = c("grey", "orange"))
                            )
dev.off()