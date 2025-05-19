### Libraries ####

library(dplyr)
library(DESeq2)
library(RColorBrewer)
library(viridis)
library(pheatmap)
library(goseq)
library(ggplot2)
library(egg)
library(patchwork)

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
multigene2graph <- function(data, genevector, coluNum){
  temp <- data[rownames(data) %in% genevector,]
  samples <- unique(str_sub(colnames(data), 1,3))
  forLineGraph <- data.frame(fullTreat = samples)
  forLineGraph$treatment <- str_sub(forLineGraph$fullTreat, 1,1)
  forLineGraph$stress <- str_sub(forLineGraph$fullTreat, 2,2)
  forLineGraph$time <- str_sub(forLineGraph$fullTreat, 3,3)
  forLineGraph$time <- factor(forLineGraph$time, levels = c("0", "8", "1", "2"))
  forLineGraph$BRHeat <- str_sub(forLineGraph$fullTreat, 1,2)
  ## Calculate means for each treatment
  medCounts <- matrix(data = NA, nrow = nrow(temp), ncol = length(samples))
  colnames(medCounts) <- samples
  for (i in 1:length(samples)){
    pattern <- samples[i]
    medCounts[,i] <- rowMeans(temp[,grep(pattern, colnames(temp))])
    print(i)
  }
  ## Calculate standard deviation for each treatment
  seCounts <- matrix(data = NA, nrow = nrow(temp), ncol = length(samples))
  colnames(seCounts) <- samples
  for (i in 1:length(samples)){
    pattern <- samples[i]
    seCounts[,i] <- matrixStats::rowSds(temp[,grep(pattern, colnames(temp))], na.rm = T)/sqrt(ncol(temp[,grep(pattern, colnames(temp))]))
    print(i)
  }
  
  temp2 <- cbind(forLineGraph, t(medCounts))
  colnames(temp2)[6:(5+length(genevector))] <- rownames(temp)
  temp3 <- temp2 %>%
    tidyr::pivot_longer(cols = starts_with("Traes"), names_to = "gene", values_to = "counts") %>%
    mutate(stage=forcats::fct_relevel(time,c("0", "8", "1", "2"))) %>%
    mutate(se = as.vector(seCounts))
  return(ggplot(temp3, aes(x = time, y = counts, group = BRHeat, col = stress)) +
           geom_point(position=position_dodge(0.1)) +
           geom_line(aes(linetype=treatment)) +
           geom_errorbar(aes(ymin=counts-se, ymax=counts+se, linetype=treatment), 
                         width=.4, position=position_dodge(0.1)) +
           scale_linetype(labels = c("Mock", "BR seed soak")) +
           facet_wrap(~gene, ncol = coluNum, scales = "free") +
           scale_color_manual(values = c("#478C0E", "#E1BE6A"), labels = c("22\u00b0", "29\u00b0")) +
           scale_x_discrete(labels=c("0" = "0", "8" = "8",
                                     "1" = "16", "2" = "24")) +
           xlab(label = "Time (hours after heat)") +
           ylab(label = expression(paste("log"[2],"(Mean normalized counts)"))) +
           egg::theme_article() + 
           theme(axis.text.x = element_text(angle=45, hjust = 1)))
  
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

### Effect of BR, no heat #####
#### DESeq2 #####

## Data subsetting
## Using only time point 0

ncol(countsRaw) == nrow(expdesignRaw)
expdesignRaw
time0 <- "0"

selected_columns1 <- grep(time0, colnames(countsRaw), value = F)
time0 <- countsRaw[,c(selected_columns1)]
time0_design <- expdesignRaw[c(selected_columns1),]

ncol(time0) == nrow(time0_design)

dds <- DESeqDataSetFromMatrix(countData=time0, colData=time0_design, design= ~ trt)
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds  <- dds[keep,]
dds <- DESeq(dds)
resultsNames(dds)

resBR <- results(dds, name="trt_ss_vs_ck", alpha = 0.05)
summary(resBR)
sum(resBR$padj < 0.05, na.rm=TRUE) #3520 DEGs due to BR in timepoint 1


t0BRdegs <-  as.data.frame(resBR)[as.data.frame(resBR)$padj < 0.05,]
t0BRdegs <- t0BRdegs[complete.cases(t0BRdegs), ]

write.csv(t0BRdegs, "hr0_BReffectDEGs.csv", quote = F)

#### GO term enrichments ####
res_de_tab_upreg <- resBR%>% 
  as.data.frame() %>%
  filter(log2FoldChange > 0) %>%
  mutate(is_de=ifelse(padj >= 0.05|is.na(padj), 0, 1),
         gene_id=rownames(.)) 

res_de_tab_dnreg <- resBR %>% 
  as.data.frame() %>%
  filter(log2FoldChange<0) %>%
  mutate(is_de=ifelse(padj>=0.05|is.na(padj), 0, 1),
         gene_id=rownames(.)) 

br_up_df <- res_de_tab_upreg %>%
  left_join(features, by=c("gene_id"="ensembl_gene_id")) %>%
  mutate(gen_len=end_position-start_position)

br_down_df <- res_de_tab_dnreg %>%
  left_join(features, by=c("gene_id"="ensembl_gene_id")) %>%
  mutate(gen_len=end_position-start_position)

brUpGOtable <- GOenrichTable(br_up_df)
brDownGOtable <- GOenrichTable(br_down_df)

brUpGOtable$GOterm <- stringr::str_wrap(brUpGOtable$term, width = 33)
go_brUpBP <- brUpGOtable %>%
  mutate(term=forcats::fct_reorder(.f = GOterm, .x = -padj_over, .fun = sort)) %>%
  filter(ontology=="BP", padj_over < 0.01, numDEInCat > 5) %>%
  slice_min(order_by = padj_over, n = 10) %>%
  ggplot(aes(x=term, y=-log10(padj_over), size=numDEInCat, color=numInCat)) +
  geom_point() +
  lims(color = c(3,1500), size = c(0,150)) +
  coord_flip() +
  labs(color = "Number of genes", size = "Number of DEGs") +
  ylab(label = expression(-log[10]*"(p-adjusted)")) +
  xlab(label = "GO term") +
  egg::theme_article(base_size = 15)  
go_brUpBP

brDownGOtable$GOterm <- stringr::str_wrap(brDownGOtable$term, width = 33)
go_brDownBP <- brDownGOtable %>%
  mutate(term=forcats::fct_reorder(.f = GOterm,.x = -padj_over, .fun = sort)) %>%
  filter(ontology=="BP", padj_over < 0.01, numDEInCat > 5) %>%
  slice_min(order_by = padj_over, n = 10) %>%
  ggplot(aes(x=term, y=-log10(padj_over), size=numDEInCat, color=numInCat)) +
  geom_point() +
  lims(color = c(3,1500), size = c(0,150)) +
  coord_flip() +
  labs(color = "Number of genes", size = "Number of DEGs") +
  ylab(label = expression(-log[10]*"(p-adjusted)")) +
  xlab(label = "GO term") +
  egg::theme_article(base_size = 15) 
go_brDownBP

#### Histone expression ####
## Examining Histone Genes ###

histones <- features %>%
  filter(go_id == "GO:0000786") %>%
  filter(str_detect(description, "Histone|histone|H2A.5|H2A.7")) # should be thorough

resBR_histones <- resBR[rownames(resBR) %in% histones$ensembl_gene_id,]
resBR_histones <- resBR_histones[-which(is.na(resBR_histones$padj)),]
dim(resBR_histones) #424 included of the 456 histones found

nrow(resBR_histones[resBR_histones$padj < 0.05,]) # 60 show differential expression
diffExpHistones <- resBR_histones[resBR_histones$padj < 0.05,] 

t0_normCounts <- countsMean[rownames(countsMean) %in% rownames(diffExpHistones),str_detect(colnames(countsMean), "0")]
head(t0_normCounts)

histoneGroups <- unique(features[features$ensembl_gene_id %in% rownames(diffExpHistones), c("ensembl_gene_id","description")])
rownames(histoneGroups) <- histoneGroups$ensembl_gene_id
histoneGroups <- histoneGroups[-1]
histoneGroups$description <- gsub("\\[[^][]*]", "", histoneGroups$description)
names(histoneGroups) <- "Histone protein"
table(histoneGroups$`Histone protein`)

temp0 <- merge(t0_normCounts, histoneGroups, by = 0, all.x = T)
temp1 <- temp0 %>%
  mutate(gene = Row.names) %>%
  pivot_longer(mock_preHeat_time0:seed_preHeat_time0, names_to = "treatment", values_to = "meanNormCounts")
head(temp1)

### These are the colors used for histones and variants
ann_colors = list(
  `Histone protein` = c(`Centromeric histone 3 ` = "red", `Histone H1 WH1A.2 ` = "yellow", 
                        `Histone H2A ` = "orange",
                        `Histone H2B ` = "green", `Histone H2B.1 ` = "darkgreen",
                        `Histone H3 ` = "blue", `Histone H3.2 ` = "lightblue",
                        `Histone H4 ` = "violet", `Histone H4 variant TH011 ` = "purple"
  ))
## adjusting the breaks for the heatmap color scale
mat_breaks <- quantile_breaks(as.matrix(t0_normCounts), n = 10)
br0heatmap <- pheatmap(t0_normCounts, color = brewer.pal(9, "Reds"),
         breaks = mat_breaks, 
         show_rownames = F,
         annotation_row = histoneGroups, 
         annotation_colors = ann_colors,
         cluster_cols = F,
         border_color = "NA"
)


### Figure 4 panels####
go_brUpBP + go_brDownBP + plot_layout(guides = 'collect', axes = 'collect') + plot_annotation(tag_levels = 'A')
br0heatmap

### Basal and Acquired thermotolerance ####

#### Thermotolerance genes from GO terms ####
hot <- features %>%
  filter(str_detect(name_1006, "response to heat"))
table(hot$description)

cellhot <- features %>%
  filter(str_detect(name_1006, "cellular response to heat")) 
table(cellhot$description)

acclim <- features %>%
  filter(str_detect(name_1006, "heat acclimation")) 
table(acclim$description)


hot <- rbind(hot, cellhot, acclim)
length(unique(hot$ensembl_gene_id)) # no gene id's are repeated. 153 total
resBR_hot <- resBR[rownames(resBR) %in% hot$ensembl_gene_id,]
resBR_hot <- resBR_hot[-which(is.na(resBR_hot$padj)),]
nrow(resBR_hot) #118 included of the 153 hot found
nrow(resBR_hot[resBR_hot$padj < 0.05,])
resBR_hot[resBR_hot$padj < 0.05,]

#### Manually annotated acquired thermotolerance genes ####
# This is the list of homologs of A thaliana thermotolerance genes identified via BLAST
acqThermo <- acqThermo <- c("TraesCS7D02G440200","TraesCS7A02G451100","TraesCS7B02G351000",
                            "TraesCS3B02G240000","TraesCS3D02G212100","TraesCS3A02G209200",
                            "TraesCS1B02G475700","TraesCS1A02G441200","TraesCS1D02G449700",
                            "TraesCS3D02G212100","TraesCS3B02G240000","TraesCS3A02G209200",
                            "TraesCS1B02G475700","TraesCS1A02G441200","TraesCS1D02G449700",
                            "TraesCS1B02G107200","TraesCS1D02G089500","TraesCS1A02G088200",
                            "TraesCS2B02G096200","TraesCS2D02G080000","TraesCS2A02G082100",
                            "TraesCS4B02G197900","TraesCS4A02G106300","TraesCS4D02G198500",
                            "TraesCS4D02G212500","TraesCS4A02G092600","TraesCS4A02G092700",
                            "TraesCS5B02G110600","TraesCS5D02G124200","TraesCS5A02G109400",
                            "TraesCS2D02G426400","TraesCS2A02G428300","TraesCS7A02G410900",
                            "TraesCS7D02G403900","TraesCS7B02G310100","TraesCS6A02G034200",
                            "TraesCS6B02G048200","TraesCS6D02G039800","TraesCS1B02G118900",
                            "TraesCS1D02G099600","TraesCS1A02G091100","TraesCS3A02G069900",
                            "TraesCS3D02G068900","TraesCS3B02G083400","TraesCS7B02G455500",
                            "TraesCS7A02G537000","TraesCS7D02G524700","TraesCS4D02G331100",
                            "TraesCS4B02G335500","TraesCS5A02G505300","TraesCS1A02G205300",
                            "TraesCS1B02G218900","TraesCS1D02G208600")

acqThermoDf <- features %>%
  filter(ensembl_gene_id %in% acqThermo) # should be thorough
length(unique(acqThermoDf$ensembl_gene_id)) 

resBR_acqThermoDf <- resBR[rownames(resBR) %in% acqThermoDf$ensembl_gene_id,]
nrow(resBR_acqThermoDf) 
nrow(resBR_acqThermoDf[resBR_acqThermoDf$padj< 0.05,])
resBR_acqThermoDf[resBR_acqThermoDf$padj< 0.05,]

multigene2graph(countsAll, rownames(resBR_acqThermoDf[resBR_acqThermoDf$padj< 0.05,]), 4)

# HsfA2 genes
hsfs <- c("TraesCS5A02G437900", "TraesCS5B02G440700", "TraesCS5D02G445100",
          "TraesCS2A02G089300", "TraesCS2B02G105100", "TraesCS2D02G087400",
          "TraesCS4A02G027700", "TraesCS4B02G278100", "TraesCS4D02G276500",
          "TraesCS5A02G383800", "TraesCS5D02G393200", "TraesCS1A02G375600",
          "TraesCS1B02G396000", "TraesCS1D02G382900", "TraesCS3A02G405200",
          "TraesCS3B02G438900")

# are any hsfA2s expressed or DEGs?
resBR_HSFsDf <- resBR[rownames(resBR) %in% hsfs,]
resBR_HSFsDf <- resBR_HSFsDf[-which(is.na(resBR_HSFsDf$padj)),]
nrow(resBR_HSFsDf)
nrow(resBR_HSFsDf[resBR_HSFsDf$padj < 0.05,]) # 0 show differential expression

####Figure S4 ####
# Patterns of HsfA2 expression
multigene2graph(countsAll, hsfs, coluNum = 4)