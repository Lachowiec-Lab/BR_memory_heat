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
library(tidyr)
library(forcats)

### Functions ####
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

quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
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

###Time 16 mock heat v rest ####
##only examining time 16, but finding genes that are DEGs only in mock heat from all others
expdesignRaw
expdesignH <- expdesignRaw %>%
  as.data.frame() %>%
  mutate(heatOnly = ifelse(stress == "heat" & trt == "ck", "heat", "else"))
expdesignH

time16 <- 16
selected_columns1 <- grep(time16, colnames(countsRaw), value = F)
time16 <- countsRaw[,c(selected_columns1)]
time16_design <- expdesignH[c(selected_columns1),]

ncol(time16) == nrow(time16_design)

dds <- DESeqDataSetFromMatrix(countData = time16, colData = time16_design, design= ~ heatOnly)
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds  <- dds[keep,]
dds <- DESeq(dds)
resultsNames(dds)

res16alt <- results(dds, name="heatOnly_heat_vs_else", alpha = 0.05)
summary(res16alt)

t16BRdegs <-  as.data.frame(res16alt)[as.data.frame(res16alt)$padj < 0.05,]
t16BRdegs <- t16BRdegs[complete.cases(t16BRdegs), ]

write.csv(t16BRdegs, "hr16_BReffectDEGs.csv", quote = F)

#### Heat response genes ####

hot <- features %>%
  filter(str_detect(name_1006, "response to heat"))
cellhot <- features %>%
  filter(str_detect(name_1006, "cellular response to heat")) 
acclim <- features %>%
  filter(str_detect(name_1006, "heat acclimation")) 
hs <- features %>%
  filter(str_detect(description, "Heat shock"))

acqThermo <- c("TraesCS7D02G440200","TraesCS7A02G451100","TraesCS7B02G351000","TraesCS3B02G240000","TraesCS3D02G212100","TraesCS3A02G209200","TraesCS1B02G475700","TraesCS1A02G441200","TraesCS1D02G449700","TraesCS3D02G212100","TraesCS3B02G240000","TraesCS3A02G209200","TraesCS1B02G475700","TraesCS1A02G441200","TraesCS1D02G449700","TraesCS1B02G107200","TraesCS1D02G089500","TraesCS1A02G088200","TraesCS2B02G096200","TraesCS2D02G080000","TraesCS2A02G082100","TraesCS4B02G197900","TraesCS4A02G106300","TraesCS4D02G198500","TraesCS4D02G212500","TraesCS4A02G092600","TraesCS4A02G092700","TraesCS5B02G110600","TraesCS5D02G124200","TraesCS5A02G109400","TraesCS2D02G426400","TraesCS2A02G428300","TraesCS7A02G410900","TraesCS7D02G403900","TraesCS7B02G310100","TraesCS6A02G034200","TraesCS6B02G048200","TraesCS6D02G039800","TraesCS1B02G118900","TraesCS1D02G099600","TraesCS1A02G091100","TraesCS3A02G069900","TraesCS3D02G068900","TraesCS3B02G083400","TraesCS7B02G455500","TraesCS7A02G537000","TraesCS7D02G524700","TraesCS4D02G331100","TraesCS4B02G335500","TraesCS5A02G505300","TraesCS1A02G205300","TraesCS1B02G218900","TraesCS1D02G208600")

hotAll <- c(hot$ensembl_gene_id, cellhot$ensembl_gene_id, acclim$ensembl_gene_id, acqThermo, hs$ensembl_gene_id) #537
length(unique(hotAll)) # 213 total unique

#remove Hsfs since already analyzed
hsfs <- c("TraesCS5A02G437900", "TraesCS5B02G440700", "TraesCS5D02G445100",
          "TraesCS2A02G089300", "TraesCS2B02G105100", "TraesCS2D02G087400",
          "TraesCS4A02G027700", "TraesCS4B02G278100", "TraesCS4D02G276500",
          "TraesCS5A02G383800", "TraesCS5D02G393200", "TraesCS1A02G375600",
          "TraesCS1B02G396000", "TraesCS1D02G382900", "TraesCS3A02G405200",
          "TraesCS3B02G438900")

hotAll <- setdiff(hotAll, hsfs)
length(hotAll) #208

res16alt_hot <- res16alt[rownames(res16alt) %in% unique(hotAll),]
#res16alt_hot <- res16alt_hot[-which(is.na(res16alt_hot$padj)),]
dim(res16alt_hot) #185 included of the 223 hot found
nrow(res16alt_hot[res16alt_hot$padj< 0.05,]) #14
nrow(res16alt_hot[res16alt_hot$padj< 0.05,])/nrow(res16alt_hot)

#####Figure 5A ####
multigene2graph(countsAll, rownames(res16alt_hot[res16alt_hot$padj< 0.05,]), 3)

#### Histones ####
histones <- features %>%
  filter(go_id == "GO:0000786") %>%
  filter(str_detect(description, "Histone|histone|H2A.5|H2A.7")) # should be thorough

length(unique(histones$ensembl_gene_id)) # no gene id's are repeated
res16alt_histones <- res16alt[rownames(res16alt) %in% histones$ensembl_gene_id,]
res16alt_histones <- res16alt_histones[-which(is.na(res16alt_histones$padj)),]
nrow(res16alt_histones) #425 included of the 456 histones found

nrow(res16alt_histones[res16alt_histones$padj < 0.05,]) # 308
nrow(res16alt_histones[res16alt_histones$log2FoldChange < 0 & res16alt_histones$padj < 0.05,]) # all 308 are downregulated

diffExpHistones <- res16alt_histones[res16alt_histones$padj < 0.05,] 
meanNormCounts <- countsMean[rownames(countsMean) %in% rownames(diffExpHistones),]
time16 <- "1"
selected_columns1 <- grep(time16, colnames(meanNormCounts), value = F)
meanNormCounts16 <- meanNormCounts[,c(selected_columns1)]

mat_breaks <- quantile_breaks(as.matrix(meanNormCounts16), n = 10)

histoneGroups <- unique(features[features$ensembl_gene_id %in% rownames(diffExpHistones), c("ensembl_gene_id","description")])
rownames(histoneGroups) <- histoneGroups$ensembl_gene_id
histoneGroups <- histoneGroups[-1]
histoneGroups$description <- gsub("\\[[^][]*]", "", histoneGroups$description)
names(histoneGroups) <- "Histone protein"

ann_colors = list(
  `Histone protein` = c(`Centromeric histone 3 ` = "red", `Histone H1 WH1A.2 ` = "yellow", 
                        `Histone H2A ` = "orange",  `Histone H2A.1 ` = "darkorange", `Protein H2A.5 ` = "darkorange2", `Protein H2A.7 ` = "darkorange4", 
                        `Histone H2B ` = "green", `Histone H2B.1 ` = "darkgreen", `Histone H2B.4 ` = "darkolivegreen",`Histone H2B.5 ` = "darkolivegreen1",
                        `Histone H3 ` = "blue", `Histone H3.2 ` = "lightblue",
                        `Histone H4 ` = "violet", `Histone H4 variant TH011 ` = "purple", `Histone H4 variant TH091 ` = "darkorchid4"
  ))


#####Figure 5B ####
pheatmap(meanNormCounts16, color = brewer.pal(9, "Reds"),
         breaks = mat_breaks, 
         show_rownames = F,
         annotation_row = histoneGroups, 
         annotation_colors = ann_colors,
         cluster_cols = F
)