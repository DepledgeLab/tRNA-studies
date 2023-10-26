library(ggplot2)
library(dplyr)
library(DESeq2)
library(EnhancedVolcano)

setwd("set/working/directory")

### Read in count files
tngidA_1 <- read.csv('./counts/Haussler_tngidA_total_tRNA_SIM-24h-1-merged.counts.txt', sep='\t', header = FALSE)
tngidA_2 <- read.csv('./counts/Haussler_tngidA_total_tRNA_SIM-24h-2-merged.counts.txt', sep='\t', header = FALSE)
tnladS_1 <- read.csv('./counts/Haussler_tnladS_total_tRNA_SIM-24h-1-merged.counts.txt', sep='\t', header = FALSE)
tnladS_2 <- read.csv('./counts/Haussler_tnladS_total_tRNA_SIM-24h-2-merged.counts.txt', sep='\t', header = FALSE)

### Set columns names
colnames(tngidA_1)<-c("tngidA_1","tRNA")
colnames(tngidA_2)<-c("tngidA_2","tRNA")
colnames(tnladS_1)<-c("tnladS_1","tRNA")
colnames(tnladS_2)<-c("tnladS_2","tRNA")

### merge dataframes into single dataframe & reorder columns
counts <- full_join(tngidA_1, tngidA_2, by = "tRNA") %>%
  full_join(tnladS_1, by = "tRNA") %>%
  full_join(tnladS_2, by = "tRNA")
col_order <- c("tRNA", "tngidA_1", "tngidA_2", "tnladS_1", "tnladS_2")
counts_fixed <- counts[, col_order]

### Drop SeC_UCA
counts_fixedAlt <- subset(counts_fixed, !grepl("SeC_UCA", counts_fixed[,1]))
countTable <- counts_fixedAlt
rownames(countTable) <- countTable$tRNA
countTable <- countTable[, -1]

### Read in experimental setup
samples <- read.table("samples.txt", header=TRUE)

### Make DeSeq dataset
Dataset <- DESeqDataSetFromMatrix(countData = countTable, colData=samples, design=~condition)

# Run DESEQ and generate a simple plot showing the distribution of regulated and unregulated genes
Dataset$condition <- relevel(Dataset$condition, "tnladS")
DatasetProcessed <- DESeq(Dataset)

### Perform contrast analyses to produce lists of differentially regulated genes between conditions (pairwise)
result <- lfcShrink(DatasetProcessed, contrast=c("condition","tngidA","tnladS"), type="ashr")
baseMean_tngidA = rowMeans(counts(DatasetProcessed,normalized=TRUE)[,DatasetProcessed$condition == "tngidA"])
baseMean_tnladS = rowMeans(counts(DatasetProcessed,normalized=TRUE)[,DatasetProcessed$condition == "tnladS"])
result = cbind(as.data.frame(result), baseMean_tngidA, baseMean_tnladS)
write.csv(result, "deseq2.output.csv", row.names=TRUE)

# Plotting
volc<-EnhancedVolcano(toptable = result, x = 'log2FoldChange', y = 'padj', 
                colAlpha = 1, 
                gridlines.major = TRUE, 
                gridlines.minor = FALSE, 
                lab = NA, 
                pointSize = 2, 
                FCcutoff = 0.2, 
                pCutoff = 0.05,
                title = 'tngidA vs. tnladS',
                hline = 0.05,
                xlim = c(-1,1),
                col=c('black', 'indianred4', 'indianred4', 'firebrick1'))

pdf(file = "volcano.pdf", width = 5, height = 7)
volc
dev.off()


