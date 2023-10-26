library(ggplot2)
library(dplyr)
library(patchwork)

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

sum_tngidA_1 <- sum(tngidA_1$tngidA_1)
sum_tngidA_2 <- sum(tngidA_2$tngidA_2)
sum_tnladS_1 <- sum(tnladS_1$tnladS_1)
sum_tnladS_2 <- sum(tnladS_2$tnladS_2)

norm_tngidA_1 <- tngidA_1 %>%
  mutate(tngidA_1_tpm = (tngidA_1 / sum_tngidA_1) * 1e6)

norm_tngidA_2 <- tngidA_2 %>%
  mutate(tngidA_2_tpm = (tngidA_2 / sum_tngidA_2) * 1e6)

norm_tnladS_1 <- tnladS_1 %>%
  mutate(tnladS_1_tpm = (tnladS_1 / sum_tnladS_1) * 1e6)

norm_tnladS_2 <- tnladS_2 %>%
  mutate(tnladS_2_tpm = (tnladS_2 / sum_tnladS_2) * 1e6)


norm_tngidA_1$id <- sub("_.*$|-GCU|-UUG", "", norm_tngidA_1$tRNA)
norm_tngidA_2$id <- sub("_.*$|-GCU|-UUG", "", norm_tngidA_2$tRNA)
norm_tnladS_1$id <- sub("_.*$|-GCU|-UUG", "", norm_tnladS_1$tRNA)
norm_tnladS_2$id <- sub("_.*$|-GCU|-UUG", "", norm_tnladS_2$tRNA)

### merge dataframes into single dataframe
counts <- full_join(norm_tngidA_1, norm_tngidA_2, by = "tRNA") %>%
  full_join(norm_tnladS_1, by = "tRNA") %>%
  full_join(norm_tnladS_2, by = "tRNA")

columns_to_drop <- c("id.y", "id.x.x", "id.y.y")
counts <- counts[, -which(names(counts) %in% columns_to_drop)]
names(counts)[names(counts) == "id.x"] <- "id"

### Reorder columns
col_order <- c("tRNA", "id", "tngidA_1", "tngidA_2", "tnladS_1", "tnladS_2", "tngidA_1_tpm", "tngidA_2_tpm", "tnladS_1_tpm", "tnladS_2_tpm")
counts_fixed <- counts[, col_order]
#counts_fixed$id <- as.factor(counts_fixed$id)
counts_fixed$id <- sub("GlT", "Glu", counts_fixed$id)

# Calculate Pearson correlation coefficient
tngidA_corr <- cor(counts_fixed$tngidA_1, counts_fixed$tngidA_2, method = "pearson")
tnladS_corr <- cor(counts_fixed$tnladS_1, counts_fixed$tnladS_2, method = "pearson")
comp_biorep1_corr <- cor(counts_fixed$tngidA_1, counts_fixed$tnladS_1, method = "pearson")
comp_biorep2_corr <- cor(counts_fixed$tngidA_2, counts_fixed$tnladS_2, method = "pearson")

# Create a scatter plot
tngidA_bioreps <- ggplot(counts_fixed, aes(x = tngidA_1_tpm, y = tngidA_2_tpm, color = id)) + geom_point(size = 2) + theme_classic() + scale_y_continuous(breaks = seq(0, 140000, by = 20000), limits=c(0,140000)) + scale_x_continuous(breaks = seq(0, 140000, by = 20000), limits=c(0,140000)) + geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") + 
      labs(x = "tngidA biorep1 normalized read count", y = "tngidA biorep2 normalized read count") + guides(color = "none")

tnladS_bioreps <- ggplot(counts_fixed, aes(x = tnladS_1_tpm, y = tnladS_2_tpm, color=id)) + geom_point(size = 2) + theme_classic() + scale_y_continuous(breaks = seq(0, 140000, by = 20000), limits=c(0,140000)) + scale_x_continuous(breaks = seq(0, 140000, by = 20000), limits=c(0,140000)) + geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  labs(x = "tnladS biorep1 normalized read count", y = "tnladS biorep2 normalized read count") + guides(color = "none")

comparison_biorep1 <- ggplot(counts_fixed, aes(x = tngidA_1_tpm, y = tnladS_2_tpm, color=id)) + geom_point(size = 2) + theme_classic() + scale_y_continuous(breaks = seq(0, 140000, by = 20000), limits=c(0,140000)) + scale_x_continuous(breaks = seq(0, 140000, by = 20000), limits=c(0,140000)) + geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") + 
  labs(x = "tngidA biorep1 normalized read count", y = "tnladS biorep1 normalized read count") + guides(color = "none")

comparison_biorep2 <- ggplot(counts_fixed, aes(x = tngidA_1_tpm, y = tnladS_2_tpm, color=id)) + geom_point(size = 2) + theme_classic() + scale_y_continuous(breaks = seq(0, 140000, by = 20000), limits=c(0,140000)) + scale_x_continuous(breaks = seq(0, 140000, by = 20000), limits=c(0,140000)) + geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") + 
  labs(x = "tngidA biorep2 normalized read count", y = "tnladS biorep2 normalized read count") + guides(color = "none")

placeholder <- ggplot(counts_fixed, aes(x = tngidA_1_tpm, y = tnladS_2_tpm, color=id)) + geom_point(size = 2) + theme_classic() + scale_y_continuous(breaks = seq(0, 140000, by = 20000), limits=c(0,140000)) + scale_x_continuous(breaks = seq(0, 140000, by = 20000), limits=c(0,140000)) + geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") + 
  labs(x = "tngidA biorep2 normalized read count", y = "tnladS biorep2 normalized read count")

final_plot <- tngidA_bioreps + tnladS_bioreps + comparison_biorep1 + comparison_biorep2 + placeholder + plot_layout(ncol = 2)
final_plot

pdf(file = "test.pdf", width = 13, height = 15)
final_plot
dev.off()


