#!/usr/bin/env Rscript
# =============================================================================
# Differential Expression Analysis with DESeq2
# =============================================================================
# Author: Waqas Yousaf
# Date: May 2026
# Usage: Rscript 02_rnaseq_deseq2.R
# =============================================================================

# Load libraries
library(DESeq2)
library(tidyverse)

# Set paths
COUNT_MATRIX <- "../04_example_data/example_count_matrix.csv"
METADATA <- "../04_example_data/example_metadata.csv"
OUTPUT_DIR <- "../results"

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Load data
count_data <- read.csv(COUNT_MATRIX, row.names = 1, check.names = FALSE)
sample_info <- read.csv(METADATA, row.names = 1)

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(
    countData = count_data,
    colData = sample_info,
    design = ~ Genotype + Treatment + Genotype:Treatment
)

# Filter low-count genes
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]

# Run DESeq2
dds <- DESeq(dds)

# Extract results for main comparison (IR64-Sub1 Combined vs Control)
res <- results(dds, contrast = list(c("Treatment_Combined_vs_Control")))

# Convert to data frame
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
res_df <- res_df[!is.na(res_df$padj), ]

# Add significance
res_df$significance <- "Not significant"
res_df$significance[res_df$padj < 0.05 & res_df$log2FoldChange > 1] <- "Upregulated"
res_df$significance[res_df$padj < 0.05 & res_df$log2FoldChange < -1] <- "Downregulated"

# Save results
write.csv(res_df, file.path(OUTPUT_DIR, "DE_results.csv"), row.names = FALSE)

# Print summary
cat("Upregulated:", sum(res_df$significance == "Upregulated"), "\n")
cat("Downregulated:", sum(res_df$significance == "Downregulated"), "\n")

# Save DESeq2 object
saveRDS(dds, file.path(OUTPUT_DIR, "deseq2_object.rds"))

cat("Analysis complete! Results saved to:", OUTPUT_DIR, "\n")
