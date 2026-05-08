#!/usr/bin/env Rscript
# =============================================================================
# Principal Component Analysis (PCA) for RNA-seq Samples
# =============================================================================
# Author: Waqas Yousaf
# Date: May 2026
# Usage: Rscript 03_pca_analysis.R
# =============================================================================

library(DESeq2)
library(ggplot2)
library(RColorBrewer)

OUTPUT_DIR <- "../results"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Load DESeq2 object
dds <- readRDS(file.path(OUTPUT_DIR, "deseq2_object.rds"))

# Variance stabilizing transformation
vsd <- vst(dds, blind = FALSE)

# Extract PCA data
pca_data <- plotPCA(vsd, intgroup = c("Genotype", "Treatment"), returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

# Create PCA plot
pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Genotype, shape = Treatment)) +
    geom_point(size = 4, alpha = 0.8) +
    stat_ellipse(aes(group = interaction(Genotype, Treatment)), linetype = 2) +
    labs(x = paste0("PC1: ", percentVar[1], "% variance"),
         y = paste0("PC2: ", percentVar[2], "% variance"),
         title = "PCA of RNA-seq Samples") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Save plot
ggsave(file.path(OUTPUT_DIR, "PCA_plot.png"), pca_plot, width = 10, height = 7, dpi = 300)

cat("PCA plot saved to:", file.path(OUTPUT_DIR, "PCA_plot.png"), "\n")
