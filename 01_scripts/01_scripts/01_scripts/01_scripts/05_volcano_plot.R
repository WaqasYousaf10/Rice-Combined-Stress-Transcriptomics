#!/usr/bin/env Rscript
# =============================================================================
# Volcano Plot for Differential Expression
# =============================================================================
# Author: Waqas Yousaf
# Date: May 2026
# Usage: Rscript 05_volcano_plot.R
# =============================================================================

library(ggplot2)
library(ggrepel)

OUTPUT_DIR <- "../results"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Load DESeq2 results
res_df <- read.csv(file.path(OUTPUT_DIR, "DE_results.csv"), row.names = 1)

# Key genes to label
key_genes <- c("OsIRT1", "OsYSL15", "SUB1A-1", "OsADH1")
res_df$label <- ifelse(rownames(res_df) %in% key_genes, rownames(res_df), "")

# Create volcano plot
volcano_plot <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("Downregulated" = "blue", "Not significant" = "grey", "Upregulated" = "red")) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_text_repel(aes(label = label), size = 4, max.overlaps = 15) +
    labs(x = expression(log[2] ~ "Fold Change"),
         y = expression(-log[10] ~ "Adjusted p-value"),
         title = "Volcano Plot: Combined Stress vs Control") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Save plot
ggsave(file.path(OUTPUT_DIR, "Volcano_plot.png"), volcano_plot, width = 10, height = 8, dpi = 300)

cat("Volcano plot saved to:", file.path(OUTPUT_DIR, "Volcano_plot.png"), "\n")
