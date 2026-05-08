#!/usr/bin/env Rscript
# =============================================================================
# Heatmap Analysis of Differentially Expressed Genes
# =============================================================================
# Author: Waqas Yousaf
# Date: May 2026
# Usage: Rscript 04_heatmap_analysis.R
# =============================================================================

library(DESeq2)
library(pheatmap)
library(RColorBrewer)

OUTPUT_DIR <- "../results"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Load DESeq2 object
dds <- readRDS(file.path(OUTPUT_DIR, "deseq2_object.rds"))
vsd <- vst(dds, blind = FALSE)

# Key genes of interest
key_genes <- c("OsIRT1", "OsYSL15", "OsIDS2", "OsNAS1", 
               "SUB1A-1", "OsADH1", "OsPDC1")

# Keep only genes present in dataset
key_genes_present <- key_genes[key_genes %in% rownames(assay(vsd))]

if (length(key_genes_present) > 0) {
    
    # Extract expression data
    expr_matrix <- assay(vsd)[key_genes_present, ]
    expr_matrix_scaled <- t(scale(t(expr_matrix)))
    
    # Sample annotations
    annotation_col <- data.frame(
        Genotype = dds$Genotype,
        Treatment = dds$Treatment
    )
    rownames(annotation_col) <- colnames(expr_matrix)
    
    # Generate heatmap
    pheatmap(expr_matrix_scaled,
        main = "Expression of Key Stress-Responsive Genes",
        clustering_method = "complete",
        show_rownames = TRUE,
        show_colnames = TRUE,
        annotation_col = annotation_col,
        color = colorRampPalette(c("blue", "white", "red"))(100),
        filename = file.path(OUTPUT_DIR, "Heatmap_key_genes.png"),
        width = 10, height = 8
    )
    
    cat("Heatmap saved to:", file.path(OUTPUT_DIR, "Heatmap_key_genes.png"), "\n")
}
