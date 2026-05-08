#!/usr/bin/env Rscript
# =============================================================================
# Gene Ontology (GO) Enrichment Analysis
# =============================================================================
# Author: Waqas Yousaf
# Date: May 2026
# Usage: Rscript 06_go_enrichment.R
# =============================================================================

library(clusterProfiler)
library(org.Osativa.eg.db)
library(ggplot2)

OUTPUT_DIR <- "../results"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Load DESeq2 results
res_df <- read.csv(file.path(OUTPUT_DIR, "DE_results.csv"), row.names = 1)

# Extract upregulated genes (log2FC > 1, padj < 0.05)
up_genes <- rownames(res_df[res_df$log2FoldChange > 1 & res_df$padj < 0.05, ])

if (length(up_genes) > 0) {
    
    # Convert to Entrez IDs
    up_entrez <- bitr(up_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Osativa.eg.db)
    
    # GO enrichment
    go_results <- enrichGO(gene = up_entrez$ENTREZID,
                           OrgDb = org.Osativa.eg.db,
                           ont = "BP",
                           pAdjustMethod = "BH",
                           pvalueCutoff = 0.05,
                           readable = TRUE)
    
    if (nrow(go_results) > 0) {
        # Save results
        write.csv(as.data.frame(go_results), file.path(OUTPUT_DIR, "GO_enrichment.csv"), row.names = FALSE)
        
        # Create dotplot
        go_plot <- dotplot(go_results, showCategory = 15, title = "GO Enrichment: Upregulated Genes")
        ggsave(file.path(OUTPUT_DIR, "GO_enrichment.png"), go_plot, width = 10, height = 8, dpi = 300)
        
        cat("GO enrichment results saved to:", OUTPUT_DIR, "\n")
    }
}
