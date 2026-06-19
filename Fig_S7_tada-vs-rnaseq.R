##############################
## Generates Figure S7 from Marshall (2026) "damidBind: an R/Bioconductor package
## for differential DamID analysis and data exploration"
##

pacman::p_load(damidBind, tximport, DESeq2, ggplot2, dplyr, AnnotationHub, ggrepel, rhdf5, apeglm, clusterProfiler, org.Dm.eg.db, patchwork)

# Set seed for reproducibility
set.seed(42)

fdr_cut=0.05

#########################
## 1. Data Preparation
##

### 1.1 damidBind analysis
# using the RNA Pol II TaDa dataset (larval NSCs vs adult neurons).

nvsn <- load_data_genes(binding_profiles_path = c("data/RNAPII_TaDa/RNAPII*"),norm_method = "quantile")

nvsn.diff <- differential_binding(nvsn, cond=c("NSCs"="Wor","Neurons"="elav"))

# As a bonus, this produces one of the volcano plot panels in Fig. 1B:
plot_volcano(nvsn.diff,label_display = list(r=0.45),label_config = list(max_overlaps=5,clean_names=TRUE), fdr_filter_threshold = 0.01,
             save = list(
                 filename="Fig_1B_RNAPII",
                 format="svg",
                 width=6,
                 height=4
                ))

# Extract the analysis table from the differential analysis
tada_df <- analysisTable(nvsn.diff) %>%
  tibble::rownames_to_column("locus") %>%
  dplyr::select(gene_id, gene_name, tada_logFC = logFC, tada_fdr = adj.P.Val)

### 1.2 Prepare RNA-seq ID Mapping
# `kallisto` outputs transcript IDs (`FBtr`). We map these to gene IDs (`FBgn`)
# to match the `damidBind` output, using the same genome version from `AnnotationHub`.

# Fetch mapping from AnnotationHub (matching the version used in damidBind, currently 113)
ah <- AnnotationHub()
query_res <- query(ah, c("EnsDb", "Drosophila melanogaster", "113"))
edb <- query_res[[1]]

# Create tx2gene table
tx <- transcripts(edb, columns = c("tx_id", "gene_id"))
tx2gene <- as.data.frame(tx) %>%
  dplyr::select(tx_id, gene_id)


#######################################
## 2. RNA-seq differential expression via DESeq2
##

### 2.1 Import kallisto files
nsc_files <- c(
    "data/RNA-seq/GSE38764/SRR513583-kallisto/abundance.tsv",
    "data/RNA-seq/GSE38764/SRR513584-kallisto/abundance.tsv",
    "data/RNA-seq/GSE38764/SRR513585-kallisto/abundance.tsv",
    "data/RNA-seq/GSE38764/SRR513586-kallisto/abundance.tsv",
    "data/RNA-seq/GSE38764/SRR513587-kallisto/abundance.tsv"
)

neuron_files <- c(
    "data/RNA-seq/GSE235989/SRR25040906-kallisto/abundance.tsv",
    "data/RNA-seq/GSE235989/SRR25040907-kallisto/abundance.tsv",
    "data/RNA-seq/GSE235989/SRR25040908-kallisto/abundance.tsv"
)

files <- c(nsc_files, neuron_files)
names(files) <- c("NSC_r1", "NSC_r2", "NSC_r3", "NSC_r4", "NSC_r5", "Neuron_r1", "Neuron_r2", "Neuron_r3")

txi <- tximport(files, type = "kallisto", tx2gene = tx2gene)

# Setup sample metadata
sample_info <- data.frame(
  condition = factor(c(rep("NSC",5),rep("Neuron",3)), levels = c("NSC", "Neuron")),
  row.names = names(files)
)


### 2.2 Run DESeq2

txi <- tximport(
  files,
  type = "kallisto",
  tx2gene = tx2gene,
  ignoreTxVersion = FALSE
)

sample_info <- data.frame(
  condition = factor(
    c(rep("NSC", 5), rep("Neuron", 3)),
    levels = c("NSC", "Neuron")
  ),
  row.names = names(files)
)

dds <- DESeqDataSetFromTximport(
  txi,
  colData = sample_info,
  design = ~ condition
)

# low-count filter
keep <- rowSums(counts(dds) >= 10) >= 2
dds <- dds[keep, ]

dds <- DESeq(dds)

resultsNames(dds)

res_neuron_vs_nsc <- results(
  dds,
  name = "condition_Neuron_vs_NSC",
  alpha = fdr_cut
)

res_shrunk <- lfcShrink(
  dds,
  coef = "condition_Neuron_vs_NSC",
  type = "apeglm"
)

# Convert to NSC-positive for comparison with TaDa
rnaseq_df <- as.data.frame(res_shrunk) %>%
  tibble::rownames_to_column("gene_id") %>%
  transmute(
    gene_id,
    rna_logFC = -log2FoldChange,
    rna_fdr = padj,
    baseMean = baseMean
  )

############################
## 3. Concordance analysis
##


### 3.1 Merge datasets
comparison_df <- tada_df %>%
  inner_join(rnaseq_df, by = "gene_id") %>%
  dplyr::filter(!is.na(rna_logFC))

# Identify significant genes by category (for highlighting)
comparison_df <- comparison_df %>%
  mutate(sig_label = case_when(
    tada_fdr < fdr_cut & rna_fdr < fdr_cut ~ "Both Significant",
    tada_fdr < fdr_cut ~ "TaDa Only",
    rna_fdr < fdr_cut ~ "RNA-seq Only",
    TRUE ~ "Non-significant"
  ))

### 3.2 Correlation plots

comparison_sig_only <- comparison_df %>%
  dplyr::filter(sig_label != "Non-significant")

spearman_rho <- cor(comparison_sig_only$tada_logFC, comparison_sig_only$rna_logFC, method = "spearman")

# The sample_labels_by_isolation() method still needs some refining, but it works pretty well!
kept_mask <- damidBind::sample_labels_by_isolation(
    df = comparison_sig_only,
    x_col = "rna_logFC",
    y_col = "tada_logFC",
    r = 0.5,
    k_priority = 30,
    scale = TRUE
)

# Ensure dpn remains after thinning, if present.  The gene dpn is *the* classic NSC marker
# so if it was a hit we don't want to accidentally thin it out with the step above
# (note -- this doesn't stop max.overlaps in ggrepel still nuking it, which actually
# happens in the main plot with these settings.  Ah well ...)
kept_mask <- kept_mask | comparison_sig_only$gene_name == "dpn"

label_df <- comparison_sig_only[kept_mask, ]

pol_overview = ggplot(comparison_sig_only, aes(x = rna_logFC, y = tada_logFC)) +
  geom_hline(yintercept = 0, linetype = "dotted", alpha = 0.4) +
  geom_vline(xintercept = 0, linetype = "dotted", alpha = 0.4) +

  # Main point layer
  geom_point(aes(colour = sig_label), alpha = 0.3, size = 0.5) +

  # Regression line per facet
  geom_smooth(method = "lm", color = "black", linetype = "dashed", linewidth = 0.5) +

  # labels
  geom_text_repel(
    data = label_df,
    aes(label = gene_name),
    size = 3,
    max.overlaps = 8,
    box.padding = 0.2,
    point.padding = 0.2,
    min.segment.length = 0,
    segment.color = "grey50",
    segment.alpha = 0.5
  ) +

  # Aesthetics
  scale_colour_manual(values = c("Both Significant" = "firebrick",
                                 "TaDa Only" = "orange",
                                 "RNA-seq Only" = "steelblue"),
                      guide = guide_legend(title = "",override.aes = list(size = 4, alpha = 1))
                      ) +
  theme_bw(base_size = 14) +
  labs(
    title = "RNA Pol II TaDa vs RNA-seq",
    subtitle = sprintf("All significant loci: Spearman's Rho = %0.2f", spearman_rho),
    x = "RNA-seq Log2 Fold Change",
    y = "RNA Pol II TaDa Log2 Fold Change"
  ) +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey95"),
    panel.grid.minor = element_blank()
  )


facet_stats <- comparison_sig_only %>%
  group_by(sig_label) %>%
  summarise(
    rho = cor(rna_logFC, tada_logFC, method = "spearman", use = "pairwise.complete.obs"),
    n = n(),
    .groups = 'drop'
  ) %>%
  mutate(label = sprintf("Spearman's Rho: %0.3f", rho))

pol_facets = ggplot(comparison_sig_only, aes(x = rna_logFC, y = tada_logFC)) +
  geom_hline(yintercept = 0, linetype = "dotted", alpha = 0.4) +
  geom_vline(xintercept = 0, linetype = "dotted", alpha = 0.4) +

  geom_point(aes(colour = sig_label), alpha = 0.5, size = 0.5) +

  # Regression line per facet
  geom_smooth(method = "lm", color = "black", linetype = "dashed", linewidth = 0.5) +

  # Statistics labels per facet
  geom_text(
    data = facet_stats,
    aes(x = -18, y = 8, label = label),
    hjust = 0, vjust = 0.5,
    size = 5, inherit.aes = FALSE
  ) +

  geom_text_repel(
    data = label_df,
    aes(label = gene_name),
    size = 3,
    max.overlaps = 8,
    box.padding = 0.2,
    point.padding = 0.2,
    min.segment.length = 0,
    segment.color = "grey50",
    segment.alpha = 0.5
  ) +

  facet_wrap(~sig_label) +

  scale_colour_manual(values = c("Both Significant" = "firebrick",
                                 "TaDa Only" = "orange",
                                 "RNA-seq Only" = "steelblue")) +
  theme_bw(base_size = 14) +
  labs(
    title = "Significant differentially expressed loci",
    subtitle = "Faceted by significance status",
    x = "RNA-seq Log2 Fold Change",
    y = "RNA Pol II TaDa Log2 Fold Change"
  ) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "grey95"),
    panel.grid.minor = element_blank()
  )


categories <- c("Both Significant", "TaDa Only", "RNA-seq Only")
universe_ids <- unique(comparison_df$gene_id)

go_results_list <- lapply(categories, function(cat_name) {

  # Extract IDs for this category
  category_ids <- comparison_df %>%
    dplyr::filter(sig_label == cat_name) %>%
    pull(gene_id) %>%
    unique()

  message(sprintf("Running GO for %s (%d genes)...", cat_name, length(category_ids)))

  # Using enrichGO with the same parameters as damidBind ...
  ego <- enrichGO(
    gene          = category_ids,
    universe      = universe_ids,
    OrgDb         = org.Dm.eg.db,
    keyType       = "ENSEMBL",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2,
    readable      = TRUE
  )

  return(ego)
})
names(go_results_list) <- categories

plot_go_category <- function(ego, title) {
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
    return(ggplot() + labs(title = paste(title, " - No terms found")) + theme_void())
  }

  # Using damidBind's plotting internals
  df <- as.data.frame(ego) %>%
    arrange(p.adjust) %>%
    slice_head(n = 15) %>%
    mutate(GeneRatio = sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1])/as.numeric(x[2])))

  damidBind:::._create_go_dotplot(df, show_category = 15, fit_labels = F, wrap_labels = F, plot_title = title, label_format_width = 30, theme_size = 14)

}

# Generate and combine the plots
p_both <- plot_go_category(go_results_list[["Both Significant"]], "Both Significant")
#p_tada <- plot_go_category(go_results_list[["TaDa Only"]], "GO: TaDa Only") # <- no enrichment
p_rna  <- plot_go_category(go_results_list[["RNA-seq Only"]], "RNA-seq Only")

combined_go = (p_both + p_rna)

top_row <- (pol_overview | combined_go) +
  plot_layout(widths = c(0.578, 0.948))

# Final plot
final_plot <- (top_row / pol_facets) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 26))

pdf("Fig_S7.pdf",width = 15, height=12)
print(final_plot)
dev.off()


### 3.3 Directional concordance stats
sig_in_any <- comparison_df %>% dplyr::filter(tada_fdr < fdr_cut | rna_fdr < fdr_cut)
concordant <- sum(sign(sig_in_any$tada_logFC) == sign(sig_in_any$rna_logFC))
concordance_pct <- (concordant / nrow(sig_in_any)) * 100

message(sprintf("Directional concordance among significant genes (in any dataset): %.1f%%, (%d / %d)", concordance_pct, concordant , nrow(sig_in_any)))

sig_in_tada <- comparison_df %>% dplyr::filter(tada_fdr < fdr_cut)
concordant <- sum(sign(sig_in_tada$tada_logFC) == sign(sig_in_tada$rna_logFC))
concordance_pct <- (concordant / nrow(sig_in_tada)) * 100

message(sprintf("Directional concordance among significant TaDa genes: %.1f%% (%d / %d)", concordance_pct, concordant , nrow(sig_in_tada)))

sig_in_rnaseq <- comparison_df %>% dplyr::filter(rna_fdr < fdr_cut)
concordant <- sum(sign(sig_in_rnaseq$tada_logFC) == sign(sig_in_rnaseq$rna_logFC))
concordance_pct <- (concordant / nrow(sig_in_rnaseq)) * 100

message(sprintf("Directional concordance among significant RNA-seq genes: %.1f%% (%d / %d)", concordance_pct, concordant , nrow(sig_in_rnaseq)))

sig_in_tada <- comparison_df %>% dplyr::filter(tada_fdr < fdr_cut & abs(tada_logFC)>2)
concordant <- sum(sign(sig_in_tada$tada_logFC) == sign(sig_in_tada$rna_logFC))
concordance_pct <- (concordant / nrow(sig_in_tada)) * 100

message(sprintf("Directional concordance among significant TaDa genes (log2FC > 2): %.1f%% (%d / %d)", concordance_pct, concordant, nrow(sig_in_tada)))

sig_in_rnaseq <- comparison_df %>% dplyr::filter(rna_fdr < fdr_cut & abs(rna_logFC)>2)
concordant <- sum(sign(sig_in_rnaseq$tada_logFC) == sign(sig_in_rnaseq$rna_logFC))
concordance_pct <- (concordant / nrow(sig_in_rnaseq)) * 100

message(sprintf("Directional concordance among significant RNA-seq genes (log2FC > 2): %.1f%% (%d / %d)", concordance_pct, concordant, nrow(sig_in_rnaseq)))
