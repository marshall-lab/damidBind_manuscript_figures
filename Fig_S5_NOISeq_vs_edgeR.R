
# If you don't have pacman installed, change the below to standard R library loading calls
pacman::p_load(damidBind, edgeR, ggplot2, ggrepel, patchwork, dplyr)

# Custom function for the proportional Venn plot (see the function code for details)
source("plot_venn_gg.R")

#' Generates Figure S5 from Marshall (2026) "damidBind: an R/Bioconductor package for differential DamID analysis and data exploration"
#' This figure compares NOISeq (via damidBind) with edgeR on CATaDa data
#' For more details, see the manuscript
#'
#' @param loaded_data The output from `load_data_peaks()`
#' @param cond Named condition vector, e.g., c("L4" = "L4", "L5" = "L5")
#' @export
generate_catada_edgeR_figure <- function(loaded_data, cond, venn_repel=FALSE, venn_bp=0.1, venn_f=0.5,venn_fp=0.5) {

  message("Preparing data and running biological analyses...")

  # Prep data using the package internal helper
  pr_data <- damidBind:::prep_data_for_differential_analysis(loaded_data, cond)
  mat_continuous <- as.matrix(pr_data$mat)
  group <- pr_data$factors$condition
  occ_df <- pr_data$occupancy_df

  if (ncol(mat_continuous) != 6) {
    stop("The figure logic requires 3 replicates per condition (6 samples total).")
  }

  # NOISeq (true biological condition) via damidBind
  res_ns_true <- differential_accessibility(loaded_data, cond = cond)
  ns_sig_loci <- unique(c(rownames(enrichedCond1(res_ns_true)),
                          rownames(enrichedCond2(res_ns_true))))
  n_ns_bio <- length(ns_sig_loci)

  # edgeR (true biological condition)
  y <- DGEList(counts = mat_continuous, group = group)
  y <- calcNormFactors(y)
  y <- estimateDisp(y)
  design <- model.matrix(~0 + group)
  colnames(design) <- levels(group)

  fit <- glmQLFit(y, design)
  contrast_str <- paste0(levels(group)[1], " - ", levels(group)[2])
  con <- makeContrasts(contrasts = contrast_str, levels = design)

  res_er_true <- topTags(glmQLFTest(fit, contrast = con), n = Inf)$table
  res_er_true$id <- rownames(res_er_true)
  res_er_true$gene_name <- occ_df$gene_name[match(res_er_true$id, rownames(occ_df))]

  er_sig_loci <- res_er_true$id[res_er_true$FDR < 0.05]
  n_er_bio <- length(er_sig_loci)

  # Hits for highlighting later
  confirmed_hits <- intersect(ns_sig_loci, er_sig_loci)
  ns_only_hits <- setdiff(ns_sig_loci, er_sig_loci)
  er_only_hits <- setdiff(er_sig_loci, ns_sig_loci)

  # Permuted null conditions
  message("Running 9 Combinatorial Null Permutations...")

  # Generate all unique splits of 6 into 2 groups of 3.
  all_combs <- combn(6, 3)
  unique_splits <- all_combs[, 1:10]

  # Identify the biological split (indices 1, 2, 3)
  bio_idx <- which(apply(unique_splits, 2, function(x) all(x <= 3)))
  null_splits <- unique_splits[, -bio_idx]

  ns_null_hits <- numeric(9)
  er_null_hits <- numeric(9)

  temp_data <- loaded_data

  # ... and we loop through all permuted null conditions:
  for (i in seq_len(ncol(null_splits))) {
    idx_A <- null_splits[, i]
    idx_B <- setdiff(1:6, idx_A)

    samples_A <- colnames(mat_continuous)[idx_A]
    samples_B <- colnames(mat_continuous)[idx_B]

    # Regex strings to force groups A and B based on the current partition for damidBind
    p_A <- paste0("^(", paste(samples_A, collapse="|"), ")$")
    p_B <- paste0("^(", paste(samples_B, collapse="|"), ")$")
    cond_perm <- c("NullA" = p_A, "NullB" = p_B)

    # Run NOISeq (via damidBind) on this null partition
    res_perm_ns <- differential_accessibility(temp_data, cond = cond_perm, regex = TRUE)
    ns_null_hits[i] <- nrow(enrichedCond1(res_perm_ns)) + nrow(enrichedCond2(res_perm_ns))

    # Run edgeR on this null partition
    perm_group <- factor(ifelse(colnames(mat_continuous) %in% samples_A, "NullA", "NullB"), levels = c("NullA", "NullB"))
    y_p <- DGEList(counts = mat_continuous, group = perm_group)
    y_p <- calcNormFactors(y_p)
    y_p <- estimateDisp(y_p)
    fit_p <- glmQLFit(y_p, model.matrix(~perm_group))
    res_perm_er <- topTags(glmQLFTest(fit_p, coef = 2), n = Inf)$table
    er_null_hits[i] <- sum(res_perm_er$FDR < 0.05)

    message(sprintf(" - Partition %d: NOISeq=%d, edgeR=%d", i, ns_null_hits[i], er_null_hits[i]))
  }

  # empirical false-positive fraction (mean hits in null / real hits)
  efpf_ns <- mean(ns_null_hits) / n_ns_bio
  efpf_er <- mean(er_null_hits) / n_er_bio

  message("\n## Empirical false-positive fraction (eFPF) summary:")
  message(sprintf("NOISeq: Bio=%d, Null Mean=%.1f (Max=%d), eFPF=%.1f%%",
                  n_ns_bio, mean(ns_null_hits), max(ns_null_hits), efpf_ns * 100))
  message(sprintf("edgeR:  Bio=%d, Null Mean=%.1f (Max=%d), eFPF=%.1f%%",
                  n_er_bio, mean(er_null_hits), max(er_null_hits), efpf_er * 100))

  df_bio <- data.frame(
    Test   = "Biological",
    Method = c("NOISeq (Non-parametric)", "edgeR (Negative Binomial)"),
    Hits   = c(n_ns_bio, n_er_bio),
    ymin   = c(n_ns_bio, n_er_bio),
    ymax   = c(n_ns_bio, n_er_bio)
  )

  df_null <- data.frame(
    Test   = "Permuted null (n=9)",
    Method = c("NOISeq (Non-parametric)", "edgeR (Negative Binomial)"),
    Hits   = c(mean(ns_null_hits), mean(er_null_hits)),
    ymin   = c(min(ns_null_hits), min(er_null_hits)),
    ymax   = c(max(ns_null_hits), max(er_null_hits))
  )

  df_bar <- rbind(df_bio, df_null)
  df_bar$Test <- factor(df_bar$Test, levels = c("Biological", "Permuted null (n=9)"))
  df_bar$show_errorbar <- df_bar$Test == "Permuted null (n=9)"


  sig_bg_col <- "#FFA500"
  highlight_col <- "darkred"

  # Panel A (bar chart of real hits / permuted null hits)
  pA <- ggplot(df_bar, aes(x = Test, y = Hits, fill = Method)) +
      geom_bar(
        stat = "identity",
        position = position_dodge(width = 0.8),
        width = 0.7
      ) +
      geom_errorbar(
        data = subset(df_bar, show_errorbar),
        aes(ymin = ymin, ymax = ymax),
        position = position_dodge(width = 0.8),
        width = 0.25
      ) +
      scale_fill_manual(values = c(
        "NOISeq (Non-parametric)" = "#56B4E9",
        "edgeR (Negative Binomial)" = "#E69F00"
      )) +
      geom_text(
        aes(label = round(Hits, 0)),
        position = position_dodge(width = 0.8),
        vjust = -1,
        size = 5
      ) +
      theme_bw(base_size = 18) +
      labs(
        title = "Significant CATaDa loci",
        y = "Number of sig loci",
        x = NULL
      ) +
      ylim(0, max(df_bar$ymax, na.rm = TRUE) * 1.1) +
      theme(
        legend.position = "bottom",
        legend.title = element_blank()
      )

  # Panel B: Venn diagram
  pB <- draw_venn_gg(list_x = er_sig_loci,
                     list_y = ns_sig_loci,
                     xtitle = "edgeR",
                     ytitle = "NOISeq",
                     palette = c("#E69F00", "#56B4E9", "#666666"),
                     nrtype = "abs",
                     set_label_position = "outside_unique",
                     set_count_nrtype = "none",
                     repel_labels = venn_repel,
                     repel_box_padding = venn_bp,
                      repel_point_padding = 0,
                      repel_force = venn_f,
                      repel_force_pull = venn_fp,
                     stroke_width = 1,
                     title = NULL,
                     label_size=6,
                     set_label_size=6)

  # Panel C: NOISeq (shared hits highlighted)
  pC <- plot_volcano(res_ns_true,
                     highlight = list("Confirmed" = confirmed_hits),
                     labels = "highlight", label_display = list(r = 0.5),
                     plot_config = list(title = "NOISeq (union with edgeR)",
                                        sig_colour = sig_bg_col, sig_alpha = 0.2, nonsig_alpha = 0.1),
                     highlight_config = list(colour = highlight_col, size = 0.8, alpha = 0.4)) + theme(legend.position = "none")

  # Panel D: edgeR (shared hits highlighted)
  df_er_d <- res_er_true %>%
      mutate(minuslogp = -log10(FDR),
             status = case_when(id %in% ns_sig_loci ~ "Highlight",
                                FDR < 0.05 ~ "Significant",
                                TRUE ~ "NS"))

  label_d_targets <- df_er_d %>% dplyr::filter(status == "Highlight")
  kept_d <- if(nrow(label_d_targets) > 0) sample_labels_by_isolation(label_d_targets, "logFC", "minuslogp", r = 0.6) else logical(0)

  pD <- ggplot(df_er_d, aes(x = logFC, y = minuslogp, colour = status)) +
    geom_point(alpha = 0.2, size = 0.8, shape = 20) +
    scale_color_manual(values = c("Highlight" = highlight_col, "Significant" = sig_bg_col, "NS" = "grey35")) +
    geom_text_repel(data = label_d_targets[kept_d, ], aes(label = gene_name), size = 3, max.overlaps = 8,
                    min.segment.length = 0, segment.alpha = 0.5, color = "black") +
    labs(title = "edgeR (union with NOISeq)", x = bquote(log[2]*"FC (L4 neurons / L5 neurons)"), y = bquote("-log"[10] ~ "FDR")) +
    theme_bw(base_size = 18) + theme(legend.position = "none")

  # Panel E: NOISeq exclusive
  pE <- plot_volcano(res_ns_true,
                     highlight = list("NOISeq Exclusive" = ns_only_hits),
                     labels = "NOISeq Exclusive", label_display = list(r = 0.5),
                     plot_config = list(title = "NOISeq (exclusive hits)",
                                        sig_colour = sig_bg_col, sig_alpha = 0.2, nonsig_alpha = 0.1),
                     highlight_config = list(colour = highlight_col, size = 1, alpha = 0.6)) + theme(legend.position = "none")

  # Panel F: edgeR exclusive
  df_er_f <- res_er_true %>%
      mutate(minuslogp = -log10(FDR),
             status = case_when(id %in% er_only_hits ~ "Highlight",
                                FDR < 0.05 ~ "Significant",
                                TRUE ~ "NS"))

  label_f_targets <- df_er_f %>% dplyr::filter(status == "Highlight")
  kept_f <- if(nrow(label_f_targets) > 0) sample_labels_by_isolation(label_f_targets, "logFC", "minuslogp", r = 0.5) else logical(0)

  pF <- ggplot(df_er_f, aes(x = logFC, y = minuslogp, colour = status)) +
    geom_point(alpha = 0.2, size = 0.8, shape = 20) +
    scale_color_manual(values = c("Highlight" = highlight_col, "Significant" = sig_bg_col, "NS" = "grey35")) +
    geom_text_repel(data = label_f_targets[kept_f, ], aes(label = gene_name), size = 3, max.overlaps = 8,
                    min.segment.length = 0, segment.alpha = 0.5, color = "black") +
    labs(title = "edgeR (exclusive hits)", x = bquote(log[2]*"FC (L4 neurons / L5 neurons)"), y = bquote("-log"[10] ~ "FDR")) +
    theme_bw(base_size = 18) + theme(legend.position = "none")

  pB2 <- pB +
      theme(
        plot.margin = margin(t = 50, r = 50, b = 50, l = 50)
      )

  # Putting all of that finally together ...
  final_plot <- ((pA | pB2) / (pC | pD) / (pE | pF)) +
    plot_annotation(tag_levels = 'A') &
    theme(plot.tag = element_text(face = "bold", size = 20))

  print(final_plot)

  if (length(confirmed_hits) > 2) {
      ns_tab <- analysisTable(res_ns_true)[confirmed_hits, ]
      er_tab <- res_er_true[confirmed_hits, ]
      cor_lfc <- cor(ns_tab$logFC, er_tab$logFC, method = "spearman")
      message(sprintf("\nSpearman's Correlation (n=%d confirmed hits) Log2FC: %0.3f",
                      length(confirmed_hits), cor_lfc))
  }

  return(invisible(list(ns_res = res_ns_true, er_res = res_er_true,
                        ns_null = ns_null_hits, er_null = er_null_hits,
                        stats = list(efpf_ns = efpf_ns, efpf_er = efpf_er))))
}

# Load up the CATaDa data via damidBind
xucat_qnorm = load_data_peaks(
    binding_profiles_path = "data/catada/*ph.gz",
    peaks_path = "data/catada/*gff.gz",
    pre_scale=FALSE,
    norm_method = "quantile"
    )

# Analyse and plot ...
png("Fig_S5.png",width=12*200,height=16*200,res=200)
catada_comparison = generate_catada_edgeR_figure(
    xucat_qnorm, c("L4 neurons"="L4","L5 neurons"="L5"),venn_repel = T,venn_f = 0.1                  )
dev.off()


