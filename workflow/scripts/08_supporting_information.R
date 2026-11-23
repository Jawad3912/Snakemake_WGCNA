source(file.path(snakemake@scriptdir, "functions.R"))
source(file.path(snakemake@scriptdir, "plot_functions.R"))
library(ggrepel)
library(openxlsx)

data <- readRDS(snakemake@input[["data"]])
meta <- readRDS(snakemake@input[["meta"]])
powers <- readRDS(snakemake@input[["powers"]])
network <- readRDS(snakemake@input[["network"]])
module_trait <- readRDS(snakemake@input[["module_trait"]])
selection <- readRDS(snakemake@input[["selection"]])
hub_table <- readRDS(snakemake@input[["hubs"]])$table
enrichment <- readRDS(snakemake@input[["enrichment"]])
mt <- module_trait$mt
ms <- module_trait$ms
colors <- network$colors
selected <- selection$selected
pos <- selected[["positive"]]
neg <- selected[["negative"]]
other <- setdiff(ms$module[ms$eligible & ms$pooled_r < 0], neg)[1]
if (is.na(other)) stop("S11 Fig needs a second eligible negative module, but none was found")
TH <- theme_figure(1)

figure_files <- output_paths("figures")
table_files <- output_paths("tables")

save_si <- function(p, name, w, h) {
  ggsave(figure_files[[name]], p, width = w, height = h, device = cairo_pdf, bg = "white")
}

tab <- function(name) read.csv(snakemake@input[[name]], check.names = FALSE)

xl <- function(name, sheets) write.xlsx(sheets, table_files[[name]], overwrite = TRUE)

selected_first <- function(modules) c(pos, neg, setdiff(modules, c(pos, neg)))

module_title <- function(m) paste(cap(m), "module")

direction_label <- function(m, direction) paste0(cap(m), " (", direction, ")")

pair_label <- function(x) paste(deparse(italic_species(x, TRUE)), collapse = "")

left_out_label <- function(x) ifelse(x == "none", "All datasets", paste0("Without ", SET_LABELS[x]))

italic_or_bold_pooled <- function(x) {
  as.expression(lapply(x, function(s) if (grepl("^Pooled", s)) bquote(bold(.(s))) else italic_species(s)))
}

two_line_natural <- function(x) {
  as.expression(lapply(x, function(s) if (s == "Natural mastitis") quote(atop("Natural", "mastitis")) else italic_species(s)))
}

two_line_without <- function(x) {
  as.expression(lapply(x, function(s) bquote(atop("Without", .(italic_species(sub("Without ", "", s)))))))
}

preservation_data <- function() {
  pd <- selection$preservation %>%
    filter(!module %in% c("grey", "gold")) %>%
    mutate(pair = paste0(SET_LABELS[reference], " vs ", SET_LABELS[test]))
  pd %>% mutate(pair = factor(sapply(pair, pair_label), sapply(sort(unique(pair)), pair_label)))
}

preservation_plot <- function(pd, y, ylab, reverse = FALSE, lines = FALSE) {
  p <- ggplot(pd, aes(module_size, .data[[y]]))
  if (lines) {
    p <- p + geom_hline(yintercept = cfg$min_zsummary, linetype = 2, colour = "red") + geom_hline(yintercept = 2, linetype = 3, colour = "blue")
  }
  p <- p +
    geom_point(aes(fill = module), shape = 21, size = 2.6, stroke = 0.4) +
    geom_text_repel(aes(label = module),
      size = 2.5, max.overlaps = Inf, min.segment.length = 0.2, segment.size = 0.2,
      segment.colour = "grey50", box.padding = 0.15, point.padding = 0.1, seed = 1
    ) +
    scale_fill_identity() +
    scale_x_log10() +
    facet_wrap(~pair, ncol = 3, labeller = label_parsed) +
    labs(x = "Module size (genes, log scale)", y = ylab) +
    TH
  if (reverse) p <- p + scale_y_reverse()
  p
}

plot_average_preservation <- function(pd) {
  average <- pd %>%
    group_by(module) %>%
    summarise(Zsummary = mean(Zsummary), medianRank = mean(medianRank), .groups = "drop") %>%
    mutate(selected = module %in% selected)
  ggplot(average, aes(Zsummary, medianRank)) +
    geom_vline(xintercept = cfg$min_zsummary, linetype = 2, colour = "red") +
    geom_point(aes(fill = module, size = selected, stroke = ifelse(selected, 1.4, 0.4)), shape = 21) +
    geom_text_repel(aes(label = module, fontface = ifelse(selected, "bold", "plain")),
      size = 3, max.overlaps = Inf,
      segment.size = 0.2, segment.colour = "grey50", box.padding = 0.25, seed = 1
    ) +
    scale_size_manual(values = c(`FALSE` = 3, `TRUE` = 4.5), guide = "none") +
    scale_fill_identity() +
    scale_y_reverse() +
    labs(x = expression("Average" ~ Z[summary]), y = "Average MedianRank") +
    TH
}

soft_threshold_panel <- function(pt, y, ylab, r2_line = FALSE) {
  chosen <- pt %>% filter(power == powers[set])
  p <- ggplot(pt, aes(power, .data[[y]], colour = dataset))
  if (r2_line) p <- p + geom_hline(yintercept = cfg$r2_cut, linetype = 2)
  p +
    geom_line() +
    geom_point(size = 0.9) +
    geom_point(data = chosen, size = 4, shape = 21, stroke = 1.3) +
    scale_colour_manual(values = dcol, name = NULL, labels = italic_labels) +
    labs(x = expression("Soft-thresholding power" ~ beta), y = ylab) +
    TH +
    theme(legend.text = element_text(size = 10.5))
}

plot_forest <- function() {
  eligible_modules <- unique(c(ms$module[ms$eligible], selected))
  pooled_label <- "Pooled (random effects)"
  forest <- mt %>%
    filter(module %in% eligible_modules) %>%
    mutate(dataset = as.character(SET_LABELS[dataset])) %>%
    bind_rows(ms %>%
      filter(module %in% eligible_modules) %>%
      transmute(module, dataset = pooled_label, r = pooled_r, ci_low = pooled_ci_low, ci_high = pooled_ci_high)) %>%
    mutate(
      pooled = dataset == pooled_label,
      dataset = factor(dataset, rev(c(SET_LABELS, pooled_label))),
      module = factor(module_title(module), module_title(selected_first(eligible_modules)))
    )
  p <- ggplot(forest, aes(r, dataset, xmin = ci_low, xmax = ci_high)) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_vline(xintercept = c(-cfg$r_cut, cfg$r_cut), linetype = 3, colour = "grey50") +
    geom_errorbar(aes(linewidth = pooled), width = 0.25) +
    geom_point(aes(shape = pooled, size = pooled)) +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 18), guide = "none") +
    scale_size_manual(values = c(`FALSE` = 2.5, `TRUE` = 4), guide = "none") +
    scale_linewidth_manual(values = c(`FALSE` = 0.5, `TRUE` = 0.9), guide = "none") +
    scale_x_continuous(limits = c(-1, 1)) +
    facet_wrap(~module, ncol = 1) +
    labs(x = "Module-trait correlation r (95% CI)", y = NULL) +
    TH +
    scale_y_discrete(labels = italic_or_bold_pooled)
  list(plot = p, height = 1.7 * n_distinct(forest$module) + 1)
}

plot_merge_height_modules <- function(merge_height) {
  ggplot(merge_height, aes(mergeCutHeight, n_modules)) +
    geom_vline(xintercept = cfg$merge_cut, linetype = 2, colour = "red") +
    geom_line() +
    geom_point(size = 2.2) +
    labs(x = "Merge cut height", y = "Number of modules") +
    TH
}

plot_merge_height_jaccard <- function(merge_height) {
  pos_label <- direction_label(pos, "positive")
  neg_label <- direction_label(neg, "negative")
  merge_height %>%
    select(mergeCutHeight, positive_jaccard, negative_jaccard) %>%
    pivot_longer(-mergeCutHeight, names_to = "module", values_to = "jaccard") %>%
    mutate(module = ifelse(module == "positive_jaccard", pos_label, neg_label)) %>%
    ggplot(aes(mergeCutHeight, jaccard, colour = module, shape = module)) +
    geom_vline(xintercept = cfg$merge_cut, linetype = 2, colour = "red") +
    geom_line() +
    geom_point(size = 2.2) +
    scale_colour_manual(values = setNames(c("black", "#1f4fd1"), c(pos_label, neg_label))) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(x = "Merge cut height", y = "Jaccard index", colour = NULL, shape = NULL) +
    TH +
    theme(legend.position = "inside", legend.position.inside = c(0.03, 0.03), legend.justification = c(0, 0), legend.background = element_rect(colour = "grey70"))
}

plot_hub_gene_significance <- function() {
  hub_gs <- hub_table %>%
    filter(hub, module %in% c(pos, neg)) %>%
    select(gene_id, module, starts_with("GS_"), -starts_with("GS_p")) %>%
    pivot_longer(starts_with("GS_"), names_to = "dataset", values_to = "GS") %>%
    mutate(
      dataset = factor(SET_LABELS[sub("GS_", "", dataset)], DS),
      panel = factor(hub_genes_title(module), hub_genes_title(c(pos, neg)))
    )
  ggplot(hub_gs, aes(dataset, GS)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_violin(fill = "grey88", colour = "grey30") +
    geom_point(position = position_jitter(width = 0.13, seed = 1), size = 0.6, alpha = 0.55) +
    facet_wrap(~panel) +
    scale_y_continuous(limits = c(-1, 1)) +
    scale_x_discrete(labels = italic_labels) +
    labs(x = NULL, y = "Gene significance (signed)") +
    TH +
    theme(axis.text.x = element_text(size = 10.5))
}

plot_expression_overlap <- function() {
  main_modules <- ms %>% filter(eligible)
  directions <- setNames(main_modules$direction, main_modules$module)
  module_label <- function(m) direction_label(m, directions[m])
  expr_modules <- tab("expression_network_module_assignment")
  expr_colors <- setNames(expr_modules$module, expr_modules$gene_id)
  overlap <- map_dfr(names(directions), function(m) {
    tb <- table(expr_colors[intersect(module_genes(colors, m), names(expr_colors))])
    data.frame(reference = module_label(m), new_module = names(tb), pct = 100 * as.integer(tb) / sum(colors == m))
  }) %>%
    group_by(new_module) %>%
    filter(max(pct) >= 5) %>%
    ungroup() %>%
    mutate(reference = factor(reference, rev(module_label(selected_first(names(directions))))))
  ggplot(overlap, aes(new_module, reference, fill = pct)) +
    geom_tile(colour = "white", linewidth = 0.8) +
    geom_text(aes(label = round(pct), colour = pct > 45), size = 3.6, show.legend = FALSE) +
    scale_colour_manual(values = c(`FALSE` = "black", `TRUE` = "white")) +
    scale_fill_gradient(low = "#fde0dd", high = "#b2182b", name = "% of genes") +
    labs(x = "Module in expression-filtered network", y = "Module in main network") +
    TH +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())
}

plot_expression_module_trait <- function() {
  expr_main <- tab("expression_network_main_modules")
  expr_ms <- tab("expression_network_pooled")
  shown <- unique(c(na.omit(expr_main$best_match), expr_ms$module[expr_ms$eligible]))
  expr_heat <- tab("expression_network_module_trait") %>%
    filter(module %in% shown) %>%
    mutate(
      dataset = factor(SET_LABELS[dataset], DS),
      module = factor(module, expr_ms$module[order(expr_ms$pooled_r)]),
      lab = paste0(sprintf("%.2f", r), ifelse(abs(r) >= cfg$r_cut & fdr < cfg$fdr, "*", ""))
    )
  ggplot(expr_heat, aes(dataset, module, fill = r)) +
    geom_tile(colour = "white", linewidth = 0.8) +
    geom_text(aes(label = lab), size = 3.6) +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0, limits = c(-1, 1), name = "r") +
    labs(x = NULL, y = "Module in expression-filtered network") +
    TH +
    scale_x_discrete(labels = two_line_natural) +
    theme(axis.text.x = element_text(size = 10), panel.grid = element_blank())
}

plot_leave_one_out_pooled <- function() {
  loo <- tab("leave_one_out_pooled_r") %>%
    mutate(
      set = factor(left_out_label(left_out), rev(left_out_label(c("none", SETS)))),
      all = left_out == "none",
      module = factor(module_title(module), module_title(selected_first(unique(module))))
    )
  ggplot(loo, aes(pooled_r, set, xmin = pooled_ci_low, xmax = pooled_ci_high)) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_vline(xintercept = c(-cfg$r_cut, cfg$r_cut), linetype = 3, colour = "grey50") +
    geom_errorbar(width = 0.25) +
    geom_point(aes(shape = all), size = 2.8) +
    scale_shape_manual(values = c(`TRUE` = 18, `FALSE` = 16), guide = "none") +
    scale_x_continuous(limits = c(-1, 1)) +
    scale_y_discrete(labels = italic_labels) +
    facet_wrap(~module, ncol = 1) +
    labs(x = "Pooled module-trait r (95% CI)", y = NULL) +
    TH
}

plot_leave_one_out_modules <- function() {
  loo_modules <- tab("leave_one_out_main_modules")
  row_order <- loo_modules %>%
    distinct(reference_module, reference_direction) %>%
    arrange(match(reference_module, c(pos, neg)), reference_module) %>%
    with(direction_label(reference_module, reference_direction))
  loo_modules <- loo_modules %>%
    mutate(
      left_out = factor(left_out_label(left_out), left_out_label(SETS)),
      module = factor(direction_label(reference_module, reference_direction), rev(row_order)),
      txt = paste0(
        cap(best_match), "\nJ = ", sprintf("%.2f", jaccard), "\nHubs: ", round(pct_hubs_still_hub), "%",
        ifelse(best_eligible, "", "\nnot eligible")
      )
    )
  ggplot(loo_modules, aes(left_out, module, fill = jaccard)) +
    geom_tile(colour = "white", linewidth = 0.8) +
    geom_text(aes(label = txt), size = 3.4, lineheight = 0.9) +
    scale_fill_gradient(low = "white", high = "steelblue", limits = c(0, 1), name = "Jaccard\nindex") +
    scale_x_discrete(labels = two_line_without) +
    labs(x = NULL, y = "Module in main network") +
    TH +
    theme(panel.grid = element_blank())
}

software_versions <- function(ensembl) {
  upstream <- cfg$upstream_software
  r_packages <- c("DESeq2", "WGCNA", "metafor", "gprofiler2")
  bind_rows(
    data.frame(step = upstream$step, software = names(upstream$versions), version = unlist(upstream$versions, use.names = FALSE)),
    data.frame(
      step = "Downstream analysis (R)",
      software = c("R", "DESeq2", "WGCNA", "metafor", paste0("gprofiler2 (g:Profiler, Ensembl ", ensembl, ")")),
      version = c(paste(R.version$major, R.version$minor, sep = "."), sapply(r_packages, function(p) as.character(packageVersion(p))))
    )
  )
}

selection_rationale <- function() {
  selection$candidates %>%
    group_by(direction) %>%
    filter(n() > 1) %>%
    summarise(text = paste0(
      tools::toTitleCase(direction[1]), " module: ", tools::toTitleCase(module[selected]), " selected over ",
      paste(tools::toTitleCase(module[!selected]), collapse = ", "),
      " (min Zsummary ", paste(sprintf("%.1f", min_Zsummary), collapse = " vs "),
      "; mean median rank ", paste(sprintf("%.1f", mean_medianRank), collapse = " vs "),
      "; I2 ", paste(sprintf("%.0f%%", I2), collapse = " vs "), ")"
    )[1]) %>%
    pull(text)
}

analysis_settings <- function(ensembl) {
  fmt <- function(x) format(x, big.mark = ",")
  n_sets <- length(SETS)
  tribble(
    ~setting, ~value,
    "Datasets", paste0(SET_LABELS, " (", dataset_cfg("bioproject"), ")", collapse = ", "),
    "Network gene set", paste0("Genes with weighted Stouffer meta-analysis FDR < ", cfg$fdr, " (", fmt(length(colors)), " genes)"),
    "Expression transformation", "Variance-stabilizing transformation (DESeq2) applied per dataset to the network genes",
    "Soft-thresholding powers", paste0(SET_LABELS[names(powers)], " = ", powers, collapse = "; "),
    "Network type", "Signed",
    "Minimum module size", as.character(cfg$min_module_size),
    "deepSplit", as.character(cfg$deep_split),
    "Module merging threshold (mergeCutHeight)", as.character(cfg$merge_cut),
    "Module-trait association criteria", paste0(
      "Same sign in all ", n_sets, " datasets; |r| >= ", cfg$r_cut, " and BH-FDR < ", cfg$fdr,
      " in all ", n_sets, " datasets; pooled random-effects FDR < ", cfg$fdr
    ),
    "Module selection", paste(c(
      paste0(
        "Strong preservation (Zsummary > ", cfg$min_zsummary, ") in all reference-test comparisons; ",
        "when several candidates remained in one direction, they were compared by minimum and median Zsummary, ",
        "mean median rank, and pooled random-effects correlation with its 95% CI and I2"
      ),
      selection_rationale()
    ), collapse = ". "),
    "Hub gene definition", paste0("Signed consensus MM > ", cfg$mm_cut, " and consensus GS > ", cfg$gs_cut, " in the direction of the module-trait association"),
    "Module preservation permutations", as.character(cfg$n_perm_preservation),
    "Module-trait permutations", paste(cfg$n_perm_trait, "random permutations of disease labels per dataset"),
    "Functional enrichment", paste0(
      "g:Profiler (Ensembl ", ensembl, "), GO:BP/MF/CC and KEGG, BH-FDR < ", cfg$fdr, "; background: the ", fmt(nrow(meta)),
      " genes tested in the meta-analysis (main analysis), with all annotated Bos taurus genes as a comparison background (S8 Fig)"
    )
  ) %>%
    as.data.frame()
}

enrichment_sheets <- function() {
  selected_sets <- c(paste0("positive_", pos, c("_hubs", "_module")), paste0("negative_", neg, c("_hubs", "_module")))
  other_modules <- setdiff(ms$module[ms$eligible], c(pos, neg))
  sheets <- list(Selected_modules = enrichment %>% filter(query_set %in% selected_sets))
  for (m in other_modules) {
    sheets[[paste0(cap(m), "_module")]] <- enrichment %>% filter(grepl(paste0("^", m, "_"), query_set))
  }
  sheets$Term_counts <- sheets$Selected_modules %>% count(query_set, background, source, name = "significant_terms")
  sheets
}

pd <- preservation_data()
save_si(preservation_plot(pd, "Zsummary", expression(Z[summary]), lines = TRUE), "S1_Fig", 11, 7.5)
save_si(preservation_plot(pd, "medianRank", "MedianRank", reverse = TRUE), "S2_Fig", 11, 7.5)
save_si(plot_average_preservation(pd), "S3_Fig", 7, 6)

pt <- tab("soft_threshold_all_powers") %>% mutate(set = dataset, dataset = factor(SET_LABELS[dataset], DS))
s4a <- soft_threshold_panel(pt, "signed_R2", expression("Scale-free fit index" ~ R^2 ~ "(signed)"), r2_line = TRUE) + labs(tag = "A")
s4b <- soft_threshold_panel(pt, "mean_k", "Mean connectivity (log scale)") + scale_y_log10() + labs(tag = "B")
save_si(s4a + s4b + plot_layout(guides = "collect") & theme(legend.position = "bottom"), "S4_Fig", 10, 4.8)

forest <- plot_forest()
save_si(forest$plot, "S5_Fig", 7, forest$height)

merge_height <- tab("merge_height_sensitivity")
s6a <- plot_merge_height_modules(merge_height) + labs(tag = "A")
s6b <- plot_merge_height_jaccard(merge_height) + labs(tag = "B")
save_si(s6a + s6b, "S6_Fig", 10, 4.3)

save_si(plot_hub_gene_significance(), "S7_Fig", 8, 4.6)

s8a <- bubble_plot(hub_enrichment(enrichment, paste0("positive_", pos), "genome"), hub_genes_title(pos), n = 10) + labs(tag = "A")
s8b <- bubble_plot(hub_enrichment(enrichment, paste0("negative_", neg), "genome"), hub_genes_title(neg), n = 10) + labs(tag = "B")
save_si(s8a | s8b, "S8_Fig", 12, 11.5)

s9a <- plot_expression_overlap() + labs(tag = "A")
s9b <- plot_expression_module_trait() + labs(tag = "B")
save_si(s9a + s9b + plot_layout(widths = c(1.4, 1)), "S9_Fig", 12, 4.8)

s10a <- plot_leave_one_out_pooled() + labs(tag = "A")
s10b <- plot_leave_one_out_modules() + labs(tag = "B")
save_si(s10a + s10b + plot_layout(widths = c(1, 1.2)), "S10_Fig", 12, 6)

s11a <- mm_gs_plot(hub_table, other) + labs(tag = "A")
s11b <- bubble_plot(hub_enrichment(enrichment, other, "tested"), hub_genes_title(other), n = 10, wrap = 38) + labs(tag = "B")
save_si(s11a + s11b + plot_layout(widths = c(1, 1.25)), "S11_Fig", 12, 6)

ensembl <- sub("^e([0-9]+)_.*", "\\1", basename(cfg$gprofiler_url))
stopifnot(grepl("^[0-9]+$", ensembl))

xl("S1_Table", list(Samples = tab("sample_characteristics"), DE_models = tab("DEGs_per_dataset")))

xl("S2_Table", list(
  Meta_analysis_all_genes = meta, Summary = tab("meta_analysis_summary"),
  Direction_concordance = tab("direction_concordance")
))

xl("S3_Table", list(Power_choice = tab("soft_threshold_choice"), All_powers = tab("soft_threshold_all_powers")))

xl("S4_Table", list(
  Module_sizes = tab("module_sizes"),
  Module_assignment = data.frame(gene_id = names(colors), symbol = gene_symbol(names(colors), data$symbol), module = colors)
))

xl("S5_Table", list(Per_dataset = mt, Pooled_random_effects = ms, Leave_one_out_pooled_r = tab("leave_one_out_pooled_r")))

xl("S6_Table", list(
  All_pairs = selection$preservation, Per_module = selection$per_module,
  Module_selection = selection$candidates
))

xl("S7_Table", list(
  Hub_genes_ranked = tab("hub_genes"), MM_GS_eligible_modules = tab("MM_GS_eligible_modules"),
  Hub_summary = tab("hub_summary"), Eigengene_correlation = tab("eigengene_correlation")
))

xl("S8_Table", enrichment_sheets())

xl("S10_Table", list(
  Merge_height = tab("merge_height_sensitivity"),
  Expression_network_main_modules = tab("expression_network_main_modules"),
  Expression_network_eligible = tab("expression_network_eligible"),
  Expression_network_powers = tab("expression_network_powers"),
  LOO_networks = tab("leave_one_out_networks"),
  LOO_main_modules = tab("leave_one_out_main_modules"),
  Power_sensitivity = tab("power_sensitivity")
))

xl("S11_Table", list(Settings = analysis_settings(ensembl), Software_versions = software_versions(ensembl)))
