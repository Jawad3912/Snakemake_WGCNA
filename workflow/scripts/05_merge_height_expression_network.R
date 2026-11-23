source(file.path(snakemake@scriptdir, "functions.R"))
allowWGCNAThreads(snakemake@threads)

data <- readRDS(snakemake@input[["data"]])
meta <- readRDS(snakemake@input[["meta"]])
multiExpr <- readRDS(snakemake@input[["multiExpr"]])
network <- readRDS(snakemake@input[["network"]])
selected <- readRDS(snakemake@input[["selection"]])$selected
ref <- readRDS(snakemake@input[["hubs"]])
colors <- network$colors
genes <- names(colors)
meta_genes <- meta$gene_id[meta$stouffer_sig]

follow_module <- function(new_colors, m) {
  tb <- table(new_colors[colors == m])
  best <- names(tb)[which.max(tb)]
  list(module = best, size = sum(new_colors == best), jaccard = round(jaccard(genes[colors == m], genes[new_colors == best]), 2))
}

merge_at_height <- function(h) {
  merged <- mergeCloseModules(multiExpr, network$unmerged, cutHeight = h, consensusQuantile = 0, relabel = FALSE, verbose = 0)
  new_colors <- setNames(as.character(matchLabels(as.character(merged$colors), colors)), genes)
  new_ms <- module_summary(module_trait(multiExpr, new_colors))
  pos <- follow_module(new_colors, selected[["positive"]])
  neg <- follow_module(new_colors, selected[["negative"]])
  data.frame(
    mergeCutHeight = h, n_modules = length(setdiff(unique(new_colors), "grey")),
    n_eligible_positive = sum(new_ms$eligible & new_ms$direction == "positive"),
    n_eligible_negative = sum(new_ms$eligible & new_ms$direction == "negative"),
    positive_becomes = pos$module, positive_size = pos$size, positive_jaccard = pos$jaccard,
    negative_becomes = neg$module, negative_size = neg$size, negative_jaccard = neg$jaccard
  )
}

expressed_genes_in <- function(s) {
  m <- data$counts[[s]]
  n_min <- min(table(data$coldata[[s]]$condition))
  rownames(m)[rowSums(m >= cfg$expression_min_count) >= n_min]
}

main_module_overlap <- function(m) {
  members <- module_genes(expr_colors, m)
  main_labels <- ref$colors[intersect(members, names(ref$colors))]
  best <- dominant_module(main_labels)
  tibble(
    size = length(members),
    in_meta_input = sum(members %in% meta_genes),
    best_main_module = best,
    pct_in_best_main = percent(sum(main_labels == best), length(members))
  )
}

merge_height <- map_dfr(unlist(cfg$merge_heights), merge_at_height)

expressed_genes <- Reduce(intersect, map(SETS, expressed_genes_in))
expr_multiExpr <- make_multiExpr(data$counts, data$coldata, SETS, expressed_genes)
expr_sft <- soft_threshold(expr_multiExpr)
expr_net <- build_network(expr_multiExpr, expr_sft$powers, snakemake@threads)
expr_colors <- network_colors(expr_net, expr_multiExpr)
expr_mt <- module_trait(expr_multiExpr, expr_colors)
expr_ms <- module_summary(expr_mt)

main_modules <- compare_modules(ref, expr_multiExpr, expr_colors, expr_ms)

expr_eligible <- expr_ms %>%
  filter(eligible) %>%
  mutate(map_dfr(module, main_module_overlap)) %>%
  mutate(pct_outside_meta_input = percent(size - in_meta_input, size)) %>%
  select(
    module, direction, size, in_meta_input, pct_outside_meta_input, best_main_module, pct_in_best_main,
    pooled_r, pooled_ci_low, pooled_ci_high, pooled_fdr, n_sets_FDR, n_sets_bonferroni
  )

expr_powers <- expr_sft$choice %>%
  transmute(dataset, n, chosen_power, reason, R2_at_chosen = R2_at_used, mean_k_at_chosen = mean_k_at_used)

write.csv(merge_height, snakemake@output[["merge_height"]], row.names = FALSE)
write.csv(main_modules, snakemake@output[["main_modules"]], row.names = FALSE)
write.csv(expr_eligible, snakemake@output[["eligible"]], row.names = FALSE)
write.csv(expr_powers, snakemake@output[["powers"]], row.names = FALSE)
write.csv(expr_mt, snakemake@output[["module_trait"]], row.names = FALSE)
write.csv(expr_ms, snakemake@output[["pooled"]], row.names = FALSE)
write.csv(data.frame(gene_id = names(expr_colors), module = expr_colors), snakemake@output[["modules"]], row.names = FALSE)
