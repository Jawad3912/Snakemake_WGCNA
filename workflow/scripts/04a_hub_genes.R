source(file.path(snakemake@scriptdir, "functions.R"))

data <- readRDS(snakemake@input[["data"]])
meta <- readRDS(snakemake@input[["meta"]])
multiExpr <- readRDS(snakemake@input[["multiExpr"]])
network <- readRDS(snakemake@input[["network"]])
ms <- readRDS(snakemake@input[["module_trait"]])$ms
colors <- network$colors
genes <- names(colors)

read_preservation <- function(file) {
  mp <- readRDS(file)
  ref <- sub("^preservation_(.*)\\.rds$", "\\1", basename(file))
  map_dfr(setdiff(SETS, ref), function(test) {
    Z <- mp$preservation$Z[[paste0("ref.", ref)]][[paste0("inColumnsAlsoPresentIn.", test)]]
    O <- mp$preservation$observed[[paste0("ref.", ref)]][[paste0("inColumnsAlsoPresentIn.", test)]]
    data.frame(
      reference = ref, test = test, module = rownames(Z), module_size = Z$moduleSize,
      Zsummary = Z$Zsummary.pres, Zdensity = Z$Zdensity.pres, Zconnectivity = Z$Zconnectivity.pres,
      medianRank = O[rownames(Z), "medianRank.pres"],
      medianRank_density = O[rownames(Z), "medianRankDensity.pres"],
      medianRank_connectivity = O[rownames(Z), "medianRankConnectivity.pres"]
    )
  })
}

membership_significance <- function(s) {
  X <- as.matrix(multiExpr[[s]]$data)
  y <- trait_of(rownames(X))
  n <- nrow(X)
  own_module <- match(paste0("ME", colors[genes]), colnames(ME[[s]]))
  ok <- !is.na(own_module)
  mm <- rep(NA_real_, length(genes))
  mm[ok] <- stats::cor(X[, ok, drop = FALSE], ME[[s]])[cbind(seq_len(sum(ok)), own_module[ok])]
  gs <- as.numeric(stats::cor(X, y))
  data.frame(gene_id = genes, dataset = s, MM = mm, MM_p = corPvalueStudent(mm, n), GS = gs, GS_p = corPvalueStudent(gs, n))
}

eigengene_correlation <- function(s) {
  cor(ME[[s]][, paste0("ME", names(modules)), drop = FALSE]) %>%
    as.table() %>%
    as.data.frame() %>%
    setNames(c("module_1", "module_2", "r")) %>%
    mutate(dataset = s)
}

preservation <- map_dfr(snakemake@input[["preservation"]], read_preservation) %>%
  mutate(
    preservation = case_when(Zsummary > cfg$min_zsummary ~ "strong", Zsummary > 2 ~ "moderate", TRUE ~ "none"),
    n_permutations = cfg$n_perm_preservation
  )

per_module <- preservation %>%
  filter(!module %in% c("grey", "gold")) %>%
  group_by(module, module_size) %>%
  summarise(
    pairs = n(), min_Zsummary = min(Zsummary), median_Zsummary = median(Zsummary),
    max_Zsummary = max(Zsummary), pairs_strong = sum(Zsummary > cfg$min_zsummary),
    pairs_moderate_or_strong = sum(Zsummary > 2), mean_medianRank = mean(medianRank),
    worst_medianRank = max(medianRank), .groups = "drop"
  ) %>%
  arrange(mean_medianRank)

candidates <- ms %>%
  filter(eligible) %>%
  left_join(per_module %>% select(module, min_Zsummary, median_Zsummary, pairs_strong, pairs, mean_medianRank), by = "module") %>%
  mutate(passes = min_Zsummary > cfg$min_zsummary) %>%
  arrange(direction, desc(passes), mean_medianRank) %>%
  group_by(direction) %>%
  mutate(selected = passes & row_number() == 1) %>%
  ungroup()

selected_in <- function(dir) {
  m <- candidates$module[candidates$selected & candidates$direction == dir]
  if (!length(m)) stop("No eligible ", dir, " module passes min_zsummary > ", cfg$min_zsummary)
  m
}
selected <- c(positive = selected_in("positive"), negative = selected_in("negative"))

candidates <- candidates %>%
  select(
    direction, module, selected, min_Zsummary, median_Zsummary, pairs_strong, pairs, mean_medianRank,
    pooled_r, pooled_ci_low, pooled_ci_high, pooled_fdr, I2, n_sets_FDR, n_sets_bonferroni, n_sets_perm
  )

eligible <- ms %>% filter(eligible)
modules <- setNames(eligible$direction, eligible$module)
module_direction <- setNames(sign(ms$pooled_r), ms$module)

ME <- module_eigengenes(multiExpr, colors)
mm_gs <- map_dfr(SETS, membership_significance) %>%
  pivot_wider(names_from = dataset, values_from = c(MM, MM_p, GS, GS_p), names_sep = "_")

MMm <- as.matrix(mm_gs[, paste0("MM_", SETS)])
GSm <- as.matrix(mm_gs[, paste0("GS_", SETS)])

hub_table <- mm_gs %>%
  mutate(
    symbol = gene_symbol(gene_id, data$symbol),
    module = unname(colors[gene_id]),
    direction = unname(modules[module]),
    selected_in_manuscript = module %in% selected,
    module_direction = unname(module_direction[module]),
    consensusMM = rowMeans(MMm),
    consensusGS = rowMeans(GSm),
    consensusMM_abs = rowMeans(abs(MMm)),
    consensusGS_abs = rowMeans(abs(GSm)),
    min_MM = apply(MMm, 1, min),
    n_sets_MM_above_cut = rowSums(MMm > cfg$mm_cut),
    n_sets_GS_in_direction = rowSums(sign(GSm) == module_direction),
    hub = consensusMM > cfg$mm_cut & consensusGS * module_direction > cfg$gs_cut,
    hub_original_definition = consensusMM_abs > cfg$mm_cut & consensusGS_abs > cfg$gs_cut,
    hub_consistent_all_datasets = hub & min_MM > 0 & n_sets_GS_in_direction == length(SETS)
  ) %>%
  filter(module %in% names(modules)) %>%
  left_join(
    meta %>% select(gene_id, starts_with("log2FC_"), starts_with("padj_"), stouffer_z, stouffer_fdr, fisher_fdr, direction_class),
    by = "gene_id"
  ) %>%
  group_by(module) %>%
  mutate(rank_in_module = rank(-(consensusMM + consensusGS * module_direction), ties.method = "first")) %>%
  ungroup() %>%
  relocate(gene_id, symbol, module, direction, selected_in_manuscript, hub, rank_in_module, consensusMM, consensusGS)

hub_genes <- hub_table %>%
  filter(hub) %>%
  arrange(direction, module, rank_in_module)

hub_summary <- hub_table %>%
  group_by(direction, module, selected_in_manuscript) %>%
  summarise(
    module_genes = n(), hubs = sum(hub), hubs_consistent_all_datasets = sum(hub_consistent_all_datasets),
    hubs_same_LFC_direction_all = sum(hub & direction_class == ALL_SAME_DIRECTION, na.rm = TRUE),
    .groups = "drop"
  )

eigengenes <- map_dfr(SETS, eigengene_correlation) %>%
  mutate(across(c(module_1, module_2), ~ sub("^ME", "", .))) %>%
  filter(module_1 < module_2)

hub_list <- map(set_names(names(modules)), ~ hub_table$gene_id[hub_table$hub & hub_table$module == .x])

genexplain <- hub_table %>%
  filter(module == selected[["positive"]], hub) %>%
  transmute(Gene = gene_id, Symbol = symbol)

saveRDS(
  list(preservation = preservation, per_module = per_module, candidates = candidates, selected = selected),
  snakemake@output[["selection"]]
)
saveRDS(list(table = hub_table, colors = colors, modules = modules, hubs = hub_list), snakemake@output[["hubs"]])
write.csv(preservation, snakemake@output[["all_pairs"]], row.names = FALSE)
write.csv(per_module, snakemake@output[["per_module"]], row.names = FALSE)
write.csv(candidates, snakemake@output[["candidates"]], row.names = FALSE)
write.csv(hub_table, snakemake@output[["mm_gs"]], row.names = FALSE)
write.csv(hub_genes, snakemake@output[["hub_genes"]], row.names = FALSE)
write.csv(hub_summary, snakemake@output[["hub_summary"]], row.names = FALSE)
write.csv(eigengenes, snakemake@output[["eigengenes"]], row.names = FALSE)
write.csv(genexplain, snakemake@output[["genexplain"]], row.names = FALSE)
