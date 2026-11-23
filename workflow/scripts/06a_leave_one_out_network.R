source(file.path(snakemake@scriptdir, "functions.R"))
allowWGCNAThreads(snakemake@threads)

data <- readRDS(snakemake@input[["data"]])
deseq2 <- readRDS(snakemake@input[["deseq2"]])
meta <- readRDS(snakemake@input[["meta"]])
ref <- readRDS(snakemake@input[["hubs"]])
left_out <- snakemake@wildcards[["dataset"]]
keep <- setdiff(SETS, left_out)
meta_genes <- meta$gene_id[meta$stouffer_sig]
n_samples <- sapply(data$coldata, nrow)

stouffer_genes <- function(sets) {
  wald <- Reduce(function(x, y) merge(x, y, by = "gene_id"), map(sets, function(s) {
    setNames(deseq2[[s]][, c("gene_id", "stat")], c("gene_id", s))
  }))
  wald <- wald[complete.cases(wald), ]
  z <- stouffer_z(as.matrix(wald[, sets]), sqrt(n_samples[sets]))
  wald$gene_id[p.adjust(two_sided_p(z), "BH") < cfg$fdr]
}

genes <- stouffer_genes(keep)
loo_multiExpr <- make_multiExpr(data$counts, data$coldata, keep, genes)
loo_sft <- soft_threshold(loo_multiExpr)
loo_net <- build_network(loo_multiExpr, loo_sft$powers, snakemake@threads)
loo_colors <- network_colors(loo_net, loo_multiExpr)
loo_ms <- module_summary(module_trait(loo_multiExpr, loo_colors), sets = keep)

main_modules <- compare_modules(ref, loo_multiExpr, loo_colors, loo_ms) %>%
  mutate(left_out = left_out, .before = 1)

network <- data.frame(
  left_out = left_out,
  input_genes = length(genes),
  shared_with_main_input = length(intersect(genes, meta_genes)),
  jaccard_with_main_input = round(jaccard(genes, meta_genes), 3),
  powers = paste(names(loo_sft$powers), loo_sft$powers, sep = "=", collapse = ", "),
  modules = length(setdiff(unique(loo_colors), "grey")),
  eligible_positive = sum(loo_ms$eligible & loo_ms$direction == "positive"),
  eligible_negative = sum(loo_ms$eligible & loo_ms$direction == "negative")
)

saveRDS(list(network = network, main_modules = main_modules), snakemake@output[[1]])
