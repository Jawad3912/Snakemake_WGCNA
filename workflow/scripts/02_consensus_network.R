source(file.path(snakemake@scriptdir, "functions.R"))
allowWGCNAThreads(snakemake@threads)

data <- readRDS(snakemake@input[["data"]])
meta <- readRDS(snakemake@input[["meta"]])
reference <- read.csv(snakemake@input[["reference"]])

meta_genes <- meta$gene_id[meta$stouffer_sig]
multiExpr <- make_multiExpr(data$counts, data$coldata, SETS, meta_genes)
genes <- colnames(multiExpr[[1]]$data)

sft <- soft_threshold(multiExpr)
powers <- sft$powers

net <- build_network(multiExpr, powers, snakemake@threads)
if (length(net$dendrograms) > 1) stop("Network was split into several blocks")

reference <- reference[match(genes, reference$gene_id), ]
match_reference_labels <- function(labels, reference_labels) {
  setNames(as.character(matchLabels(as.character(labels), reference_labels)), genes)
}
colors <- match_reference_labels(net$colors, reference$module)
unmerged <- match_reference_labels(net$unmergedColors, reference$module_unmerged)

network <- list(
  colors = colors, unmerged = unmerged, dendro = net$dendrograms[[1]],
  blockGenes = net$blockGenes[[1]], powers = powers
)
sizes <- sort(table(colors), decreasing = TRUE)

set.seed(cfg$seed)
mt <- module_trait(multiExpr, colors, n_perm = cfg$n_perm_trait)
ms <- module_summary(mt)

saveRDS(multiExpr, snakemake@output[["multiExpr"]])
saveRDS(powers, snakemake@output[["powers"]])
saveRDS(network, snakemake@output[["network"]])
saveRDS(list(mt = mt, ms = ms), snakemake@output[["module_trait"]])
write.csv(sft$all, snakemake@output[["all_powers"]], row.names = FALSE)
write.csv(sft$choice, snakemake@output[["choice"]], row.names = FALSE)
write.csv(data.frame(module = names(sizes), n_genes = as.integer(sizes)), snakemake@output[["sizes"]], row.names = FALSE)
write.csv(data.frame(gene_id = genes, module = colors, module_unmerged = unmerged), snakemake@output[["modules"]], row.names = FALSE)
write.csv(mt, snakemake@output[["per_dataset"]], row.names = FALSE)
write.csv(ms, snakemake@output[["pooled"]], row.names = FALSE)
