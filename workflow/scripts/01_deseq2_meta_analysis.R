source(file.path(snakemake@scriptdir, "functions.R"))

samples <- read.csv(snakemake@input[["samples"]], na.strings = c("", "NA"))
count_tables <- map(setNames(snakemake@input[["counts"]], SETS), read.delim, check.names = FALSE)

count_matrix <- function(x) {
  m <- round(as.matrix(x[, -(1:2)]))
  rownames(m) <- x$gene_id
  m
}

gene_symbols <- function(count_tables) {
  map(count_tables, ~ select(.x, gene_id, gene_name)) %>%
    bind_rows() %>%
    distinct(gene_id, .keep_all = TRUE) %>%
    mutate(gene_name = ifelse(is.na(gene_name) | gene_name == "", gene_id, gene_name)) %>%
    with(setNames(gene_name, gene_id))
}

sample_coldata <- function(s, sample_ids) {
  info <- samples %>% filter(dataset == s)
  info <- info[match(sample_ids, info$sample), ]
  stopifnot(all(info$sample == sample_ids))
  cd <- info %>%
    transmute(sample, condition = factor(condition, levels = c("Control", "Mastitis")), cow, parity, lactation, farm) %>%
    select(sample, condition, where(~ !anyNA(.x))) %>%
    mutate(across(-c(sample, condition), ~ factor(as.character(.x)))) %>%
    column_to_rownames("sample")
  stopifnot(all(trait_of(info$sample) == (cd$condition == "Mastitis")))
  cd
}

design_of <- function(s) as.formula(cfg$datasets[[s]]$design)

run_deseq2 <- function(s) {
  dds <- DESeqDataSetFromMatrix(counts[[s]], coldata[[s]], design_of(s))
  res <- results(DESeq(dds, quiet = TRUE), contrast = c("condition", "Mastitis", "Control"))
  data.frame(gene_id = rownames(res), as.data.frame(res), row.names = NULL)
}

de_summary <- function(s) {
  deseq2[[s]] %>%
    summarise(
      dataset = s,
      n_control = sum(coldata[[s]]$condition == "Control"),
      n_mastitis = sum(coldata[[s]]$condition == "Mastitis"),
      design = paste(deparse(design_of(s)), collapse = ""),
      genes_tested = sum(!is.na(padj)),
      "DEGs_padj_{cfg$fdr}" := sum(padj < cfg$fdr, na.rm = TRUE),
      "DEGs_padj_{cfg$fdr}_absLFC_1" := sum(padj < cfg$fdr & abs(log2FoldChange) >= 1, na.rm = TRUE),
      up = sum(padj < cfg$fdr & log2FoldChange >= 1, na.rm = TRUE),
      down = sum(padj < cfg$fdr & log2FoldChange <= -1, na.rm = TRUE)
    )
}

per_dataset_columns <- function(s) {
  deseq2[[s]] %>%
    select(gene_id, log2FC = log2FoldChange, SE = lfcSE, stat, p = pvalue, padj) %>%
    rename_with(~ paste0(.x, "_", s), -gene_id)
}

stat_matrix <- function(column, genes) {
  sapply(SETS, function(s) deseq2[[s]][[column]][match(genes, deseq2[[s]]$gene_id)])
}

counts <- map(count_tables, count_matrix)
symbol <- gene_symbols(count_tables)
coldata <- map(set_names(SETS), function(s) sample_coldata(s, colnames(counts[[s]])))
n_samples <- sapply(coldata, nrow)

deseq2 <- map(set_names(SETS), run_deseq2)
degs <- map_dfr(SETS, de_summary)

genes <- Reduce(intersect, map(deseq2, function(d) d$gene_id[!is.na(d$pvalue)]))
pvalues <- stat_matrix("pvalue", genes)
wald <- stat_matrix("stat", genes)

fisher_chi2 <- -2 * rowSums(log(pvalues))
fisher_p <- pchisq(fisher_chi2, 2 * length(SETS), lower.tail = FALSE)
stouffer <- stouffer_z(wald, sqrt(n_samples[SETS]))
stouffer_p <- two_sided_p(stouffer)

meta <- data.frame(
  gene_id = genes, stouffer_z = stouffer, stouffer_p = stouffer_p,
  stouffer_fdr = p.adjust(stouffer_p, "BH"), fisher_chi2 = fisher_chi2,
  fisher_p = fisher_p, fisher_fdr = p.adjust(fisher_p, "BH")
) %>%
  left_join(reduce(map(SETS, per_dataset_columns), inner_join, by = "gene_id"), by = "gene_id")

lfc <- as.matrix(meta[, paste0("log2FC_", SETS)])
pv <- as.matrix(meta[, paste0("p_", SETS)])

meta <- meta %>%
  mutate(
    symbol = gene_symbol(gene_id, symbol),
    n_up = rowSums(lfc > 0, na.rm = TRUE),
    n_down = rowSums(lfc < 0, na.rm = TRUE),
    direction_class = case_when(
      rowSums(is.na(lfc)) > 0 ~ "missing in a dataset",
      pmax(n_up, n_down) == length(SETS) ~ ALL_SAME_DIRECTION,
      TRUE ~ paste(pmax(n_up, n_down), "of", length(SETS), "agree")
    ),
    significant_conflict = rowSums(lfc > 0 & pv < cfg$fdr, na.rm = TRUE) > 0 & rowSums(lfc < 0 & pv < cfg$fdr, na.rm = TRUE) > 0,
    meta_direction = ifelse(stouffer_z > 0, "up", "down"),
    fisher_sig = fisher_fdr < cfg$fdr,
    stouffer_sig = stouffer_fdr < cfg$fdr
  ) %>%
  relocate(gene_id, symbol)

fisher_genes <- meta$gene_id[meta$fisher_sig]
stouffer_genes <- meta$gene_id[meta$stouffer_sig]
meta_summary <- data.frame(
  gene_set = c("Fisher", "Stouffer", "Fisher and Stouffer", "Fisher only", "Stouffer only"),
  n_genes = c(
    length(fisher_genes), length(stouffer_genes), length(intersect(fisher_genes, stouffer_genes)),
    length(setdiff(fisher_genes, stouffer_genes)), length(setdiff(stouffer_genes, fisher_genes))
  )
)

direction <- bind_rows(Fisher = filter(meta, fisher_sig), Stouffer = filter(meta, stouffer_sig), .id = "gene_set") %>%
  count(gene_set, direction_class, significant_conflict) %>%
  group_by(gene_set) %>%
  mutate(percent = round(100 * n / sum(n), 1)) %>%
  ungroup()

saveRDS(list(counts = counts, coldata = coldata, symbol = symbol), snakemake@output[["data"]])
saveRDS(deseq2, snakemake@output[["deseq2"]])
saveRDS(meta, snakemake@output[["meta_rds"]])
write.csv(samples, snakemake@output[["samples"]], row.names = FALSE)
write.csv(degs, snakemake@output[["degs"]], row.names = FALSE)
write.csv(meta, snakemake@output[["meta"]], row.names = FALSE)
write.csv(meta_summary, snakemake@output[["meta_summary"]], row.names = FALSE)
write.csv(direction, snakemake@output[["direction"]], row.names = FALSE)
