if (length(snakemake@log) > 0) {
  log_con <- file(snakemake@log[[1]], open = "wt")
  sink(log_con)
  sink(log_con, type = "message")
}

suppressPackageStartupMessages({
  library(DESeq2)
  library(WGCNA)
  library(metafor)
  library(tidyverse)
})

select <- dplyr::select
filter <- dplyr::filter
count <- dplyr::count
rename <- dplyr::rename
options(stringsAsFactors = FALSE)

cfg <- snakemake@config
SETS <- names(cfg$datasets)
dataset_cfg <- function(field) sapply(cfg$datasets, `[[`, field)
SET_LABELS <- dataset_cfg("label")
ALL_SAME_DIRECTION <- paste("same direction in all", length(SETS))

trait_of <- function(samples) as.numeric(grepl("^TREATMENT", samples))

jaccard <- function(a, b) length(intersect(a, b)) / length(union(a, b))

faq_power <- function(n) case_when(n < 20 ~ 18, n < 30 ~ 16, n < 40 ~ 14, TRUE ~ 12)

fisher_z <- function(r) atanh(pmin(pmax(r, -0.9999), 0.9999))

stouffer_z <- function(stats, weights) as.numeric(stats %*% weights) / sqrt(sum(weights^2))

two_sided_p <- function(z) 2 * pnorm(-abs(z))

gene_symbol <- function(ids, symbol) {
  out <- unname(symbol[ids])
  out[is.na(out)] <- ids[is.na(out)]
  out
}

module_genes <- function(colors, modules) names(colors)[colors %in% modules]

dominant_module <- function(labels) {
  counts <- sort(table(labels[labels != "grey"]), decreasing = TRUE)
  if (length(counts)) names(counts)[1] else NA_character_
}

vst_matrix <- function(counts, coldata) {
  dds <- DESeqDataSetFromMatrix(counts, coldata[, "condition", drop = FALSE], ~condition)
  assay(vst(dds))
}

make_multiExpr <- function(counts, coldata, sets, genes) {
  vst_data <- map(set_names(sets), function(s) {
    vst_matrix(counts[[s]][intersect(genes, rownames(counts[[s]])), ], coldata[[s]])
  })
  common <- Reduce(intersect, map(vst_data, rownames))
  multiExpr <- map(vst_data, function(v) list(data = as.data.frame(t(v[common, ]))))
  gsg <- goodSamplesGenesMS(multiExpr, verbose = 0)
  if (!gsg$allOK) {
    multiExpr <- Map(function(set, good_samples) list(data = set$data[good_samples, gsg$goodGenes]), multiExpr, gsg$goodSamples)
  }
  multiExpr
}

scale_free_fit <- function(multiExpr) {
  map_dfr(names(multiExpr), function(s) {
    X <- as.matrix(multiExpr[[s]]$data)
    fit <- pickSoftThreshold(X,
      powerVector = 1:30, networkType = "signed", corFnc = "cor",
      corOptions = list(use = "p"), RsquaredCut = cfg$r2_cut, verbose = 0
    )$fitIndices
    data.frame(
      dataset = s, n = nrow(X), power = fit$Power, slope = fit$slope, unsigned_R2 = fit$SFT.R.sq,
      signed_R2 = -sign(fit$slope) * fit$SFT.R.sq, mean_k = fit$mean.k., median_k = fit$median.k.
    )
  })
}

smallest_power_above <- function(power, r2, cut) {
  passing <- power[which(r2 >= cut)]
  if (length(passing)) min(passing) else NA_real_
}

soft_threshold <- function(multiExpr) {
  all_powers <- scale_free_fit(multiExpr)
  choice <- all_powers %>%
    group_by(dataset, n) %>%
    summarise(
      smallest_power_signed_R2 = smallest_power_above(power, signed_R2, cfg$r2_cut),
      max_signed_R2 = max(signed_R2),
      power_at_max = power[which.max(signed_R2)],
      .groups = "drop"
    ) %>%
    mutate(
      faq_power = faq_power(n),
      chosen_power = ifelse(is.na(smallest_power_signed_R2), faq_power, smallest_power_signed_R2),
      reason = ifelse(is.na(smallest_power_signed_R2),
        paste("signed R2 never >=", cfg$r2_cut, ": WGCNA FAQ power"),
        paste("smallest power with signed R2 >=", cfg$r2_cut)
      )
    )
  powers <- setNames(choice$chosen_power, choice$dataset)[names(multiExpr)]
  choice <- choice %>%
    mutate(power_used = powers[dataset]) %>%
    left_join(
      all_powers %>% select(dataset, power, R2_at_used = signed_R2, mean_k_at_used = mean_k),
      by = c("dataset", "power_used" = "power")
    )
  list(all = all_powers, choice = choice, powers = powers)
}

build_network <- function(multiExpr, powers, threads) {
  tmp <- tempfile("network_")
  dir.create(tmp)
  old <- setwd(tmp)
  on.exit({
    setwd(old)
    unlink(tmp, recursive = TRUE)
  })
  blockwiseConsensusModules(multiExpr,
    power = unname(powers), networkType = "signed", corType = "pearson",
    corOptions = list(use = "p"), minModuleSize = cfg$min_module_size,
    deepSplit = cfg$deep_split, pamRespectsDendro = FALSE, mergeCutHeight = cfg$merge_cut,
    numericLabels = FALSE, maxBlockSize = 50000, saveIndividualTOMs = FALSE,
    nThreads = threads, verbose = 0
  )
}

network_colors <- function(net, multiExpr) setNames(as.character(net$colors), colnames(multiExpr[[1]]$data))

module_eigengenes <- function(multiExpr, colors) {
  genes <- colnames(multiExpr[[1]]$data)
  MEs <- multiSetMEs(multiExpr, universalColors = colors[genes], excludeGrey = TRUE, verbose = 0)
  setNames(map(MEs, function(x) as.matrix(x$data)), names(multiExpr))
}

perm_p <- function(x, y, n_perm) {
  r0 <- abs(stats::cor(x, y))
  rp <- abs(as.numeric(stats::cor(x, replicate(n_perm, sample(y)))))
  (sum(rp >= r0 - 1e-10) + 1) / (n_perm + 1)
}

module_trait <- function(multiExpr, colors, n_perm = 0) {
  MEs <- module_eigengenes(multiExpr, colors)
  map_dfr(names(multiExpr), function(s) {
    me <- MEs[[s]]
    y <- trait_of(rownames(multiExpr[[s]]$data))
    n <- length(y)
    r <- as.numeric(stats::cor(me, y))
    se <- 1 / sqrt(n - 3)
    data.frame(
      dataset = s, module = sub("^ME", "", colnames(me)), n = n, r = r, p = corPvalueStudent(r, n),
      ci_low = tanh(fisher_z(r) - 1.96 * se), ci_high = tanh(fisher_z(r) + 1.96 * se),
      perm_p = if (n_perm > 0) apply(me, 2, perm_p, y = y, n_perm = n_perm) else NA_real_
    )
  }) %>%
    mutate(fdr = p.adjust(p, "BH"), bonferroni = p.adjust(p, "bonferroni"))
}

pool_correlations <- function(r, n) {
  yi <- fisher_z(r)
  vi <- 1 / (n - 3)
  fit <- tryCatch(
    rma(yi = yi, vi = vi, method = "REML"),
    error = function(e) rma(yi = yi, vi = vi, method = "DL")
  )
  data.frame(
    pooled_r = tanh(fit$beta[1]), pooled_ci_low = tanh(fit$ci.lb), pooled_ci_high = tanh(fit$ci.ub),
    pooled_p = fit$pval, I2 = fit$I2, tau2 = fit$tau2
  )
}

module_summary <- function(mt, sets = SETS, min_sets = length(sets)) {
  mt <- mt %>%
    filter(dataset %in% sets) %>%
    mutate(strong = abs(r) >= cfg$r_cut)
  pooled <- mt %>%
    group_by(module) %>%
    summarise(pool_correlations(r, n), .groups = "drop") %>%
    mutate(pooled_fdr = p.adjust(pooled_p, "BH"))
  mt %>%
    group_by(module) %>%
    summarise(
      n_positive = sum(r > 0),
      n_negative = sum(r < 0),
      n_sets_nominal = sum(strong & p < cfg$fdr),
      n_sets_FDR = sum(strong & fdr < cfg$fdr),
      n_sets_bonferroni = sum(strong & bonferroni < cfg$fdr),
      n_sets_perm = if (all(is.na(perm_p))) NA_integer_ else sum(strong & perm_p < cfg$fdr),
      .groups = "drop"
    ) %>%
    left_join(pooled, by = "module") %>%
    mutate(
      same_sign_all = n_positive == length(sets) | n_negative == length(sets),
      eligible = same_sign_all & n_sets_FDR >= min_sets & pooled_fdr < cfg$fdr,
      eligible_bonferroni = same_sign_all & n_sets_bonferroni >= min_sets & pooled_fdr < cfg$fdr,
      direction = ifelse(pooled_r > 0, "positive", "negative")
    ) %>%
    arrange(pooled_p)
}

hub_calc <- function(multiExpr, colors, modules, directions) {
  MEs <- module_eigengenes(multiExpr, colors)
  genes <- module_genes(colors, modules[!is.na(modules)])
  map_dfr(names(multiExpr), function(s) {
    X <- as.matrix(multiExpr[[s]]$data[, genes, drop = FALSE])
    y <- trait_of(rownames(X))
    own_module <- match(paste0("ME", colors[genes]), colnames(MEs[[s]]))
    data.frame(
      gene_id = genes,
      MM = stats::cor(X, MEs[[s]])[cbind(seq_along(genes), own_module)],
      GS = as.numeric(stats::cor(X, y))
    )
  }) %>%
    group_by(gene_id) %>%
    summarise(consensusMM = mean(MM), consensusGS = mean(GS), .groups = "drop") %>%
    mutate(
      module = unname(colors[gene_id]),
      direction = unname(directions[module]),
      hub = consensusMM > cfg$mm_cut & consensusGS * direction > cfg$gs_cut
    )
}

best_match <- function(ref, new_colors) {
  sapply(names(ref$modules), function(m) {
    shared <- intersect(module_genes(ref$colors, m), names(new_colors))
    dominant_module(new_colors[shared])
  })
}

destination_summary <- function(labels, k = 4) {
  counts <- head(sort(table(labels[labels != "grey"]), decreasing = TRUE), k)
  if (length(counts)) paste0(names(counts), " (", as.integer(counts), ")", collapse = ", ") else ""
}

percent <- function(part, total = 1) round(100 * part / total, 1)

compare_modules <- function(ref, multiExpr, new_colors, new_ms) {
  matches <- best_match(ref, new_colors)
  matched <- unique(na.omit(matches))
  new_directions <- setNames(sign(new_ms$pooled_r[match(matched, new_ms$module)]), matched)
  new_hubs <- hub_calc(multiExpr, new_colors, matched, new_directions)
  map_dfr(names(ref$modules), function(m) {
    ref_genes <- module_genes(ref$colors, m)
    shared <- intersect(ref_genes, names(new_colors))
    new_labels <- new_colors[shared]
    best <- matches[[m]]
    best_genes <- module_genes(new_colors, best)
    best_ms <- new_ms[match(best, new_ms$module), ]
    ref_hubs <- ref$hubs[[m]]
    ref_sign <- ifelse(ref$modules[m] == "positive", 1, -1)
    best_hubs <- new_hubs$gene_id[new_hubs$hub & new_hubs$module %in% best]
    data.frame(
      reference_module = m,
      reference_direction = unname(ref$modules[m]),
      reference_size = length(ref_genes),
      in_new_network = length(shared),
      pct_in_new_network = percent(length(shared), length(ref_genes)),
      best_match = best,
      best_match_size = length(best_genes),
      pct_reference_in_best = percent(length(intersect(ref_genes, best_genes)), length(ref_genes)),
      jaccard = round(jaccard(ref_genes, best_genes), 3),
      pct_reference_to_grey = percent(sum(new_labels == "grey"), length(ref_genes)),
      best_pooled_r = best_ms$pooled_r,
      best_pooled_fdr = best_ms$pooled_fdr,
      best_eligible = isTRUE(best_ms$eligible),
      same_direction = isTRUE(sign(best_ms$pooled_r) == ref_sign),
      reference_hubs = length(ref_hubs),
      pct_hubs_in_best = percent(mean(new_colors[ref_hubs] %in% best)),
      pct_hubs_still_hub = percent(mean(ref_hubs %in% best_hubs)),
      destinations = destination_summary(new_labels)
    )
  })
}
