source(file.path(snakemake@scriptdir, "functions.R"))
library(gprofiler2)

data <- readRDS(snakemake@input[["data"]])
meta <- readRDS(snakemake@input[["meta"]])
selected <- readRDS(snakemake@input[["selection"]])$selected
hubs <- readRDS(snakemake@input[["hubs"]])
colors <- hubs$colors
modules <- hubs$modules
genes <- names(colors)

query_gprofiler <- function(genes, background) {
  gost(
    query = genes, organism = "btaurus", correction_method = "fdr", user_threshold = cfg$fdr,
    sources = c("GO:BP", "GO:MF", "GO:CC", "KEGG"), evcodes = TRUE, significant = TRUE,
    domain_scope = if (is.null(background)) "annotated" else "custom_annotated",
    custom_bg = background
  )
}

gost_with_retry <- function(genes, background, label, attempts = 6) {
  for (attempt in seq_len(attempts)) {
    result <- tryCatch(
      list(value = query_gprofiler(genes, background)),
      error = function(e) list(error = conditionMessage(e))
    )
    if (is.null(result$error)) {
      return(result$value)
    }
    message("g:Profiler attempt ", attempt, " failed: ", result$error)
    Sys.sleep(20 * attempt)
  }
  stop("g:Profiler failed for ", label)
}

run_gost <- function(genes, background, label, bg_label) {
  res <- gost_with_retry(genes, background, label)
  if (is.null(res)) {
    return(NULL)
  }
  res$result %>%
    select(-any_of(c("parents", "evidence_codes", "query", "significant"))) %>%
    mutate(
      query_set = label, background = bg_label,
      intersection_symbols = map_chr(strsplit(intersection, ","), ~ paste(gene_symbol(.x, data$symbol), collapse = ","))
    )
}

order_modules <- c(
  names(modules)[modules == "positive" & names(modules) %in% selected],
  names(modules)[modules == "negative" & names(modules) %in% selected]
)
order_modules <- c(order_modules, setdiff(names(modules), order_modules))

queries <- list()
for (m in order_modules) {
  label <- if (m %in% selected) paste0(modules[m], "_", m) else m
  queries[[paste0(label, "_hubs")]] <- hubs$hubs[[m]]
  queries[[paste0(label, "_module")]] <- genes[colors == m]
}

set_base_url(cfg$gprofiler_url)
enrichment <- map_dfr(names(queries), function(q) {
  bind_rows(
    run_gost(queries[[q]], NULL, q, "genome"),
    run_gost(queries[[q]], meta$gene_id, q, "tested")
  )
})

term_counts <- enrichment %>% count(query_set, background, source, name = "significant_terms")

saveRDS(enrichment, snakemake@output[["enrichment"]])
write.csv(enrichment, snakemake@output[["enrichment_csv"]], row.names = FALSE)
write.csv(term_counts, snakemake@output[["term_counts"]], row.names = FALSE)
