source(file.path(snakemake@scriptdir, "functions.R"))
allowWGCNAThreads(snakemake@threads)

multiExpr <- readRDS(snakemake@input[["multiExpr"]])
powers <- readRDS(snakemake@input[["powers"]])
ref <- readRDS(snakemake@input[["hubs"]])

power_label <- function(p) paste(names(p), p, sep = "=", collapse = ", ")

alternative <- dataset_cfg("alternative_power")
alt_net <- build_network(multiExpr, alternative, snakemake@threads)
alt_colors <- network_colors(alt_net, multiExpr)
alt_ms <- module_summary(module_trait(multiExpr, alt_colors))

power_sensitivity <- compare_modules(ref, multiExpr, alt_colors, alt_ms) %>%
  mutate(powers_main = power_label(powers), powers_alternative = power_label(alternative), .before = 1)

write.csv(power_sensitivity, snakemake@output[[1]], row.names = FALSE)
