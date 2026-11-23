source(file.path(snakemake@scriptdir, "functions.R"))

mt_main <- readRDS(snakemake@input[["module_trait"]])$mt
ref <- readRDS(snakemake@input[["hubs"]])
loo <- map(snakemake@input[["networks"]], readRDS)

pooled_without <- function(d) module_summary(mt_main, sets = setdiff(SETS, d)) %>% mutate(left_out = d)

pooled <- bind_rows(
  module_summary(mt_main) %>% mutate(left_out = "none"),
  map_dfr(SETS, pooled_without)
) %>%
  filter(module %in% names(ref$modules)) %>%
  select(left_out, module, direction, pooled_r, pooled_ci_low, pooled_ci_high, pooled_p, I2, n_sets_FDR, eligible)

write.csv(pooled, snakemake@output[["loo_pooled"]], row.names = FALSE)
write.csv(map_dfr(loo, "network"), snakemake@output[["loo_networks"]], row.names = FALSE)
write.csv(map_dfr(loo, "main_modules"), snakemake@output[["loo_main_modules"]], row.names = FALSE)
