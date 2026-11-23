source(file.path(snakemake@scriptdir, "functions.R"))
allowWGCNAThreads(snakemake@threads)

multiExpr <- readRDS(snakemake@input[["multiExpr"]])
network <- readRDS(snakemake@input[["network"]])
reference <- snakemake@wildcards[["dataset"]]

colors <- network$colors
colorList <- setNames(rep(list(colors), length(SETS)), SETS)
max_size <- max(table(colors[colors != "grey"]))

mp <- modulePreservation(multiExpr, colorList,
  referenceNetworks = match(reference, SETS),
  nPermutations = cfg$n_perm_preservation, networkType = "signed",
  maxModuleSize = max_size, randomSeed = cfg$seed, quickCor = 0,
  savePermutedStatistics = FALSE, verbose = 0
)

saveRDS(mp, snakemake@output[[1]])
