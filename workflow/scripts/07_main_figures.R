source(file.path(snakemake@scriptdir, "functions.R"))
source(file.path(snakemake@scriptdir, "plot_functions.R"))
library(pheatmap)
library(png)
library(grid)

data <- readRDS(snakemake@input[["data"]])
multiExpr <- readRDS(snakemake@input[["multiExpr"]])
powers <- readRDS(snakemake@input[["powers"]])
network <- readRDS(snakemake@input[["network"]])
module_trait <- readRDS(snakemake@input[["module_trait"]])
selected <- readRDS(snakemake@input[["selection"]])$selected
hub_table <- readRDS(snakemake@input[["hubs"]])$table
enrichment <- readRDS(snakemake@input[["enrichment"]])
mt <- module_trait$mt
ms <- module_trait$ms
colors <- network$colors
genes <- names(colors)
pos <- selected[["positive"]]
neg <- selected[["negative"]]

panel_files <- output_paths("panels")
pdf_files <- output_paths("pdf")
jpg_files <- output_paths("jpg")
has_species <- grepl("E\\. coli|S\\. aureus", DS)
DS_EXPR <- as.expression(lapply(seq_along(DS), function(i) if (has_species[i]) italic_species(DS[i]) else DS[i]))
DS_PARSE <- setNames(ifelse(has_species, sapply(DS, function(x) paste(deparse(italic_species(x)), collapse = "")), sprintf("'%s'", DS)), DS)
th <- theme_bw(base_size = 12) + theme(plot.title = element_text(face = "bold", hjust = 0.5))
heat_palette <- colorRampPalette(c("#2166ac", "white", "#b2182b"))(100)

save_png <- function(name, w, h, expr) {
  f <- panel_files[[name]]
  png(f, width = w, height = h, units = "in", res = 300)
  expr
  dev.off()
  wrap_elements(full = rasterGrob(readPNG(f), interpolate = TRUE))
}

add_tags <- function(p) p + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "bold", size = 20))

save_figure <- function(p, name, w, h, jpg_from_pdf = TRUE) {
  pdf_file <- pdf_files[[name]]
  jpg_file <- jpg_files[[name]]
  ggsave(pdf_file, p, width = w, height = h, device = cairo_pdf, bg = "white")
  if (jpg_from_pdf) {
    system2("pdftoppm", c(
      "-jpeg", "-jpegopt", "quality=95", "-singlefile", "-scale-to-x", w * 300, "-scale-to-y", h * 300,
      shQuote(pdf_file), shQuote(tools::file_path_sans_ext(jpg_file))
    ))
  } else {
    ggsave(jpg_file, p, width = w, height = h, dpi = 300, quality = 95, bg = "white")
  }
}

combined_vst <- function() {
  vst_all <- map(SETS, function(s) {
    cnt <- data$counts[[s]][rowSums(data$counts[[s]] >= cfg$expression_min_count) >= 3, ]
    vst_matrix(cnt, data$coldata[[s]])
  })
  n_per_set <- map_int(vst_all, ncol)
  common <- Reduce(intersect, map(vst_all, rownames))
  expr <- do.call(cbind, map(vst_all, function(v) v[common, ]))
  colnames(expr) <- paste0(rep(SETS, n_per_set), "_", colnames(expr))
  ann <- data.frame(
    Dataset = factor(rep(DS, n_per_set), DS),
    Condition = unlist(map(SETS, function(s) as.character(data$coldata[[s]]$condition))),
    row.names = colnames(expr)
  )
  top <- head(order(apply(expr, 1, var), decreasing = TRUE), 1000)
  list(expr = expr[top, ], ann = ann)
}

plot_pca <- function(expr, ann) {
  pca <- prcomp(t(expr))
  ve <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
  pca_data <- data.frame(pca$x[, 1:2], ann) %>%
    mutate(point_fill = ifelse(Condition == "Mastitis", dcol[as.character(Dataset)], "white"))
  ggplot(pca_data, aes(PC1, PC2, colour = Dataset, shape = Condition, fill = point_fill)) +
    geom_point(size = 3.2, stroke = 1.1) +
    scale_colour_manual(values = dcol, breaks = DS, labels = DS_EXPR) +
    scale_shape_manual(values = c(Control = 21, Mastitis = 22)) +
    scale_fill_identity() +
    guides(
      colour = guide_legend(order = 1, override.aes = list(shape = 15, size = 4)),
      shape = guide_legend(order = 2, override.aes = list(fill = c("white", "grey30"), colour = "grey30", size = 3.5))
    ) +
    labs(x = paste0("PC1 (", ve[1], "% of variance)"), y = paste0("PC2 (", ve[2], "% of variance)")) +
    th +
    theme(legend.text = element_text(hjust = 0))
}

plot_sample_correlation <- function(expr, ann) {
  heatmap <- pheatmap(cor(expr),
    annotation_row = ann, annotation_col = ann,
    annotation_colors = list(Dataset = dcol, Condition = c(Control = "grey80", Mastitis = "black")),
    show_rownames = FALSE, show_colnames = FALSE, annotation_names_row = FALSE,
    color = colorRampPalette(c("blue", "white", "red"))(100), border_color = NA,
    main = "Sample-to-sample Pearson correlation (r)", fontsize = 11,
    treeheight_row = 30, treeheight_col = 30, silent = TRUE
  )
  wrap_elements(full = heatmap$gtable) + theme(plot.margin = margin(t = 25, l = 15))
}

plot_soft_threshold <- function() {
  pt <- read.csv(snakemake@input[["all_powers"]]) %>%
    mutate(set = dataset, dataset = factor(SET_LABELS[dataset], DS), dsp = factor(DS_PARSE[as.character(dataset)], DS_PARSE))
  chosen <- pt %>% filter(power == powers[set])
  ggplot(pt, aes(power, signed_R2)) +
    geom_hline(yintercept = cfg$r2_cut, linetype = 2, colour = "red") +
    geom_line() +
    geom_point(size = 0.9) +
    geom_point(data = chosen, colour = "red", size = 3.2, shape = 21, stroke = 1.2) +
    facet_wrap(~dsp, ncol = 1, labeller = label_parsed) +
    scale_y_continuous(breaks = seq(-0.5, 1, 0.25)) +
    labs(x = expression("Soft-thresholding power (" * beta * ")"), y = expression("Scale-free topology fit (" * R^2 * ", signed)")) +
    th
}

plot_gene_dendrogram <- function() {
  block <- genes[network$blockGenes]
  plotDendroAndColors(network$dendro, cbind(network$unmerged[block], colors[block]),
    c("Dynamic Tree Cut", "Merged modules"),
    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE,
    guideHang = 0.05, main = "Consensus gene dendrogram and module colors", marAll = c(1, 9, 3, 1),
    ylab = "Height (1 - TOM)", cex.colorLabels = 1, cex.main = 1.3, cex.lab = 1.1, cex.axis = 1
  )
}

mean_eigengene_correlation <- function() {
  ME <- module_eigengenes(multiExpr, colors)
  mods <- sub("^ME", "", colnames(ME[[1]]))
  C <- Reduce(`+`, map(ME, cor)) / length(ME)
  dimnames(C) <- list(mods, mods)
  C
}

plot_eigengene_dendrogram <- function(hc) {
  par(mar = c(1, 5, 3, 1))
  plot(hc,
    labels = paste0("ME", hc$labels), xlab = "", sub = "", main = "Consensus module eigengene dendrogram",
    ylab = "Height (1 - |r|)", hang = -1, cex = 0.85, cex.lab = 1.1, cex.main = 1.2
  )
}

plot_eigengene_heatmap <- function(C, hc) {
  o <- hc$order
  mods <- rownames(C)[o]
  n <- length(o)
  layout(matrix(1:2, 1), widths = c(5.5, 1))
  par(mar = c(2, 2, 3, 1))
  image(1:n, 1:n, C[o, rev(o)],
    col = heat_palette, zlim = c(-1, 1), axes = FALSE, xlab = "", ylab = "",
    main = paste0("Eigengene correlation (mean of ", c("two", "three", "four", "five")[length(SETS) - 1], " datasets)")
  )
  for (i in seq_along(o)) {
    rect(i - 0.5, -0.6, i + 0.5, 0.45, col = mods[i], border = NA, xpd = TRUE)
    rect(-0.6, n + 1 - i - 0.5, 0.45, n + 1 - i + 0.5, col = mods[i], border = NA, xpd = TRUE)
  }
  par(mar = c(8, 1, 8, 4))
  image(1, seq(-1, 1, length.out = 100), t(seq(-1, 1, length.out = 100)), col = heat_palette, axes = FALSE, xlab = "", ylab = "")
  axis(4, las = 1)
  mtext("Pearson r", side = 3, line = 0.8, cex = 0.9)
}

plot_module_sizes <- function() {
  sizes <- sort(table(colors), decreasing = TRUE)
  ggplot(data.frame(module = names(sizes), n = as.integer(sizes)), aes(reorder(module, -n), n, fill = module)) +
    geom_col(colour = "black", linewidth = 0.2) +
    geom_text(aes(label = n), vjust = -0.3, size = 3) +
    scale_fill_identity() +
    labs(x = "Module", y = "Number of genes") +
    th +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
}

plot_module_trait_heatmap <- function() {
  heat <- mt %>%
    mutate(
      dataset = factor(SET_LABELS[dataset], DS),
      module = factor(module, ms$module[order(ms$pooled_r)]),
      star = ifelse(abs(r) >= cfg$r_cut & fdr < cfg$fdr, "*", ""),
      lab = sprintf("%.2f (%s)", r, formatC(fdr, format = "e", digits = 1))
    )
  ggplot(heat, aes(dataset, module, fill = r)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = lab), size = 3) +
    geom_text(aes(label = star), size = 6, nudge_x = 0.42, nudge_y = -0.15) +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0, limits = c(-1, 1), name = "Pearson r") +
    scale_x_discrete(labels = DS_EXPR, position = "bottom") +
    labs(x = NULL, y = "Module") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(face = "bold"), axis.text.y = element_text(face = "bold"), panel.grid = element_blank())
}

plot_hub_enrichment <- function(query, m, k) {
  bubble_plot(hub_enrichment(enrichment, query, "tested"), hub_genes_title(m), n = 10, wrap = 46, k = k, pub = TRUE)
}

vst_top <- combined_vst()
fig1 <- add_tags(plot_pca(vst_top$expr, vst_top$ann) + plot_sample_correlation(vst_top$expr, vst_top$ann) +
  plot_layout(widths = c(1, 1.05)))
save_figure(fig1, "Fig1", 13, 6)

C <- mean_eigengene_correlation()
hc <- hclust(as.dist(1 - abs(C)), method = "average")
fig2a <- plot_soft_threshold()
fig2b <- save_png("Fig2B", 10, 5.5, plot_gene_dendrogram())
fig2c_dendro <- save_png("Fig2C_dendrogram", 7.5, 4.5, plot_eigengene_dendrogram(hc))
fig2c_heatmap <- save_png("Fig2C_heatmap", 7.5, 6.5, plot_eigengene_heatmap(C, hc))
fig2c <- wrap_elements(full = fig2c_dendro / fig2c_heatmap + plot_layout(heights = c(1, 1.1)))
fig2d <- plot_module_sizes()
fig2 <- add_tags((fig2a + fig2b + plot_layout(widths = c(1, 2.2))) /
  (fig2c + fig2d + plot_layout(widths = c(0.85, 1.4))) +
  plot_layout(heights = c(1, 1.1)))
save_figure(fig2, "Fig2", 14, 13)

save_figure(plot_module_trait_heatmap(), "Fig3", 7, 10)

k <- 1.25
fig4 <- wrap_elements(full = (mm_gs_plot(hub_table, pos, k) + labs(tag = "A")) + (mm_gs_plot(hub_table, neg, k) + labs(tag = "B"))) /
  ((plot_hub_enrichment(paste0("positive_", pos), pos, k) + labs(tag = "C")) +
    (plot_hub_enrichment(paste0("negative_", neg), neg, k) + labs(tag = "D"))) +
  plot_layout(heights = c(1, 2.9))
save_figure(fig4, "Fig4", 14, 17.5, jpg_from_pdf = FALSE)
