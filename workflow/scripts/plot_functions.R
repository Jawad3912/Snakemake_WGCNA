library(patchwork)

DS <- unname(SET_LABELS)
dcol <- setNames(c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e")[seq_along(DS)], DS)

cap <- function(x) paste0(toupper(substring(x, 1, 1)), substring(x, 2))

italic_species <- function(s, bold = FALSE) {
  it <- if (bold) "bolditalic" else "italic"
  pl <- if (bold) "bold" else "plain"
  pieces <- regmatches(s, gregexpr("E\\. coli|S\\. aureus", s), invert = NA)[[1]]
  out <- character()
  for (i in seq_along(pieces)) {
    q <- trimws(pieces[i])
    if (nzchar(q)) out <- c(out, sprintf('%s("%s")', if (i %% 2 == 0) it else pl, q))
  }
  parse(text = paste(out, collapse = "~"))[[1]]
}

italic_labels <- function(x) as.expression(lapply(x, italic_species))

hub_genes_title <- function(m) paste(cap(m), "module hub genes")

hub_enrichment <- function(enrichment, query, bg) {
  enrichment %>% filter(query_set == paste0(query, "_hubs"), background == bg)
}

output_paths <- function(key) {
  paths <- snakemake@output[[key]]
  setNames(paths, tools::file_path_sans_ext(basename(paths)))
}

theme_figure <- function(k = 1) {
  theme_bw(base_size = 11 * k) +
    theme(
      axis.text = element_text(colour = "black", size = 9.5 * k),
      axis.title = element_text(colour = "black", size = 11 * k),
      strip.text = element_text(face = "bold", size = 10 * k),
      strip.background = element_rect(fill = "grey90"),
      plot.tag = element_text(face = "bold", size = 16 * k),
      plot.title = element_text(face = "bold", size = 11.5 * k, hjust = 0.5),
      legend.text = element_text(size = 9.5 * k),
      legend.title = element_text(size = 10 * k),
      panel.grid.minor = element_blank()
    )
}

bubble_plot <- function(d, title, n = 8, wrap = 42, k = 1, pub = FALSE) {
  d <- d %>%
    mutate(source = sub("GO:", "", source)) %>%
    group_by(source) %>%
    arrange(p_value, .by_group = TRUE) %>%
    slice_head(n = n) %>%
    ungroup() %>%
    mutate(
      rich = intersection_size / term_size, lp = -log10(p_value),
      source = factor(source, c("BP", "CC", "MF", "KEGG")), term = str_wrap(cap(term_name), wrap)
    ) %>%
    arrange(source, rich) %>%
    mutate(term = factor(term, unique(term)))
  p <- ggplot(d, aes(rich, term, size = intersection_size, colour = lp)) +
    geom_point(alpha = 0.9) +
    facet_grid(source ~ ., scales = "free_y", space = "free_y") +
    scale_colour_gradient(low = "#00BFC4", high = "#F8696B", name = expression(-log[10] * "(FDR)")) +
    scale_size_continuous(range = c(2, 7), name = "Gene count") +
    scale_x_continuous(expand = expansion(mult = 0.12), breaks = scales::pretty_breaks(4)) +
    labs(x = "Rich factor", y = NULL, title = title) +
    theme_figure(k) +
    theme(
      axis.text.y = element_text(size = 9.3 * k, lineheight = 0.8, margin = margin(r = 3)),
      strip.text.y = element_text(angle = -90)
    )
  if (pub) {
    p <- p +
      guides(
        colour = guide_colourbar(title.position = "top", barwidth = unit(4, "cm"), barheight = unit(0.45, "cm"), order = 1),
        size = guide_legend(title.position = "top", nrow = 1, order = 2)
      ) +
      scale_size_continuous(range = c(2, 7), name = "Gene count", breaks = function(l) {
        b <- pretty(l, 3)
        b[b >= l[1] & b <= l[2]]
      }) +
      theme(
        axis.text.y = element_text(size = 8.6 * k, face = "bold", lineheight = 0.85, margin = margin(r = 3)),
        axis.text.x = element_text(face = "bold"), axis.title.x = element_text(face = "bold"),
        strip.text.y = element_text(angle = -90, size = 9.5 * k), strip.background = element_rect(fill = "grey85"),
        legend.position = "bottom", legend.box = "horizontal",
        legend.title = element_text(face = "bold", size = 9.5 * k, vjust = 0.8),
        legend.text = element_text(face = "bold", size = 8.6 * k), legend.spacing.x = unit(0.6, "cm"),
        panel.spacing = unit(0.25, "cm"), plot.margin = margin(5, 20, 5, 5)
      )
  }
  p
}

mm_gs_plot <- function(hub_table, m, k = 1) {
  d <- hub_table %>%
    filter(module == m) %>%
    mutate(GSd = consensusGS * module_direction, hub = consensusMM > cfg$mm_cut & GSd > cfg$gs_cut)
  ggplot(d, aes(consensusMM, GSd)) +
    geom_point(aes(colour = hub), size = 1.2 * k) +
    geom_smooth(method = "lm", formula = y ~ x, colour = "black", linewidth = 0.6) +
    geom_vline(xintercept = cfg$mm_cut, linetype = "dotted", colour = "red") +
    geom_hline(yintercept = cfg$gs_cut, linetype = "dotted", colour = "red") +
    scale_colour_manual(values = c(`TRUE` = "red", `FALSE` = m), guide = "none") +
    coord_cartesian(xlim = c(min(d$consensusMM), 1), ylim = c(0, 1)) +
    labs(title = paste(cap(m), "module"), x = "Module membership (MM)", y = "Gene significance (GS)\nin module direction") +
    theme_figure(k)
}
