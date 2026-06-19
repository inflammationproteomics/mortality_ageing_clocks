################################################################################
# 03_biological_pathway_visualization.R
#
# Biological pathways represented by mortality-associated protein signatures
#
# Purpose:
#   1. Read sex-specific LASSO-selected proteins from Part 2.
#   2. Generate Venn diagrams for male/female overlap:
#        Figure 1C: Olink Target Inflammation panel
#        Figure 1D: Olink Explore proteomic panel
#   3. Perform KEGG pathway enrichment using clusterProfiler on the union of
#      proteins selected in men and women within each panel.
#   4. Generate bubble + Sankey plots:
#        Figure 1G: KEGG enrichment for the union of male/female proteins from
#                   the Olink Target Inflammation panel
#        Figure 1H: KEGG enrichment for the union of male/female proteins from
#                   the Olink Explore proteomic panel
################################################################################


# ==============================================================================
# 0. Packages
# ==============================================================================

cran_packages <- c(
  "dplyr", "tidyr", "purrr", "stringr", "readr", "tibble",
  "ggplot2", "ggrepel", "ggsankey", "VennDiagram",
  "grid", "svglite", "patchwork"
)

bioc_packages <- c(
  "clusterProfiler", "org.Hs.eg.db"
)

install_cran_if_missing <- function(pkgs) {
  missing_pkgs <- pkgs[!pkgs %in% rownames(installed.packages())]
  if (length(missing_pkgs) > 0) install.packages(missing_pkgs)
  invisible(lapply(pkgs, library, character.only = TRUE))
}

install_bioc_if_missing <- function(pkgs) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  missing_pkgs <- pkgs[!pkgs %in% rownames(installed.packages())]
  if (length(missing_pkgs) > 0) {
    BiocManager::install(missing_pkgs, update = FALSE, ask = FALSE)
  }
  invisible(lapply(pkgs, library, character.only = TRUE))
}

install_cran_if_missing(cran_packages)
install_bioc_if_missing(bioc_packages)


# ==============================================================================
# 1. User configuration
# ==============================================================================

set.seed(20240618)

file_selected_proteins <- "results/selected_proteins/selected_proteins_all_panels.csv"

dir_results <- "results"
dir_tables  <- file.path(dir_results, "tables")
dir_figures <- file.path(dir_results, "figures")

dir.create(dir_tables,  recursive = TRUE, showWarnings = FALSE)
dir.create(dir_figures, recursive = TRUE, showWarnings = FALSE)

panel_target_inflammation <- "target_inflammation"
panel_explore             <- "explore"

panel_labels <- c(
  target_inflammation = "Olink Target Inflammation panel",
  explore = "Olink Explore proteomic panel"
)

sex_order <- c("Men", "Women")

base_family <- "sans"

venn_colors <- list(
  target_inflammation = c("#00BFC4", "#C77CFF"),
  explore = c("#F8766D", "#7CAE00")
)

kegg_class_colors <- c(
  "Environmental Information Processing" = "#E88E8F",
  "Organismal Systems" = "#FFA61D",
  "Cellular Processes" = "#9A99E1",
  "Human Diseases" = "#0099B4",
  "Metabolism" = "#56B4E9",
  "Genetic Information Processing" = "#66C2A5",
  "Unclassified" = "#BDBDBD"
)

bubble_colors <- grDevices::colorRampPalette(
  c("#FEE08B", "#F46D43", "#A50026")
)(100)

fdr_cutoff <- 0.05
organism_kegg <- "hsa"


# ==============================================================================
# 2. Helper functions
# ==============================================================================

require_file <- function(path, description) {
  if (!file.exists(path)) {
    stop(
      description, " not found: ", path, "\n",
      "Please create this file before running the script."
    )
  }
}

normalize_panel_name <- function(x) {
  x <- as.character(x)
  dplyr::case_when(
    x %in% c("target_inflammation", "Target Inflammation", "inflammation", "Inflammation") ~ "target_inflammation",
    x %in% c("explore", "Explore", "proteomic", "proteomics", "Olink Explore") ~ "explore",
    TRUE ~ x
  )
}

normalize_sex_name <- function(x) {
  x <- as.character(x)
  dplyr::case_when(
    x %in% c("Men", "Male", "M", "1") ~ "Men",
    x %in% c("Women", "Female", "F", "0", "2") ~ "Women",
    TRUE ~ x
  )
}

read_selected_proteins <- function(path) {
  require_file(path, "Selected protein file from Part 2")
  
  selected <- readr::read_csv(path, show_col_types = FALSE)
  
  required_cols <- c("panel", "sex", "protein")
  missing_cols <- setdiff(required_cols, names(selected))
  if (length(missing_cols) > 0) {
    stop("Selected protein file is missing columns: ", paste(missing_cols, collapse = ", "))
  }
  
  selected %>%
    dplyr::mutate(
      panel = normalize_panel_name(panel),
      sex = normalize_sex_name(sex),
      protein = as.character(protein)
    ) %>%
    dplyr::filter(
      panel %in% c(panel_target_inflammation, panel_explore),
      sex %in% sex_order,
      !is.na(protein),
      protein != ""
    ) %>%
    dplyr::distinct(panel, sex, protein)
}

get_panel_protein_sets <- function(selected, panel_name) {
  panel_dat <- selected %>% dplyr::filter(panel == panel_name)
  
  men <- panel_dat %>%
    dplyr::filter(sex == "Men") %>%
    dplyr::pull(protein) %>%
    unique() %>%
    sort()
  
  women <- panel_dat %>%
    dplyr::filter(sex == "Women") %>%
    dplyr::pull(protein) %>%
    unique() %>%
    sort()
  
  if (length(men) == 0 || length(women) == 0) {
    stop("No selected proteins found for both sexes in panel: ", panel_name)
  }
  
  list(Men = men, Women = women)
}

make_venn_overlap_table <- function(protein_sets, panel_name) {
  men <- protein_sets$Men
  women <- protein_sets$Women
  
  tibble::tibble(
    panel = panel_name,
    category = c("Men only", "Women only", "Overlap", "Union"),
    n = c(
      length(setdiff(men, women)),
      length(setdiff(women, men)),
      length(intersect(men, women)),
      length(union(men, women))
    ),
    proteins = c(
      paste(setdiff(men, women), collapse = "; "),
      paste(setdiff(women, men), collapse = "; "),
      paste(intersect(men, women), collapse = "; "),
      paste(union(men, women), collapse = "; ")
    )
  )
}

plot_venn_diagram <- function(protein_sets, panel_name) {
  VennDiagram::venn.diagram(
    x = list(
      Men = protein_sets$Men,
      Women = protein_sets$Women
    ),
    filename = NULL,
    fill = venn_colors[[panel_name]],
    alpha = 0.45,
    lwd = 0,
    cex = 1.4,
    cat.cex = 1.2,
    cat.fontface = "bold",
    cat.pos = c(-25, 25),
    cat.dist = 0.05,
    margin = 0.08
  )
}

save_venn_outputs <- function(venn_grob, file_prefix, width = 5.2, height = 4.6) {
  pdf(file.path(dir_figures, paste0(file_prefix, ".pdf")), width = width, height = height)
  grid::grid.newpage()
  grid::grid.draw(venn_grob)
  dev.off()
  
  png(
    file.path(dir_figures, paste0(file_prefix, ".png")),
    width = width,
    height = height,
    units = "in",
    res = 600
  )
  grid::grid.newpage()
  grid::grid.draw(venn_grob)
  dev.off()
  
  svglite::svglite(
    file.path(dir_figures, paste0(file_prefix, ".svg")),
    width = width,
    height = height
  )
  grid::grid.newpage()
  grid::grid.draw(venn_grob)
  dev.off()
}

map_symbols_to_entrez <- function(symbols, panel_name) {
  mapped <- clusterProfiler::bitr(
    symbols,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db::org.Hs.eg.db
  ) %>%
    dplyr::distinct(SYMBOL, ENTREZID)
  
  unmapped <- setdiff(symbols, mapped$SYMBOL)
  
  readr::write_csv(
    tibble::tibble(panel = panel_name, unmapped_symbol = unmapped),
    file.path(dir_tables, paste0("unmapped_symbols_", panel_name, ".csv"))
  )
  
  if (nrow(mapped) == 0) {
    stop("No proteins could be mapped to ENTREZ IDs for panel: ", panel_name)
  }
  
  mapped
}

run_kegg_enrichment <- function(symbols, panel_name) {
  mapped <- map_symbols_to_entrez(symbols, panel_name)
  
  ekegg <- clusterProfiler::enrichKEGG(
    gene = unique(mapped$ENTREZID),
    organism = organism_kegg,
    keyType = "kegg",
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    qvalueCutoff = 1
  )
  
  if (is.null(ekegg) || nrow(as.data.frame(ekegg)) == 0) {
    warning("No KEGG enrichment results found for panel: ", panel_name)
    return(
      list(
        mapped = mapped,
        enrichment = tibble::tibble()
      )
    )
  }
  
  kegg_df <- as.data.frame(ekegg) %>%
    tibble::as_tibble() %>%
    dplyr::mutate(
      panel = panel_name,
      KEGG_ID = ID,
      Pathway = Description,
      GeneRatio_numeric = purrr::map_dbl(GeneRatio, function(x) {
        parts <- strsplit(x, "/", fixed = TRUE)[[1]]
        as.numeric(parts[1]) / as.numeric(parts[2])
      }),
      Count = as.numeric(Count),
      logFDR = -log10(p.adjust),
      geneID = as.character(geneID)
    )
  
  symbol_lookup <- mapped %>%
    dplyr::distinct(ENTREZID, SYMBOL)
  
  kegg_df <- kegg_df %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      Gene_Names = paste(
        symbol_lookup$SYMBOL[
          match(
            strsplit(geneID, "/", fixed = TRUE)[[1]],
            symbol_lookup$ENTREZID
          )
        ],
        collapse = "; "
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      Gene_Names = stringr::str_replace_all(Gene_Names, "NA;\\s*|;\\s*NA|NA", "")
    ) %>%
    dplyr::arrange(p.adjust, dplyr::desc(Count))
  
  list(
    mapped = mapped,
    enrichment = kegg_df
  )
}

get_kegg_class <- function(kegg_ids) {
  # Optional KEGG class annotation.
  # To avoid unstable internet-dependent KEGGREST calls, pathways are left as
  # Unclassified by default. If you have a curated KEGG class table, merge it here.
  rep("Unclassified", length(kegg_ids))
}

prepare_kegg_table_for_plot <- function(kegg_df, panel_name) {
  if (nrow(kegg_df) == 0) return(kegg_df)
  
  kegg_df %>%
    dplyr::filter(!is.na(p.adjust), p.adjust < fdr_cutoff) %>%
    dplyr::mutate(
      Class = get_kegg_class(KEGG_ID),
      logFDR = -log10(p.adjust),
      GeneRatio = GeneRatio_numeric
    ) %>%
    dplyr::arrange(p.adjust, dplyr::desc(GeneRatio))
}

prepare_kegg_for_sankey <- function(kegg) {
  kegg_long <- kegg %>%
    dplyr::select(KEGG_ID, Pathway, Gene_Names, Class, GeneRatio, pvalue, p.adjust, logFDR) %>%
    tidyr::separate_rows(Gene_Names, sep = "\\s*[,;/|]\\s*") %>%
    dplyr::mutate(
      Gene_Names = stringr::str_trim(Gene_Names),
      Pathway_label = paste0(KEGG_ID, ": ", Pathway)
    ) %>%
    dplyr::filter(Gene_Names != "")
  
  pathway_order <- kegg %>%
    dplyr::arrange(logFDR) %>%
    dplyr::mutate(Pathway_label = paste0(KEGG_ID, ": ", Pathway)) %>%
    dplyr::pull(Pathway_label) %>%
    rev()
  
  gene_order <- kegg_long %>%
    dplyr::count(Gene_Names, sort = TRUE) %>%
    dplyr::pull(Gene_Names)
  
  list(
    long = kegg_long,
    pathway_order = pathway_order,
    gene_order = gene_order
  )
}

plot_kegg_bubble <- function(kegg, title = NULL) {
  if (nrow(kegg) == 0) {
    return(
      ggplot2::ggplot() +
        ggplot2::theme_void(base_family = base_family) +
        ggplot2::labs(title = paste0(title, ": no FDR-significant KEGG pathways"))
    )
  }
  
  kegg_plot <- kegg %>%
    dplyr::mutate(
      KEGG_ID = factor(KEGG_ID, levels = KEGG_ID[order(logFDR, decreasing = TRUE)])
    )
  
  ggplot2::ggplot(kegg_plot, ggplot2::aes(x = GeneRatio, y = logFDR)) +
    ggplot2::geom_point(ggplot2::aes(size = Count, color = logFDR), shape = 16, alpha = 0.90) +
    ggplot2::scale_size_continuous(range = c(2.5, 8), name = "Count") +
    ggplot2::scale_color_gradientn(
      colors = bubble_colors,
      name = expression(-log[10](FDR))
    ) +
    ggrepel::geom_text_repel(
      ggplot2::aes(label = KEGG_ID),
      size = 3.0,
      family = base_family,
      color = "black",
      box.padding = 0.45,
      max.overlaps = 50,
      min.segment.length = 0
    ) +
    ggplot2::labs(
      x = "Gene ratio",
      y = expression(-log[10](FDR)),
      title = title
    ) +
    ggplot2::theme_classic(base_family = base_family, base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 12),
      legend.position = c(0.02, 0.98),
      legend.justification = c(0, 1),
      legend.background = ggplot2::element_blank(),
      legend.key = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(size = 8),
      legend.text = ggplot2::element_text(size = 8),
      axis.text = ggplot2::element_text(color = "black"),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.6)
    )
}

plot_kegg_sankey <- function(kegg) {
  if (nrow(kegg) == 0) {
    return(
      ggplot2::ggplot() +
        ggplot2::theme_void(base_family = base_family) +
        ggplot2::labs(title = "No FDR-significant KEGG pathways")
    )
  }
  
  prepared <- prepare_kegg_for_sankey(kegg)
  sankey <- prepared$long
  
  sankey <- sankey %>%
    dplyr::mutate(
      Pathway_label = factor(Pathway_label, levels = prepared$pathway_order),
      Gene_Names = factor(Gene_Names, levels = prepared$gene_order)
    )
  
  pathway_colors <- sankey %>%
    dplyr::distinct(Pathway_label, Class) %>%
    dplyr::mutate(color = kegg_class_colors[Class]) %>%
    dplyr::mutate(
      color = dplyr::if_else(
        is.na(color),
        kegg_class_colors["Unclassified"],
        color
      )
    )
  
  gene_colors <- tibble::tibble(
    Gene_Names = levels(sankey$Gene_Names),
    color = "#BDBDBD"
  )
  
  sankey_long <- sankey %>%
    ggsankey::make_long(Pathway_label, Gene_Names) %>%
    dplyr::mutate(
      node_chr = as.character(node),
      color = dplyr::case_when(
        node_chr %in% as.character(gene_colors$Gene_Names) ~
          gene_colors$color[match(node_chr, as.character(gene_colors$Gene_Names))],
        node_chr %in% as.character(pathway_colors$Pathway_label) ~
          pathway_colors$color[match(node_chr, as.character(pathway_colors$Pathway_label))],
        TRUE ~ "#BDBDBD"
      )
    )
  
  ggplot2::ggplot(
    sankey_long,
    ggplot2::aes(
      x = x,
      next_x = next_x,
      node = node,
      next_node = next_node,
      fill = color,
      label = node
    )
  ) +
    ggplot2::scale_fill_identity() +
    ggsankey::geom_sankey(flow.alpha = 0.42, smooth = 8, width = 0.08) +
    ggsankey::geom_sankey_text(
      size = 2.7,
      color = "black",
      family = base_family,
      hjust = 0,
      nudge_x = 0.015
    ) +
    ggplot2::theme_void(base_family = base_family) +
    ggplot2::theme(
      legend.position = "none",
      plot.margin = ggplot2::margin(2, 18, 2, 2, "pt")
    ) +
    ggplot2::scale_x_discrete(expand = ggplot2::expansion(mult = c(0.10, 0.18)))
}

combine_bubble_sankey <- function(bubble_plot, sankey_plot) {
  bubble_plot + sankey_plot +
    patchwork::plot_layout(widths = c(0.38, 0.62))
}

save_plot_all_formats <- function(plot, file_prefix, width = 10.5, height = 6.0) {
  ggplot2::ggsave(
    filename = file.path(dir_figures, paste0(file_prefix, ".pdf")),
    plot = plot,
    width = width,
    height = height,
    device = cairo_pdf
  )
  
  ggplot2::ggsave(
    filename = file.path(dir_figures, paste0(file_prefix, ".png")),
    plot = plot,
    width = width,
    height = height,
    dpi = 600
  )
  
  svglite::svglite(
    filename = file.path(dir_figures, paste0(file_prefix, ".svg")),
    width = width,
    height = height,
    system_fonts = list(sans = base_family)
  )
  print(plot)
  dev.off()
}

run_clusterprofiler_kegg_visualization <- function(symbols, panel_name, figure_prefix) {
  enrichment_result <- run_kegg_enrichment(symbols, panel_name)
  
  mapped_file <- file.path(dir_tables, paste0(figure_prefix, "_symbol_to_entrez_mapping.csv"))
  readr::write_csv(enrichment_result$mapped, mapped_file)
  
  full_file <- file.path(dir_tables, paste0(figure_prefix, "_kegg_clusterProfiler_full.csv"))
  readr::write_csv(enrichment_result$enrichment, full_file)
  
  kegg_sig <- prepare_kegg_table_for_plot(enrichment_result$enrichment, panel_name)
  
  sig_file <- file.path(dir_tables, paste0(figure_prefix, "_kegg_clusterProfiler_FDR_significant.csv"))
  readr::write_csv(kegg_sig, sig_file)
  
  bubble <- plot_kegg_bubble(kegg_sig)
  sankey <- plot_kegg_sankey(kegg_sig)
  combined <- combine_bubble_sankey(bubble, sankey)
  
  save_plot_all_formats(
    combined,
    file_prefix = paste0(figure_prefix, "_kegg_bubble_sankey"),
    width = 10.5,
    height = 6.0
  )
  
  invisible(
    list(
      mapped = enrichment_result$mapped,
      full = enrichment_result$enrichment,
      significant = kegg_sig,
      bubble = bubble,
      sankey = sankey,
      combined = combined
    )
  )
}


# ==============================================================================
# 3. Read selected proteins and generate Venn diagrams
# ==============================================================================

selected_proteins <- read_selected_proteins(file_selected_proteins)

target_sets <- get_panel_protein_sets(selected_proteins, panel_target_inflammation)
explore_sets <- get_panel_protein_sets(selected_proteins, panel_explore)

target_overlap <- make_venn_overlap_table(target_sets, panel_target_inflammation)
readr::write_csv(
  target_overlap,
  file.path(dir_tables, "Figure_1C_venn_target_inflammation_overlap.csv")
)

venn_target <- plot_venn_diagram(target_sets, panel_target_inflammation)
save_venn_outputs(
  venn_target,
  file_prefix = "Figure_1C_venn_target_inflammation"
)

explore_overlap <- make_venn_overlap_table(explore_sets, panel_explore)
readr::write_csv(
  explore_overlap,
  file.path(dir_tables, "Figure_1D_venn_explore_overlap.csv")
)

venn_explore <- plot_venn_diagram(explore_sets, panel_explore)
save_venn_outputs(
  venn_explore,
  file_prefix = "Figure_1D_venn_explore"
)


# ==============================================================================
# 4. Create union protein lists for KEGG enrichment
# ==============================================================================

union_target <- tibble::tibble(
  panel = panel_target_inflammation,
  protein = sort(union(target_sets$Men, target_sets$Women))
)

union_explore <- tibble::tibble(
  panel = panel_explore,
  protein = sort(union(explore_sets$Men, explore_sets$Women))
)

readr::write_csv(
  union_target,
  file.path(dir_tables, "Figure_1G_clusterProfiler_input_target_inflammation_union_proteins.csv")
)

readr::write_csv(
  union_explore,
  file.path(dir_tables, "Figure_1H_clusterProfiler_input_explore_union_proteins.csv")
)


# ==============================================================================
# 5. Run clusterProfiler KEGG enrichment and generate plots
# ==============================================================================

kegg_target <- run_clusterprofiler_kegg_visualization(
  symbols = union_target$protein,
  panel_name = panel_target_inflammation,
  figure_prefix = "Figure_1G_target_inflammation"
)

kegg_explore <- run_clusterprofiler_kegg_visualization(
  symbols = union_explore$protein,
  panel_name = panel_explore,
  figure_prefix = "Figure_1H_explore"
)


# ==============================================================================
# 6. Console summary
# ==============================================================================

cat("\nBiological pathway visualization completed successfully.\n")

cat("\nFigure 1C/D overlap summary:\n")
print(dplyr::bind_rows(target_overlap, explore_overlap))

cat("\nKEGG enrichment summary:\n")
cat("  Target Inflammation panel:\n")
cat("    Input proteins: ", nrow(union_target), "\n", sep = "")
cat("    Mapped proteins: ", nrow(kegg_target$mapped), "\n", sep = "")
cat("    FDR-significant KEGG pathways: ", nrow(kegg_target$significant), "\n", sep = "")

cat("  Olink Explore proteomic panel:\n")
cat("    Input proteins: ", nrow(union_explore), "\n", sep = "")
cat("    Mapped proteins: ", nrow(kegg_explore$mapped), "\n", sep = "")
cat("    FDR-significant KEGG pathways: ", nrow(kegg_explore$significant), "\n", sep = "")

cat("\nGenerated outputs:\n")
cat("  - Venn diagrams: ", dir_figures, "\n", sep = "")
cat("  - KEGG bubble + Sankey plots: ", dir_figures, "\n", sep = "")
cat("  - Processed overlap, mapping, and KEGG tables: ", dir_tables, "\n", sep = "")

cat("\nPackage versions:\n")
cat("  R: ", R.version.string, "\n", sep = "")
cat("  clusterProfiler: ", as.character(utils::packageVersion("clusterProfiler")), "\n", sep = "")
cat("  org.Hs.eg.db: ", as.character(utils::packageVersion("org.Hs.eg.db")), "\n", sep = "")
