################################################################################
# 03_biological_pathway_visualization.R
#
# Biological pathways represented by mortality-associated protein signatures
#
# Manuscript section:
#   Biological pathways represented by mortality-associated protein signatures
#
# Purpose:
#   1. Read sex-specific LASSO-selected proteins from Part 2.
#   2. Generate Venn diagrams for male/female overlap:
#        Figure 1C: Olink Target Inflammation panel
#        Figure 1D: Olink Explore proteomic panel
#   3. Read KEGG pathway enrichment results generated from DAVID.
#   4. Generate bubble + Sankey plots:
#        Figure 1G: KEGG enrichment for the union of male/female proteins from
#                   the Olink Target Inflammation panel
#        Figure 1H: KEGG enrichment for the union of male/female proteins from
#                   the Olink Explore proteomic panel
#
# Required inputs from Part 2:
#   results/selected_proteins/selected_proteins_all_panels.csv
#
# Required DAVID KEGG enrichment input files:
#   data/kegg_target_inflammation_union.csv
#   data/kegg_explore_union.csv
#
# Expected KEGG input columns:
#   Description   KEGG pathway name, preferably including hsa ID
#   Gene_Names    contributing proteins separated by comma, semicolon, slash, or pipe
#   GeneRatio     numeric gene ratio or percentage; if absent, calculated from Count / universe
#   pvalue        enrichment P value
#   Class         optional KEGG class for Sankey pathway coloring
#
# Example Description:
#   hsa04060:Cytokine-cytokine receptor interaction
#
# Outputs:
#   results/figures/
#     Figure_1C_venn_target_inflammation.pdf/png/svg
#     Figure_1D_venn_explore.pdf/png/svg
#     Figure_1G_kegg_target_inflammation_bubble_sankey.pdf/png/svg
#     Figure_1H_kegg_explore_bubble_sankey.pdf/png/svg
#
#   results/tables/
#     Figure_1C_venn_target_inflammation_overlap.csv
#     Figure_1D_venn_explore_overlap.csv
#     Figure_1G_kegg_target_inflammation_processed.csv
#     Figure_1H_kegg_explore_processed.csv
################################################################################


# ==============================================================================
# 0. Packages
# ==============================================================================

required_packages <- c(
  "dplyr", "tidyr", "purrr", "stringr", "readr", "tibble",
  "ggplot2", "ggrepel", "ggsankey", "VennDiagram",
  "grid", "svglite", "patchwork"
)

install_if_missing <- function(pkgs) {
  missing_pkgs <- pkgs[!pkgs %in% rownames(installed.packages())]
  if (length(missing_pkgs) > 0) install.packages(missing_pkgs)
  invisible(lapply(pkgs, library, character.only = TRUE))
}

install_if_missing(required_packages)


# ==============================================================================
# 1. User configuration
# ==============================================================================

set.seed(20240618)

# ---- Inputs ----
file_selected_proteins <- "results/selected_proteins/selected_proteins_all_panels.csv"

file_kegg_target_inflammation <- "data/kegg_target_inflammation_union.csv"
file_kegg_explore             <- "data/kegg_explore_union.csv"

# ---- Output directories ----
dir_results <- "results"
dir_tables  <- file.path(dir_results, "tables")
dir_figures <- file.path(dir_results, "figures")

dir.create(dir_tables,  recursive = TRUE, showWarnings = FALSE)
dir.create(dir_figures, recursive = TRUE, showWarnings = FALSE)

# ---- Panel naming from Part 2 ----
panel_target_inflammation <- "target_inflammation"
panel_explore             <- "explore"

# ---- Display labels ----
panel_labels <- c(
  target_inflammation = "Olink Target Inflammation panel",
  explore = "Olink Explore proteomic panel"
)

sex_order <- c("Men", "Women")

# ---- Figure aesthetics ----
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

bubble_colors <- grDevices::colorRampPalette(c("#FEE08B", "#F46D43", "#A50026"))(100)


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
    dplyr::filter(panel %in% c(panel_target_inflammation, panel_explore),
                  sex %in% sex_order,
                  !is.na(protein),
                  protein != "") %>%
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

  png(file.path(dir_figures, paste0(file_prefix, ".png")),
      width = width, height = height, units = "in", res = 600)
  grid::grid.newpage()
  grid::grid.draw(venn_grob)
  dev.off()

  svglite::svglite(file.path(dir_figures, paste0(file_prefix, ".svg")),
                   width = width, height = height)
  grid::grid.newpage()
  grid::grid.draw(venn_grob)
  dev.off()
}

clean_kegg_description <- function(x) {
  stringr::str_replace(x, "^hsa[0-9]+:", "")
}

extract_kegg_id <- function(x) {
  out <- stringr::str_extract(x, "hsa[0-9]+")
  ifelse(is.na(out), x, out)
}

read_kegg_results <- function(path, panel_name) {
  require_file(path, paste0("KEGG enrichment file for ", panel_name))

  kegg <- readr::read_csv(path, show_col_types = FALSE)

  # Flexible column handling for DAVID exports.
  if (!"Description" %in% names(kegg)) {
    possible <- intersect(c("Term", "term", "Pathway", "pathway", "Description"), names(kegg))
    if (length(possible) == 0) stop("KEGG file must contain Description or Term/Pathway column.")
    kegg <- kegg %>% dplyr::rename(Description = dplyr::all_of(possible[1]))
  }

  if (!"Gene_Names" %in% names(kegg)) {
    possible <- intersect(c("Genes", "genes", "Gene", "gene", "geneID", "gene_id"), names(kegg))
    if (length(possible) == 0) stop("KEGG file must contain Gene_Names or Genes/geneID column.")
    kegg <- kegg %>% dplyr::rename(Gene_Names = dplyr::all_of(possible[1]))
  }

  if (!"pvalue" %in% names(kegg)) {
    possible <- intersect(c("PValue", "P.Value", "P_value", "p_value", "P"), names(kegg))
    if (length(possible) == 0) stop("KEGG file must contain pvalue or PValue column.")
    kegg <- kegg %>% dplyr::rename(pvalue = dplyr::all_of(possible[1]))
  }

  if (!"Class" %in% names(kegg)) {
    kegg$Class <- "Unclassified"
  }

  if (!"GeneRatio" %in% names(kegg)) {
    if ("Count" %in% names(kegg)) {
      kegg$GeneRatio <- as.numeric(kegg$Count)
    } else {
      kegg$GeneRatio <- stringr::str_count(as.character(kegg$Gene_Names), "[,;/|]") + 1
    }
  }

  kegg %>%
    dplyr::mutate(
      panel = panel_name,
      Description = as.character(Description),
      Pathway = clean_kegg_description(Description),
      KEGG_ID = extract_kegg_id(Description),
      Gene_Names = as.character(Gene_Names),
      GeneRatio = as.numeric(GeneRatio),
      pvalue = as.numeric(pvalue),
      Class = dplyr::if_else(is.na(Class) | Class == "", "Unclassified", as.character(Class)),
      logP = -log10(pvalue)
    ) %>%
    dplyr::filter(!is.na(pvalue), !is.na(GeneRatio), pvalue > 0) %>%
    dplyr::arrange(pvalue, dplyr::desc(GeneRatio))
}

prepare_kegg_for_sankey <- function(kegg) {
  kegg_long <- kegg %>%
    dplyr::select(Description, Pathway, KEGG_ID, Gene_Names, Class, GeneRatio, pvalue, logP) %>%
    tidyr::separate_rows(Gene_Names, sep = "\\s*[,;/|]\\s*") %>%
    dplyr::mutate(
      Gene_Names = stringr::str_trim(Gene_Names),
      Pathway_label = paste0(KEGG_ID, ": ", Pathway)
    ) %>%
    dplyr::filter(Gene_Names != "")

  pathway_order <- kegg %>%
    dplyr::arrange(logP) %>%
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
  kegg_plot <- kegg %>%
    dplyr::mutate(
      KEGG_ID = factor(KEGG_ID, levels = KEGG_ID[order(logP, decreasing = TRUE)])
    )

  ggplot2::ggplot(kegg_plot, ggplot2::aes(x = GeneRatio, y = logP)) +
    ggplot2::geom_point(ggplot2::aes(size = logP, color = logP), shape = 16, alpha = 0.90) +
    ggplot2::scale_size_continuous(range = c(2.5, 8)) +
    ggplot2::scale_color_gradientn(colors = bubble_colors, name = expression(-log[10](P))) +
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
      y = expression(-log[10](P)),
      title = title
    ) +
    ggplot2::guides(size = "none") +
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
    dplyr::mutate(color = dplyr::if_else(is.na(color), kegg_class_colors["Unclassified"], color))

  gene_colors <- tibble::tibble(
    Gene_Names = levels(sankey$Gene_Names),
    color = "#BDBDBD"
  )

  sankey_long <- sankey %>%
    ggsankey::make_long(Pathway_label, Gene_Names) %>%
    dplyr::mutate(
      node_chr = as.character(node),
      color = dplyr::case_when(
        node_chr %in% as.character(gene_colors$Gene_Names) ~ gene_colors$color[match(node_chr, as.character(gene_colors$Gene_Names))],
        node_chr %in% as.character(pathway_colors$Pathway_label) ~ pathway_colors$color[match(node_chr, as.character(pathway_colors$Pathway_label))],
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

run_kegg_visualization <- function(kegg_file, panel_name, figure_prefix) {
  kegg <- read_kegg_results(kegg_file, panel_name)

  processed_file <- file.path(dir_tables, paste0(figure_prefix, "_processed.csv"))
  readr::write_csv(kegg, processed_file)

  bubble <- plot_kegg_bubble(kegg)
  sankey <- plot_kegg_sankey(kegg)
  combined <- combine_bubble_sankey(bubble, sankey)

  save_plot_all_formats(
    combined,
    file_prefix = paste0(figure_prefix, "_bubble_sankey"),
    width = 10.5,
    height = 6.0
  )

  invisible(list(kegg = kegg, bubble = bubble, sankey = sankey, combined = combined))
}


# ==============================================================================
# 3. Read selected proteins and generate Venn diagrams
# ==============================================================================

selected_proteins <- read_selected_proteins(file_selected_proteins)

target_sets <- get_panel_protein_sets(selected_proteins, panel_target_inflammation)
explore_sets <- get_panel_protein_sets(selected_proteins, panel_explore)

# ---- Figure 1C: Target Inflammation Venn ----
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

# ---- Figure 1D: Explore Venn ----
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

# ---- Union protein lists for DAVID input and reproducibility ----
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
  file.path(dir_tables, "Figure_1G_DAVID_input_target_inflammation_union_proteins.csv")
)

readr::write_csv(
  union_explore,
  file.path(dir_tables, "Figure_1H_DAVID_input_explore_union_proteins.csv")
)


# ==============================================================================
# 4. Generate KEGG bubble and Sankey plots from DAVID output
# ==============================================================================

# ---- Figure 1G: Target Inflammation KEGG enrichment ----
kegg_target <- run_kegg_visualization(
  kegg_file = file_kegg_target_inflammation,
  panel_name = panel_target_inflammation,
  figure_prefix = "Figure_1G_kegg_target_inflammation"
)

# ---- Figure 1H: Explore KEGG enrichment ----
kegg_explore <- run_kegg_visualization(
  kegg_file = file_kegg_explore,
  panel_name = panel_explore,
  figure_prefix = "Figure_1H_kegg_explore"
)


# ==============================================================================
# 5. Console summary
# ==============================================================================

cat("\nBiological pathway visualization completed successfully.\n")
cat("\nFigure 1C/D overlap summary:\n")
print(dplyr::bind_rows(target_overlap, explore_overlap))

cat("\nGenerated outputs:\n")
cat("  - Venn diagrams: ", dir_figures, "\n", sep = "")
cat("  - KEGG bubble + Sankey plots: ", dir_figures, "\n", sep = "")
cat("  - Processed overlap and KEGG tables: ", dir_tables, "\n", sep = "")
