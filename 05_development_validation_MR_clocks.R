################################################################################
# 05_01_build_mr_clocks.R
################################################################################

library(dplyr)
library(readr)
library(glmnet)

dir_data <- "results/analysis_datasets"
dir_tables <- "results/tables"
dir_models <- "results/models"

ukb_derivation <- readRDS(file.path(dir_data, "ukb_derivation.rds"))
ukb_internal <- readRDS(file.path(dir_data, "ukb_internal_validation.rds"))
esther_TI <- readRDS(file.path(dir_data, "esther_target_inflammation.rds"))
esther_HT <- readRDS(file.path(dir_data, "esther_explore_ht.rds"))

coef_tbl <- read_csv(file.path(dir_tables,
                               "Supplementary_Table_10_clock_coefficients.csv"))

predict_raw_age <- function(df, coef_df){

  intercept <- coef_df %>%
    filter(Variable == "Intercept") %>%
    pull(Coefficient)

  beta_df <- coef_df %>%
    filter(Variable != "Intercept")

  x <- as.matrix(df[, beta_df$Variable])

  intercept + x %*% beta_df$Coefficient

}

calibrate_age <- function(raw_age, coef_df){

  a <- coef_df %>%
    filter(Variable == "Calibration_Intercept") %>%
    pull(Coefficient)

  b <- coef_df %>%
    filter(Variable == "Calibration_Slope") %>%
    pull(Coefficient)

  a + b * raw_age

}

################################################################################
# 05_02_generate_figure2.R
#
# Figure 2
#
# Distributions of chronological age and mortality risk ageing clocks
# in the UK Biobank and ESTHER cohorts
#
# Corresponding manuscript components:
#   Figure 2A–F
#
# Inputs:
#   results/tables/mr_ages_derivation.csv
#   results/tables/mr_ages_internal.csv
#   results/tables/mr_ages_esther_TI.csv
#   results/tables/mr_ages_esther_HT.csv
#
# Outputs:
#   results/figures/Figure2.pdf
#   results/figures/Figure2.png
#
################################################################################

library(dplyr)
library(readr)
library(ggplot2)
library(patchwork)

dir_tables <- "results/tables"
dir_figures <- "results/figures"

ukb_derivation <- read_csv(file.path(dir_tables, "mr_ages_derivation.csv"))
ukb_internal <- read_csv(file.path(dir_tables, "mr_ages_internal.csv"))
esther_TI <- read_csv(file.path(dir_tables, "mr_ages_esther_TI.csv"))
esther_HT <- read_csv(file.path(dir_tables, "mr_ages_esther_HT.csv"))

plot_histogram <- function(df, clock, title){

  ggplot(df) +
    geom_histogram(aes(age_share),
                   fill = "#3B7DDD",
                   alpha = 0.5,
                   bins = 30) +
    geom_histogram(aes(.data[[clock]]),
                   fill = "#D84B4B",
                   alpha = 0.5,
                   bins = 30) +
    labs(title = title,
         x = "Age (years)",
         y = "Participants") +
    theme_bw()

}

################################################################################
# 05_03_generate_figure3.R
#
# Figure 3
#
# Comparison of the C-statistics of chronological age and the
# Inflammation-MR-Clocks and Proteomics-MR-Clocks for prediction
# of 10-year all-cause mortality
#
# Corresponding manuscript components:
#   Figure 3
#
# Inputs:
#   results/tables/Figure3_Cstatistics.csv
#
# Outputs:
#   results/figures/Figure3.pdf
#   results/figures/Figure3.png
#   results/figures/Figure3.pptx
#
################################################################################

library(dplyr)
library(readr)
library(ggplot2)
library(officer)
library(rvg)

dir_tables <- "results/tables"
dir_figures <- "results/figures"

################################################################################
# Load C-statistics
################################################################################

df_cstat <- read_csv(file.path(dir_tables,
                               "Figure3_Cstatistics.csv"))

################################################################################
# Format data
################################################################################

df_cstat <- df_cstat %>%
  mutate(
    Dataset = factor(
      Dataset,
      levels = c(
        "UKB 70%",
        "UKB 30%",
        "ESTHER"
      )
    ),
    Group = factor(
      Group,
      levels = c(
        "Men",
        "Women",
        "Overall"
      )
    ),
    Age_type = factor(
      Age_type,
      levels = c(
        "Chronological age",
        "Inflammation-MR-Clock1",
        "Inflammation-MR-Clock2",
        "Proteomics-MR-Clock1",
        "Proteomics-MR-Clock2"
      )
    ),
    label = sprintf(
      "%.3f\n(%.3f–%.3f)",
      Cstat,
      Cstat_lower,
      Cstat_upper
    )
  )

################################################################################
# Generate heatmap
################################################################################

p_heatmap <- ggplot(
  df_cstat,
  aes(
    x = Group,
    y = Age_type,
    fill = Cstat
  )
) +
  geom_tile(
    color = "white"
  ) +
  geom_text(
    aes(label = label),
    size = 3.4,
    color = "black",
    lineheight = 0.95
  ) +
  scale_fill_gradientn(
    colours = c(
      "#f7fbff",
      "#6baed6",
      "#08306b"
    ),
    limits = c(0.65, 0.80),
    oob = scales::squish,
    name = "C-statistic"
  ) +
  facet_wrap(
    ~Dataset,
    ncol = 3
  ) +
  scale_y_discrete(
    limits = rev(
      levels(df_cstat$Age_type)
    )
  ) +
  theme_minimal(
    base_size = 13
  ) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(
      fill = "grey90",
      color = NA
    ),
    strip.text = element_text(
      size = 13,
      face = "bold"
    ),
    axis.title = element_blank(),
    axis.text.x = element_text(
      size = 11
    ),
    axis.text.y = element_text(
      size = 11
    ),
    legend.position = "right"
  )

################################################################################
# Export Figure 3
################################################################################

ggsave(file.path(dir_figures,
                 "Figure3.pdf"),
       p_heatmap,
       width = 11,
       height = 5)

ggsave(file.path(dir_figures,
                 "Figure3.png"),
       p_heatmap,
       dpi = 600,
       width = 11,
       height = 5)

doc <- read_pptx()
doc <- add_slide(doc,
                 layout = "Blank",
                 master = "Office Theme")

doc <- ph_with(
  doc,
  dml(ggobj = p_heatmap),
  location = ph_location(
    left = 0.3,
    top = 0.3,
    width = 12.8,
    height = 6.8
  )
)

print(doc,
      target = file.path(dir_figures,
                         "Figure3.pptx"))

################################################################################
# 05_04_generate_supplementary_figures6_7.R
################################################################################

library(dplyr)
library(readr)
library(ggplot2)
library(patchwork)

dir_tables <- "results/tables"
dir_figures <- "results/figures"

ukb_derivation <- read_csv(file.path(dir_tables, "mr_ages_derivation.csv"))
ukb_internal <- read_csv(file.path(dir_tables, "mr_ages_internal.csv"))
esther_TI <- read_csv(file.path(dir_tables, "mr_ages_esther_TI.csv"))
esther_HT <- read_csv(file.path(dir_tables, "mr_ages_esther_HT.csv"))

plot_cor <- function(df, clock, sex_name, title){

  if(sex_name != "Overall"){
    df <- filter(df, sex_share == sex_name)
  }

  r <- cor(df$age_share,
           df[[clock]],
           method = "pearson",
           use = "complete.obs")

  mae <- mean(abs(df$age_share - df[[clock]]),
              na.rm = TRUE)

  ggplot(df,
         aes(age_share,
             .data[[clock]])) +
    geom_point(alpha = 0.30,
               size = 0.6,
               color = "#3B7DDD") +
    geom_smooth(method = "lm",
                se = FALSE,
                color = "#D84B4B",
                linewidth = 0.7) +
    annotate("text",
             x = -Inf,
             y = Inf,
             hjust = -0.1,
             vjust = 1.2,
             label = paste0("r = ",
                            round(r, 3),
                            "\nMAE = ",
                            round(mae, 2)),
             size = 3.2) +
    labs(title = title,
         x = "Chronological age",
         y = "Estimated age") +
    theme_bw() +
    theme(
      plot.title = element_text(size = 10,
                                face = "bold",
                                hjust = 0.5),
      axis.title = element_text(size = 9),
      axis.text = element_text(size = 8),
      panel.grid = element_blank()
    )

}



################################################################################
# Supplementary Figure 6
################################################################################

plots_s6 <- list()

for(clock in c("Inflammation-MR-Clock1",
               "Inflammation-MR-Clock2")){

  for(sex in c("Men",
               "Women",
               "Overall")){

    plots_s6[[length(plots_s6)+1]] <-
      plot_cor(ukb_derivation,
               clock,
               sex,
               paste0(sex,"\nUKB 70%"))

    plots_s6[[length(plots_s6)+1]] <-
      plot_cor(ukb_internal,
               clock,
               sex,
               paste0(sex,"\nUKB 30%"))

    plots_s6[[length(plots_s6)+1]] <-
      plot_cor(esther_TI,
               clock,
               sex,
               paste0(sex,"\nESTHER"))

  }

}


FigureS6 <- wrap_plots(
  plots_s6,
  ncol = 3
)


ggsave(
  file.path(
    dir_figures,
    "Supplementary_Figure6.pdf"
  ),
  FigureS6,
  width = 12,
  height = 18
)


ggsave(
  file.path(
    dir_figures,
    "Supplementary_Figure6.png"
  ),
  FigureS6,
  dpi = 600,
  width = 12,
  height = 18
)


################################################################################
# Supplementary Figure 7
################################################################################

plots_s7 <- list()

for(clock in c("Proteomics-MR-Clock1",
               "Proteomics-MR-Clock2")){

  for(sex in c("Men",
               "Women",
               "Overall")){

    plots_s7[[length(plots_s7)+1]] <-
      plot_cor(ukb_derivation,
               clock,
               sex,
               paste0(sex,"\nUKB 70%"))

    plots_s7[[length(plots_s7)+1]] <-
      plot_cor(ukb_internal,
               clock,
               sex,
               paste0(sex,"\nUKB 30%"))

    plots_s7[[length(plots_s7)+1]] <-
      plot_cor(esther_HT,
               clock,
               sex,
               paste0(sex,"\nESTHER"))

  }

}


FigureS7 <- wrap_plots(
  plots_s7,
  ncol = 3
)


ggsave(
  file.path(
    dir_figures,
    "Supplementary_Figure7.pdf"
  ),
  FigureS7,
  width = 12,
  height = 18
)


ggsave(
  file.path(
    dir_figures,
    "Supplementary_Figure7.png"
  ),
  FigureS7,
  dpi = 600,
  width = 12,
  height = 18
)


cat("\nSupplementary Figures 6–7 completed.\n")
