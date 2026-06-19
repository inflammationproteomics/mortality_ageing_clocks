################################################################################
# 01_data_extraction_baseline_characteristics.R
#
# Data extraction and organization; Table 1; Supplementary Figure 1;
# Supplementary Table 1
#
# Manuscript section:
#   Baseline characteristics and conventional mortality risk factors
#
# Outputs:
#   results/analysis_datasets/
#     ukb_overall.rds
#     ukb_derivation.rds
#     ukb_internal_validation.rds
#     esther_target_inflammation.rds
#     esther_explore_ht.rds
#
#   results/tables/
#     Table_1_baseline_characteristics.csv
#     Supplementary_Table_1_conventional_risk_factors.csv
#     Supplementary_Table_1_conventional_risk_factors_formatted.csv
#
#   results/figures/
#     Supplementary_Figure_1_participant_selection.pdf
#     Supplementary_Figure_1_participant_selection.png
#     Supplementary_Figure_1_counts.csv
################################################################################

# ==============================================================================
# 0. Packages
# ==============================================================================

required_packages <- c(
  "dplyr", "tidyr", "purrr", "stringr", "readr",
  "tibble", "survival", "broom", "ggplot2"
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

# ---- Input RDS files ----
file_ukb_overall    <- "data/ukb_combined_analysis_imputed.rds"
file_ukb_derivation <- "data/ukb_train.rds"
file_ukb_internal   <- "data/ukb_test.rds"
file_esther_ti      <- "data/es_combined_analysis_imputed.rds"
file_esther_ht      <- "data/es_ht_combined_analysis_imputed.rds"

# ---- Output directories ----
dir_results <- "results"
dir_tables  <- file.path(dir_results, "tables")
dir_figures <- file.path(dir_results, "figures")
dir_dataout <- file.path(dir_results, "analysis_datasets")

dir.create(dir_tables,  recursive = TRUE, showWarnings = FALSE)
dir.create(dir_figures, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_dataout, recursive = TRUE, showWarnings = FALSE)

# ---- Core variables ----
id_var    <- "id"
time_var  <- "followup_year"
event_var <- "mortality10y"
alternative_time_vars <- c("followup_year", "followup")

# ---- Baseline variables and conventional risk factors ----
baseline_vars <- c(
  "age_share",
  "sex_share",
  "education_share",
  "smoking_share",
  "activity_share",
  "alcohol_share",
  "BMI_share",
  "hypertension_share",
  "diabetes_share",
  "dyslipidemia_share",
  "CVD_share",
  "cancer_share"
)

conventional_risk_factors <- baseline_vars

variable_labels <- c(
  age_share          = "Age, years",
  sex_share          = "Sex",
  education_share    = "Education",
  smoking_share      = "Smoking status",
  activity_share     = "Physical activity",
  alcohol_share      = "Alcohol consumption",
  BMI_share          = "Body mass index",
  hypertension_share = "Hypertension",
  diabetes_share     = "Diabetes",
  dyslipidemia_share = "Dyslipidemia",
  CVD_share          = "Cardiovascular disease",
  cancer_share       = "History of cancer"
)

continuous_vars <- c("age_share")

cohort_order <- c(
  "UK Biobank overall",
  "UK Biobank derivation set",
  "UK Biobank internal validation set",
  "ESTHER Target Inflammation",
  "ESTHER Explore HT"
)

# ==============================================================================
# 2. Helper functions
# ==============================================================================

read_rds_if_exists <- function(path, object_name) {
  if (!file.exists(path)) {
    warning(sprintf("File not found for %s: %s", object_name, path))
    return(NULL)
  }
  readRDS(path)
}

detect_time_var <- function(data, preferred = time_var, alternatives = alternative_time_vars) {
  candidates <- unique(c(preferred, alternatives))
  hit <- candidates[candidates %in% names(data)]
  if (length(hit) == 0) {
    stop("No follow-up time variable found. Checked: ", paste(candidates, collapse = ", "))
  }
  hit[1]
}

standardize_core_names <- function(data) {
  if (is.null(data)) return(NULL)
  detected_time <- detect_time_var(data)
  data %>%
    dplyr::rename(
      followup_time = dplyr::all_of(detected_time),
      death_10y = dplyr::all_of(event_var)
    ) %>%
    dplyr::mutate(
      followup_time = as.numeric(followup_time),
      death_10y = as.integer(as.numeric(as.character(death_10y)) == 1)
    )
}

keep_existing_vars <- function(vars, data) {
  intersect(vars, names(data))
}

format_count <- function(x) {
  formatC(x, format = "d", big.mark = ",")
}

format_mean_sd <- function(x, digits = 1) {
  sprintf(paste0("%.", digits, "f (%.", digits, "f)"),
          mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
}

format_n_percent <- function(n, denom, digits = 1) {
  sprintf("%s (%.1f%%)", format_count(n), 100 * n / denom)
}

is_binary_numeric <- function(x) {
  ux <- sort(unique(stats::na.omit(as.numeric(as.character(x)))))
  length(ux) <= 2 && all(ux %in% c(0, 1))
}

prepare_covariates <- function(data, include_sex = TRUE) {
  covars <- conventional_risk_factors
  if (!include_sex) covars <- setdiff(covars, "sex_share")
  covars <- keep_existing_vars(covars, data)

  for (v in covars) {
    if (v == "age_share") {
      data[[v]] <- as.numeric(data[[v]]) / 10
    } else {
      data[[v]] <- as.factor(data[[v]])
    }
  }

  if ("BMI_share" %in% names(data)) {
    data$BMI_share <- as.factor(data$BMI_share)
    if (nlevels(data$BMI_share) >= 2) {
      data$BMI_share <- stats::relevel(data$BMI_share, ref = levels(data$BMI_share)[2])
    }
  }

  if ("alcohol_share" %in% names(data)) {
    data$alcohol_share <- as.factor(data$alcohol_share)
    if (nlevels(data$alcohol_share) >= 2) {
      data$alcohol_share <- stats::relevel(data$alcohol_share, ref = tail(levels(data$alcohol_share), 1))
    }
  }

  list(data = droplevels(data), covariates = covars)
}

summarise_one_dataset <- function(data, cohort_name, variables = baseline_vars) {
  variables <- keep_existing_vars(variables, data)
  total_n <- nrow(data)

  purrr::map_dfr(variables, function(v) {
    x <- data[[v]]
    label <- dplyr::coalesce(unname(variable_labels[v]), v)

    if (v %in% continuous_vars || (is.numeric(x) && !is_binary_numeric(x))) {
      tibble::tibble(
        Cohort = cohort_name,
        Variable = label,
        Variable_raw = v,
        Category = "",
        Summary = format_mean_sd(as.numeric(x))
      )
    } else {
      x_fac <- as.factor(x)
      tab <- table(x_fac, useNA = "no")
      tibble::tibble(
        Cohort = cohort_name,
        Variable = label,
        Variable_raw = v,
        Category = names(tab),
        Summary = purrr::map_chr(as.integer(tab), ~ format_n_percent(.x, total_n))
      )
    }
  })
}

make_table1 <- function(dataset_list) {
  long_table <- purrr::imap_dfr(dataset_list, summarise_one_dataset)

  long_table %>%
    dplyr::mutate(
      Cohort = factor(Cohort, levels = cohort_order),
      Variable = factor(Variable, levels = unique(unname(variable_labels[baseline_vars])))
    ) %>%
    dplyr::arrange(Variable, Category, Cohort) %>%
    dplyr::select(Variable, Category, Cohort, Summary) %>%
    tidyr::pivot_wider(names_from = Cohort, values_from = Summary)
}

extract_sample_counts <- function(dataset_list) {
  purrr::imap_dfr(dataset_list, function(data, cohort_name) {
    tibble::tibble(
      Cohort = cohort_name,
      N = nrow(data),
      Deaths_10y = sum(data$death_10y == 1, na.rm = TRUE),
      Deaths_10y_percent = 100 * mean(data$death_10y == 1, na.rm = TRUE)
    )
  }) %>%
    dplyr::mutate(Cohort = factor(Cohort, levels = cohort_order)) %>%
    dplyr::arrange(Cohort)
}

plot_participant_selection <- function(counts_df) {
  counts_df <- counts_df %>%
    dplyr::mutate(Cohort = factor(Cohort, levels = cohort_order))

  ggplot2::ggplot(counts_df, ggplot2::aes(x = Cohort, y = N)) +
    ggplot2::geom_col(width = 0.65) +
    ggplot2::geom_text(
      ggplot2::aes(label = paste0("n = ", format_count(N), "\n",
                                  format_count(Deaths_10y), " deaths")),
      vjust = -0.25,
      size = 3.6
    ) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.15))) +
    ggplot2::labs(
      x = NULL,
      y = "Participants",
      title = "Participant selection and 10-year deaths"
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 30, hjust = 1),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0)
    )
}

fit_conventional_risk_cox <- function(data, cohort_name, include_sex = TRUE) {
  prep <- prepare_covariates(data, include_sex = include_sex)
  dat <- prep$data
  covars <- prep$covariates

  analysis_vars <- c("followup_time", "death_10y", covars)
  analysis_vars <- keep_existing_vars(analysis_vars, dat)

  dat <- dat[, analysis_vars, drop = FALSE] %>%
    dplyr::filter(!is.na(followup_time), followup_time > 0, !is.na(death_10y)) %>%
    tidyr::drop_na() %>%
    droplevels()

  if (sum(dat$death_10y == 1) == 0) {
    warning("No events in ", cohort_name, ". Cox model skipped.")
    return(NULL)
  }

  covars_final <- covars[purrr::map_lgl(covars, function(v) {
    if (is.factor(dat[[v]])) nlevels(dat[[v]]) >= 2 else stats::sd(dat[[v]], na.rm = TRUE) > 0
  })]

  f <- as.formula(
    paste0("survival::Surv(followup_time, death_10y) ~ ",
           paste(covars_final, collapse = " + "))
  )

  fit <- survival::coxph(
    f,
    data = dat,
    x = TRUE,
    y = TRUE,
    model = TRUE,
    singular.ok = TRUE,
    control = survival::coxph.control(timefix = FALSE)
  )

  tidy_fit <- broom::tidy(fit, conf.int = TRUE, exponentiate = TRUE)

  zph <- tryCatch(survival::cox.zph(fit), error = function(e) NULL)
  zph_table <- NULL
  if (!is.null(zph)) {
    zph_table <- as.data.frame(zph$table)
    zph_table$term <- rownames(zph_table)
  }

  tidy_fit %>%
    dplyr::transmute(
      Cohort = cohort_name,
      Term = term,
      HR = estimate,
      CI_lower = conf.low,
      CI_upper = conf.high,
      P_value = p.value,
      N = nrow(dat),
      Events = sum(dat$death_10y == 1),
      Schoenfeld_P = if (!is.null(zph_table)) zph_table$p[match(term, zph_table$term)] else NA_real_
    )
}

format_hr_ci <- function(hr, lower, upper, digits = 2) {
  sprintf(paste0("%.", digits, "f (%.", digits, "f–%.", digits, "f)"),
          hr, lower, upper)
}

make_supplementary_table1 <- function(dataset_list) {
  cox_results <- purrr::imap_dfr(dataset_list, function(data, cohort_name) {
    fit_conventional_risk_cox(
      data = data,
      cohort_name = cohort_name,
      include_sex = TRUE
    )
  })

  cox_results <- cox_results %>%
    dplyr::group_by(Cohort) %>%
    dplyr::mutate(FDR = p.adjust(P_value, method = "BH")) %>%
    dplyr::ungroup()

  formatted <- cox_results %>%
    dplyr::mutate(
      HR_95CI = format_hr_ci(HR, CI_lower, CI_upper),
      P_value = signif(P_value, 3),
      FDR = signif(FDR, 3),
      Schoenfeld_P = signif(Schoenfeld_P, 3)
    )

  list(raw = cox_results, formatted = formatted)
}

# ==============================================================================
# 3. Load and organize analysis datasets
# ==============================================================================

analysis_datasets <- list(
  "UK Biobank overall" = standardize_core_names(read_rds_if_exists(file_ukb_overall, "UKB overall")),
  "UK Biobank derivation set" = standardize_core_names(read_rds_if_exists(file_ukb_derivation, "UKB derivation")),
  "UK Biobank internal validation set" = standardize_core_names(read_rds_if_exists(file_ukb_internal, "UKB internal validation")),
  "ESTHER Target Inflammation" = standardize_core_names(read_rds_if_exists(file_esther_ti, "ESTHER Target Inflammation")),
  "ESTHER Explore HT" = standardize_core_names(read_rds_if_exists(file_esther_ht, "ESTHER Explore HT"))
)

analysis_datasets <- analysis_datasets[!vapply(analysis_datasets, is.null, logical(1))]

# Save organized datasets for downstream scripts.
if ("UK Biobank overall" %in% names(analysis_datasets)) {
  readr::write_rds(analysis_datasets[["UK Biobank overall"]],
                   file.path(dir_dataout, "ukb_overall.rds"))
}
if ("UK Biobank derivation set" %in% names(analysis_datasets)) {
  readr::write_rds(analysis_datasets[["UK Biobank derivation set"]],
                   file.path(dir_dataout, "ukb_derivation.rds"))
}
if ("UK Biobank internal validation set" %in% names(analysis_datasets)) {
  readr::write_rds(analysis_datasets[["UK Biobank internal validation set"]],
                   file.path(dir_dataout, "ukb_internal_validation.rds"))
}
if ("ESTHER Target Inflammation" %in% names(analysis_datasets)) {
  readr::write_rds(analysis_datasets[["ESTHER Target Inflammation"]],
                   file.path(dir_dataout, "esther_target_inflammation.rds"))
}
if ("ESTHER Explore HT" %in% names(analysis_datasets)) {
  readr::write_rds(analysis_datasets[["ESTHER Explore HT"]],
                   file.path(dir_dataout, "esther_explore_ht.rds"))
}

# ==============================================================================
# 4. Supplementary Figure 1: participant selection and deaths
# ==============================================================================

participant_counts <- extract_sample_counts(analysis_datasets)
readr::write_csv(participant_counts, file.path(dir_figures, "Supplementary_Figure_1_counts.csv"))

p_flow <- plot_participant_selection(participant_counts)

ggplot2::ggsave(
  filename = file.path(dir_figures, "Supplementary_Figure_1_participant_selection.pdf"),
  plot = p_flow,
  width = 8.5,
  height = 5.0
)

ggplot2::ggsave(
  filename = file.path(dir_figures, "Supplementary_Figure_1_participant_selection.png"),
  plot = p_flow,
  width = 8.5,
  height = 5.0,
  dpi = 600
)

# ==============================================================================
# 5. Table 1: baseline characteristics
# ==============================================================================

table1 <- make_table1(analysis_datasets)
readr::write_csv(table1, file.path(dir_tables, "Table_1_baseline_characteristics.csv"))

# ==============================================================================
# 6. Supplementary Table 1: conventional risk factors and 10-year mortality
# ==============================================================================

supp_table1 <- make_supplementary_table1(analysis_datasets)

readr::write_csv(
  supp_table1$raw,
  file.path(dir_tables, "Supplementary_Table_1_conventional_risk_factors.csv")
)

readr::write_csv(
  supp_table1$formatted,
  file.path(dir_tables, "Supplementary_Table_1_conventional_risk_factors_formatted.csv")
)

# ==============================================================================
# 7. Console summary
# ==============================================================================

cat("\nAnalysis completed successfully.\n")
cat("Generated outputs:\n")
cat("  - Table 1: ", file.path(dir_tables, "Table_1_baseline_characteristics.csv"), "\n", sep = "")
cat("  - Supplementary Figure 1 counts: ", file.path(dir_figures, "Supplementary_Figure_1_counts.csv"), "\n", sep = "")
cat("  - Supplementary Figure 1 PDF/PNG: ", dir_figures, "\n", sep = "")
cat("  - Supplementary Table 1: ", file.path(dir_tables, "Supplementary_Table_1_conventional_risk_factors_formatted.csv"), "\n", sep = "")
