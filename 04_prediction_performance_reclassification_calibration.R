################################################################################
# 04_prediction_performance_reclassification_calibration.R
#
# Inflammatory and proteomic biomarkers improve prediction of mortality risk
#
# Manuscript section:
#   Inflammatory and proteomic biomarkers improve prediction of mortality risk
#
# Purpose:
#   1. Read organized datasets from Part 1 and selected proteins from Part 2.
#   2. Apply derivation-set protein standardization parameters to internal and
#      external validation datasets.
#   3. Generate Supplementary Figure 2: Spearman correlations between selected
#      proteins and conventional risk factors.
#   4. Fit panel- and sex-specific Cox prediction models in the UKB derivation set:
#        a) conventional risk factors only
#        b) selected proteins only
#        c) conventional risk factors + selected proteins
#      and evaluate them without re-estimation in internal and external validation.
#   5. Generate:
#        Table 2: 10-year mortality prediction
#        Supplementary Table 3: 5-year mortality prediction
#        Supplementary Figure 3: individual protein incremental discrimination
#        Supplementary Tables 4–5: sex-specific 10-year prediction
#        Supplementary Figures 4–5: calibration plots
#        Supplementary Tables 6–9: categorical NRI cross-tabulations
#
# Statistical approach:
#   - Cox proportional hazards models are fit in the UKB derivation set.
#   - Predicted absolute risk at 5 or 10 years is calculated from Cox baseline
#     hazards and linear predictors.
#   - C-statistics are calculated as time-horizon ROC AUCs for 5- or 10-year
#     all-cause mortality using predicted absolute risks.
#   - Continuous NRI, categorical NRI, and IDI compare combined models against
#     conventional risk factor models.
#   - Bootstrap confidence intervals are provided for C-statistics, NRI, and IDI.
#
# Notes:
#   - This script intentionally avoids fitting models in validation datasets.
#   - Overall analyses use sex-specific biomarker models and combine predicted
#     risks across men and women.
#   - Conventional risk factor models for the overall population include sex.
#   - Sex-stratified models exclude sex.
################################################################################


# ==============================================================================
# 0. Packages
# ==============================================================================

required_packages <- c(
  "dplyr", "tidyr", "purrr", "stringr", "readr", "tibble",
  "survival", "broom", "pROC", "ggplot2", "ggrepel",
  "scales", "patchwork"
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

# ---- Inputs from previous parts ----
file_ukb_derivation <- "results/analysis_datasets/ukb_derivation.rds"
file_ukb_internal   <- "results/analysis_datasets/ukb_internal_validation.rds"
file_esther_ti      <- "results/analysis_datasets/esther_target_inflammation.rds"
file_esther_ht      <- "results/analysis_datasets/esther_explore_ht.rds"

file_selected_proteins <- "results/selected_proteins/selected_proteins_all_panels.csv"

file_scaling_ti <- "results/analysis_datasets/protein_standardization_parameters_target_inflammation.csv"
file_scaling_ex <- "results/analysis_datasets/protein_standardization_parameters_explore.csv"

# ---- Output directories ----
dir_results <- "results"
dir_tables  <- file.path(dir_results, "tables")
dir_figures <- file.path(dir_results, "figures")
dir_models  <- file.path(dir_results, "models")

dir.create(dir_tables,  recursive = TRUE, showWarnings = FALSE)
dir.create(dir_figures, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_models,  recursive = TRUE, showWarnings = FALSE)

# ---- Core variables ----
sex_var   <- "sex_share"
time_var  <- "followup_time"
event_var <- "death_10y"

male_labels   <- c("Male", "Men", "M", "1", 1)
female_labels <- c("Female", "Women", "F", "0", 0, "2", 2)

# ---- Conventional risk factors ----
risk_factors <- c(
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

continuous_covariates <- c("age_share")

risk_factor_labels <- c(
  age_share          = "Age",
  sex_share          = "Sex",
  education_share    = "Education",
  smoking_share      = "Smoking",
  activity_share     = "Physical activity",
  alcohol_share      = "Alcohol",
  BMI_share          = "BMI",
  hypertension_share = "Hypertension",
  diabetes_share     = "Diabetes",
  dyslipidemia_share = "Dyslipidemia",
  CVD_share          = "CVD",
  cancer_share       = "Cancer"
)

# ---- Panels and datasets ----
panel_names <- c("target_inflammation", "explore")

panel_labels <- c(
  target_inflammation = "Inflammation panel",
  explore = "Proteomic panel"
)

dataset_order <- c("derivation", "internal_validation", "external_validation")

dataset_labels <- c(
  derivation = "UK Biobank derivation set",
  internal_validation = "UK Biobank internal validation set",
  external_validation = "ESTHER"
)

# ---- Risk prediction parameters ----
horizons <- c(10, 5)
risk_categories <- c(0, 0.05, 0.10, 1)
n_bootstrap_metrics <- 1000
bootstrap_seed <- 20240618


# ==============================================================================
# 2. General helper functions
# ==============================================================================

require_file <- function(path, description) {
  if (!file.exists(path)) {
    stop(description, " not found: ", path)
  }
}

format_count <- function(x) {
  formatC(x, format = "d", big.mark = ",")
}

standardize_sex <- function(x) {
  x_chr <- as.character(x)
  dplyr::case_when(
    x_chr %in% as.character(male_labels) ~ "Men",
    x_chr %in% as.character(female_labels) ~ "Women",
    TRUE ~ x_chr
  )
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

make_horizon_endpoint <- function(data, horizon) {
  data %>%
    dplyr::mutate(
      horizon_time = pmin(.data[[time_var]], horizon),
      horizon_event = as.integer(.data[[event_var]] == 1 & .data[[time_var]] <= horizon)
    )
}

keep_existing_vars <- function(vars, data) {
  intersect(vars, names(data))
}

read_selected_proteins <- function(path) {
  require_file(path, "Selected protein file from Part 2")

  selected <- readr::read_csv(path, show_col_types = FALSE)

  selected %>%
    dplyr::mutate(
      panel = normalize_panel_name(panel),
      sex = normalize_sex_name(sex),
      protein = as.character(protein)
    ) %>%
    dplyr::filter(panel %in% panel_names, sex %in% c("Men", "Women"), protein != "") %>%
    dplyr::distinct(panel, sex, protein)
}

get_selected_proteins <- function(selected, panel_name, sex_level = NULL) {
  out <- selected %>% dplyr::filter(panel == panel_name)
  if (!is.null(sex_level)) out <- out %>% dplyr::filter(sex == sex_level)
  out %>% dplyr::pull(protein) %>% unique()
}

apply_protein_standardization <- function(data, scaling_file, panel_name) {
  require_file(scaling_file, paste0("Protein scaling parameters for ", panel_name))

  scaling <- readr::read_csv(scaling_file, show_col_types = FALSE)
  proteins <- intersect(scaling$protein, names(data))

  out <- data
  for (p in proteins) {
    center_p <- scaling$center[match(p, scaling$protein)]
    scale_p  <- scaling$scale[match(p, scaling$protein)]
    if (is.na(scale_p) || scale_p == 0) scale_p <- 1
    out[[p]] <- (as.numeric(out[[p]]) - center_p) / scale_p
  }

  out
}

prepare_covariates <- function(data, include_sex = TRUE, training_levels = NULL) {
  covars <- risk_factors
  if (!include_sex) covars <- setdiff(covars, "sex_share")
  covars <- keep_existing_vars(covars, data)

  for (v in covars) {
    if (v %in% continuous_covariates) {
      data[[v]] <- as.numeric(data[[v]])
      if (v == "age_share") data[[v]] <- data[[v]] / 10
    } else {
      data[[v]] <- as.factor(data[[v]])
    }
  }

  if ("BMI_share" %in% names(data)) {
    data$BMI_share <- as.factor(data$BMI_share)
    if (is.null(training_levels) && nlevels(data$BMI_share) >= 2) {
      data$BMI_share <- stats::relevel(data$BMI_share, ref = levels(data$BMI_share)[2])
    }
  }

  if ("alcohol_share" %in% names(data)) {
    data$alcohol_share <- as.factor(data$alcohol_share)
    if (is.null(training_levels) && nlevels(data$alcohol_share) >= 2) {
      data$alcohol_share <- stats::relevel(data$alcohol_share, ref = tail(levels(data$alcohol_share), 1))
    }
  }

  if (!is.null(training_levels)) {
    for (v in names(training_levels)) {
      if (v %in% names(data)) {
        data[[v]] <- factor(as.character(data[[v]]), levels = training_levels[[v]])
      }
    }
  }

  data <- droplevels(data)
  list(data = data, covariates = covars)
}

capture_factor_levels <- function(data, covariates) {
  covariates <- keep_existing_vars(covariates, data)
  factor_covars <- covariates[vapply(data[covariates], is.factor, logical(1))]
  purrr::map(data[factor_covars], levels)
}

prepare_analysis_data <- function(data, proteins, include_sex, horizon,
                                  training_levels = NULL, sex_level = NULL) {
  dat <- data %>%
    dplyr::mutate(sex_standardized = standardize_sex(.data[[sex_var]]))

  if (!is.null(sex_level)) {
    dat <- dat %>% dplyr::filter(sex_standardized == sex_level)
  }

  dat <- make_horizon_endpoint(dat, horizon)

  for (p in proteins) {
    if (p %in% names(dat)) dat[[p]] <- as.numeric(dat[[p]])
  }

  prep <- prepare_covariates(dat, include_sex = include_sex, training_levels = training_levels)
  dat <- prep$data
  covars <- prep$covariates

  required_vars <- c("horizon_time", "horizon_event", covars, proteins)
  required_vars <- keep_existing_vars(required_vars, dat)

  dat <- dat[, unique(required_vars), drop = FALSE] %>%
    dplyr::filter(!is.na(horizon_time), horizon_time > 0, !is.na(horizon_event)) %>%
    tidyr::drop_na() %>%
    droplevels()

  covars_final <- covars[purrr::map_lgl(covars, function(v) {
    if (!v %in% names(dat)) return(FALSE)
    if (is.factor(dat[[v]])) nlevels(dat[[v]]) >= 2 else stats::sd(dat[[v]], na.rm = TRUE) > 0
  })]

  proteins_final <- proteins[proteins %in% names(dat)]

  list(data = dat, covariates = covars_final, proteins = proteins_final)
}

fit_cox_prediction_model <- function(data, predictors) {
  predictors <- unique(predictors)
  predictors <- keep_existing_vars(predictors, data)

  if (length(predictors) == 0) stop("No predictors available for Cox model.")

  f <- as.formula(
    paste0("survival::Surv(horizon_time, horizon_event) ~ ",
           paste(predictors, collapse = " + "))
  )

  survival::coxph(
    f,
    data = data,
    x = TRUE,
    y = TRUE,
    model = TRUE,
    singular.ok = TRUE,
    control = survival::coxph.control(timefix = FALSE)
  )
}

predict_absolute_risk <- function(fit, newdata, horizon) {
  bh <- survival::basehaz(fit, centered = FALSE)

  if (nrow(bh) == 0) stop("Baseline hazard is empty.")

  hazard_t <- bh$hazard[max(which(bh$time <= horizon))]
  if (length(hazard_t) == 0 || is.na(hazard_t)) hazard_t <- 0

  lp <- as.numeric(stats::predict(fit, newdata = newdata, type = "lp"))
  risk <- 1 - exp(-hazard_t * exp(lp))
  pmin(pmax(risk, 0), 1)
}

safe_auc <- function(event, risk) {
  if (length(unique(event)) < 2) return(NA_real_)
  roc_obj <- tryCatch(
    pROC::roc(response = event, predictor = risk, quiet = TRUE, direction = "<"),
    error = function(e) NULL
  )
  if (is.null(roc_obj)) return(NA_real_)
  as.numeric(pROC::auc(roc_obj))
}

bootstrap_auc_ci <- function(event, risk, n_boot = n_bootstrap_metrics, seed = bootstrap_seed) {
  set.seed(seed)
  complete <- is.finite(event) & is.finite(risk)
  event <- event[complete]
  risk <- risk[complete]

  point <- safe_auc(event, risk)

  boot_values <- replicate(n_boot, {
    idx <- sample(seq_along(event), replace = TRUE)
    safe_auc(event[idx], risk[idx])
  })

  ci <- stats::quantile(boot_values, probs = c(0.025, 0.975), na.rm = TRUE)

  tibble::tibble(
    C_statistic = point,
    CI_lower = as.numeric(ci[1]),
    CI_upper = as.numeric(ci[2])
  )
}

calculate_binary_nri_idi <- function(event, risk_old, risk_new, categories = risk_categories) {
  complete <- is.finite(event) & is.finite(risk_old) & is.finite(risk_new)
  event <- as.integer(event[complete] == 1)
  risk_old <- risk_old[complete]
  risk_new <- risk_new[complete]

  if (length(unique(event)) < 2) {
    return(tibble::tibble(
      continuous_NRI = NA_real_,
      categorical_NRI = NA_real_,
      IDI = NA_real_,
      event_up = NA_real_,
      event_down = NA_real_,
      nonevent_up = NA_real_,
      nonevent_down = NA_real_
    ))
  }

  event_idx <- event == 1
  nonevent_idx <- event == 0

  delta <- risk_new - risk_old

  event_up <- mean(delta[event_idx] > 0)
  event_down <- mean(delta[event_idx] < 0)
  nonevent_up <- mean(delta[nonevent_idx] > 0)
  nonevent_down <- mean(delta[nonevent_idx] < 0)

  continuous_nri <- (event_up - event_down) + (nonevent_down - nonevent_up)

  old_cat <- cut(risk_old, breaks = categories, include.lowest = TRUE, right = FALSE)
  new_cat <- cut(risk_new, breaks = categories, include.lowest = TRUE, right = FALSE)

  old_num <- as.numeric(old_cat)
  new_num <- as.numeric(new_cat)

  event_up_cat <- mean(new_num[event_idx] > old_num[event_idx], na.rm = TRUE)
  event_down_cat <- mean(new_num[event_idx] < old_num[event_idx], na.rm = TRUE)
  nonevent_up_cat <- mean(new_num[nonevent_idx] > old_num[nonevent_idx], na.rm = TRUE)
  nonevent_down_cat <- mean(new_num[nonevent_idx] < old_num[nonevent_idx], na.rm = TRUE)

  categorical_nri <- (event_up_cat - event_down_cat) + (nonevent_down_cat - nonevent_up_cat)

  idi <- (mean(risk_new[event_idx]) - mean(risk_new[nonevent_idx])) -
    (mean(risk_old[event_idx]) - mean(risk_old[nonevent_idx]))

  tibble::tibble(
    continuous_NRI = continuous_nri,
    categorical_NRI = categorical_nri,
    IDI = idi,
    event_up = event_up,
    event_down = event_down,
    nonevent_up = nonevent_up,
    nonevent_down = nonevent_down
  )
}

bootstrap_nri_idi_ci <- function(event, risk_old, risk_new, n_boot = n_bootstrap_metrics,
                                 seed = bootstrap_seed) {
  set.seed(seed)
  complete <- is.finite(event) & is.finite(risk_old) & is.finite(risk_new)

  event <- event[complete]
  risk_old <- risk_old[complete]
  risk_new <- risk_new[complete]

  point <- calculate_binary_nri_idi(event, risk_old, risk_new)

  boot <- replicate(n_boot, {
    idx <- sample(seq_along(event), replace = TRUE)
    calculate_binary_nri_idi(event[idx], risk_old[idx], risk_new[idx]) %>%
      dplyr::select(continuous_NRI, categorical_NRI, IDI) %>%
      as.numeric()
  })

  boot <- t(boot)
  ci <- apply(boot, 2, stats::quantile, probs = c(0.025, 0.975), na.rm = TRUE)
  colnames(ci) <- c("continuous_NRI", "categorical_NRI", "IDI")

  tibble::tibble(
    continuous_NRI = point$continuous_NRI,
    continuous_NRI_lower = ci[1, "continuous_NRI"],
    continuous_NRI_upper = ci[2, "continuous_NRI"],
    categorical_NRI = point$categorical_NRI,
    categorical_NRI_lower = ci[1, "categorical_NRI"],
    categorical_NRI_upper = ci[2, "categorical_NRI"],
    IDI = point$IDI,
    IDI_lower = ci[1, "IDI"],
    IDI_upper = ci[2, "IDI"]
  )
}

make_category_crosstab <- function(event, risk_old, risk_new, categories = risk_categories) {
  complete <- is.finite(event) & is.finite(risk_old) & is.finite(risk_new)

  df <- tibble::tibble(
    event = as.integer(event[complete] == 1),
    old_category = cut(risk_old[complete], breaks = categories, include.lowest = TRUE, right = FALSE),
    new_category = cut(risk_new[complete], breaks = categories, include.lowest = TRUE, right = FALSE)
  )

  df %>%
    dplyr::count(event, old_category, new_category, name = "n") %>%
    dplyr::group_by(event) %>%
    dplyr::mutate(percent = 100 * n / sum(n)) %>%
    dplyr::ungroup()
}


# ==============================================================================
# 3. Load data and apply protein standardization
# ==============================================================================

for (f in c(file_ukb_derivation, file_ukb_internal, file_esther_ti,
            file_esther_ht, file_selected_proteins, file_scaling_ti, file_scaling_ex)) {
  require_file(f, f)
}

ukb_derivation_raw <- readr::read_rds(file_ukb_derivation)
ukb_internal_raw   <- readr::read_rds(file_ukb_internal)
esther_ti_raw      <- readr::read_rds(file_esther_ti)
esther_ht_raw      <- readr::read_rds(file_esther_ht)

selected_proteins <- read_selected_proteins(file_selected_proteins)

datasets_by_panel <- list(
  target_inflammation = list(
    derivation = apply_protein_standardization(ukb_derivation_raw, file_scaling_ti, "target_inflammation"),
    internal_validation = apply_protein_standardization(ukb_internal_raw, file_scaling_ti, "target_inflammation"),
    external_validation = apply_protein_standardization(esther_ti_raw, file_scaling_ti, "target_inflammation")
  ),
  explore = list(
    derivation = apply_protein_standardization(ukb_derivation_raw, file_scaling_ex, "explore"),
    internal_validation = apply_protein_standardization(ukb_internal_raw, file_scaling_ex, "explore"),
    external_validation = apply_protein_standardization(esther_ht_raw, file_scaling_ex, "explore")
  )
)


# ==============================================================================
# 4. Supplementary Figure 2: correlation heatmaps
# ==============================================================================

make_spearman_correlation_table <- function(data, proteins, covariates) {
  vars <- unique(c(proteins, covariates))
  vars <- keep_existing_vars(vars, data)

  dat <- data[, vars, drop = FALSE] %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), ~ as.numeric(as.character(.x))))

  cor_mat <- stats::cor(dat, method = "spearman", use = "pairwise.complete.obs")

  as.data.frame(as.table(cor_mat)) %>%
    tibble::as_tibble() %>%
    dplyr::rename(Variable_1 = Var1, Variable_2 = Var2, Spearman_r = Freq)
}

plot_correlation_heatmap <- function(cor_table, title) {
  ggplot2::ggplot(cor_table, ggplot2::aes(x = Variable_1, y = Variable_2, fill = Spearman_r)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.15) +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-1, 1),
      name = "Spearman r"
    ) +
    ggplot2::coord_fixed() +
    ggplot2::labs(x = NULL, y = NULL, title = title) +
    ggplot2::theme_classic(base_size = 9) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
      axis.text.y = ggplot2::element_text(size = 7),
      axis.line = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", size = 11)
    )
}

corr_outputs <- purrr::map(panel_names, function(panel_name) {
  proteins_union <- get_selected_proteins(selected_proteins, panel_name)
  covariates <- keep_existing_vars(risk_factors, datasets_by_panel[[panel_name]]$derivation)

  cor_table <- make_spearman_correlation_table(
    data = datasets_by_panel[[panel_name]]$derivation,
    proteins = proteins_union,
    covariates = covariates
  ) %>%
    dplyr::mutate(panel = panel_name)

  readr::write_csv(
    cor_table,
    file.path(dir_tables, paste0("Supplementary_Figure_2_correlation_", panel_name, ".csv"))
  )

  p <- plot_correlation_heatmap(cor_table, panel_labels[[panel_name]])

  ggplot2::ggsave(
    file.path(dir_figures, paste0("Supplementary_Figure_2_correlation_", panel_name, ".pdf")),
    p, width = 8, height = 7
  )
  ggplot2::ggsave(
    file.path(dir_figures, paste0("Supplementary_Figure_2_correlation_", panel_name, ".png")),
    p, width = 8, height = 7, dpi = 600
  )

  cor_table
})

cor_all <- dplyr::bind_rows(corr_outputs)

cor_summary <- cor_all %>%
  dplyr::mutate(abs_r = abs(Spearman_r)) %>%
  dplyr::group_by(panel) %>%
  dplyr::summarise(max_abs_spearman = max(abs_r[Variable_1 != Variable_2], na.rm = TRUE), .groups = "drop")

readr::write_csv(cor_summary, file.path(dir_tables, "Supplementary_Figure_2_max_correlations.csv"))


# ==============================================================================
# 5. Fit derivation Cox models and evaluate prediction performance
# ==============================================================================

fit_prediction_models_for_group <- function(panel_name, sex_level, horizon) {
  derivation <- datasets_by_panel[[panel_name]]$derivation
  proteins <- get_selected_proteins(selected_proteins, panel_name, sex_level)

  train_cov <- prepare_analysis_data(
    data = derivation,
    proteins = character(0),
    include_sex = FALSE,
    horizon = horizon,
    sex_level = sex_level
  )

  training_levels <- capture_factor_levels(train_cov$data, train_cov$covariates)

  train_full <- prepare_analysis_data(
    data = derivation,
    proteins = proteins,
    include_sex = FALSE,
    horizon = horizon,
    training_levels = training_levels,
    sex_level = sex_level
  )

  fit_covariates <- fit_cox_prediction_model(train_cov$data, train_cov$covariates)
  fit_proteins <- fit_cox_prediction_model(train_full$data, train_full$proteins)
  fit_combined <- fit_cox_prediction_model(train_full$data, c(train_full$covariates, train_full$proteins))

  list(
    panel = panel_name,
    sex = sex_level,
    horizon = horizon,
    proteins = proteins,
    covariates = train_full$covariates,
    training_levels = training_levels,
    fit_covariates = fit_covariates,
    fit_proteins = fit_proteins,
    fit_combined = fit_combined
  )
}

predict_group_models <- function(model_object, data, dataset_name) {
  sex_level <- model_object$sex
  horizon <- model_object$horizon

  dat_cov <- prepare_analysis_data(
    data = data,
    proteins = character(0),
    include_sex = FALSE,
    horizon = horizon,
    training_levels = model_object$training_levels,
    sex_level = sex_level
  )$data

  dat_full <- prepare_analysis_data(
    data = data,
    proteins = model_object$proteins,
    include_sex = FALSE,
    horizon = horizon,
    training_levels = model_object$training_levels,
    sex_level = sex_level
  )$data

  # Use complete cases shared by the combined model.
  shared_n <- nrow(dat_full)

  tibble::tibble(
    panel = model_object$panel,
    sex = sex_level,
    dataset = dataset_name,
    horizon = horizon,
    time = dat_full$horizon_time,
    event = dat_full$horizon_event,
    risk_covariates = predict_absolute_risk(model_object$fit_covariates, dat_full, horizon),
    risk_proteins = predict_absolute_risk(model_object$fit_proteins, dat_full, horizon),
    risk_combined = predict_absolute_risk(model_object$fit_combined, dat_full, horizon),
    n = shared_n
  )
}

evaluate_prediction_table <- function(prediction_df) {
  auc_cov <- bootstrap_auc_ci(prediction_df$event, prediction_df$risk_covariates)
  auc_pro <- bootstrap_auc_ci(prediction_df$event, prediction_df$risk_proteins)
  auc_com <- bootstrap_auc_ci(prediction_df$event, prediction_df$risk_combined)

  nri <- bootstrap_nri_idi_ci(
    event = prediction_df$event,
    risk_old = prediction_df$risk_covariates,
    risk_new = prediction_df$risk_combined
  )

  tibble::tibble(
    C_covariates = auc_cov$C_statistic,
    C_covariates_lower = auc_cov$CI_lower,
    C_covariates_upper = auc_cov$CI_upper,
    C_proteins = auc_pro$C_statistic,
    C_proteins_lower = auc_pro$CI_lower,
    C_proteins_upper = auc_pro$CI_upper,
    C_combined = auc_com$C_statistic,
    C_combined_lower = auc_com$CI_lower,
    C_combined_upper = auc_com$CI_upper,
    Delta_C_combined_vs_covariates = auc_com$C_statistic - auc_cov$C_statistic
  ) %>%
    dplyr::bind_cols(nri)
}

all_model_objects <- list()
all_predictions <- list()

for (panel_name in panel_names) {
  for (horizon in horizons) {
    for (sex_level in c("Men", "Women")) {
      model_key <- paste(panel_name, horizon, sex_level, sep = "_")
      model_object <- fit_prediction_models_for_group(panel_name, sex_level, horizon)
      all_model_objects[[model_key]] <- model_object

      for (dataset_name in dataset_order) {
        pred <- predict_group_models(
          model_object = model_object,
          data = datasets_by_panel[[panel_name]][[dataset_name]],
          dataset_name = dataset_name
        )
        all_predictions[[paste(model_key, dataset_name, sep = "_")]] <- pred
      }
    }
  }
}

readr::write_rds(all_model_objects, file.path(dir_models, "prediction_cox_models_all_panels.rds"))

prediction_long <- dplyr::bind_rows(all_predictions)
readr::write_csv(prediction_long, file.path(dir_tables, "predicted_risks_all_models.csv"))


# ==============================================================================
# 6. Table 2 and Supplementary Table 3: overall prediction performance
# ==============================================================================

# Overall prediction combines sex-specific risks within each panel.
overall_prediction <- prediction_long %>%
  dplyr::group_by(panel, dataset, horizon) %>%
  dplyr::group_modify(~ .x %>% dplyr::ungroup()) %>%
  dplyr::ungroup()

overall_results <- overall_prediction %>%
  dplyr::group_by(panel, dataset, horizon) %>%
  dplyr::group_modify(~ evaluate_prediction_table(.x)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    population = "Overall",
    panel_label = panel_labels[panel],
    dataset_label = dataset_labels[dataset]
  ) %>%
  dplyr::select(population, panel, panel_label, dataset, dataset_label, horizon, dplyr::everything())

table2 <- overall_results %>% dplyr::filter(horizon == 10)
supp_table3 <- overall_results %>% dplyr::filter(horizon == 5)

readr::write_csv(table2, file.path(dir_tables, "Table_2_prediction_performance_10_year.csv"))
readr::write_csv(supp_table3, file.path(dir_tables, "Supplementary_Table_3_prediction_performance_5_year.csv"))


# ==============================================================================
# 7. Supplementary Tables 4–5: sex-specific 10-year prediction performance
# ==============================================================================

sex_specific_results <- prediction_long %>%
  dplyr::filter(horizon == 10) %>%
  dplyr::group_by(panel, sex, dataset, horizon) %>%
  dplyr::group_modify(~ evaluate_prediction_table(.x)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    panel_label = panel_labels[panel],
    dataset_label = dataset_labels[dataset]
  )

supp_table4_men <- sex_specific_results %>% dplyr::filter(sex == "Men")
supp_table5_women <- sex_specific_results %>% dplyr::filter(sex == "Women")

readr::write_csv(supp_table4_men, file.path(dir_tables, "Supplementary_Table_4_prediction_performance_men_10_year.csv"))
readr::write_csv(supp_table5_women, file.path(dir_tables, "Supplementary_Table_5_prediction_performance_women_10_year.csv"))


# ==============================================================================
# 8. Supplementary Table 2: model coefficients
# ==============================================================================

extract_model_coefficients <- function(model_object) {
  extract_one <- function(fit, model_type) {
    broom::tidy(fit) %>%
      dplyr::transmute(
        panel = model_object$panel,
        sex = model_object$sex,
        horizon = model_object$horizon,
        model = model_type,
        term = term,
        coefficient = estimate,
        SE = std.error,
        statistic = statistic,
        P_value = p.value
      )
  }

  dplyr::bind_rows(
    extract_one(model_object$fit_covariates, "conventional_risk_factors"),
    extract_one(model_object$fit_proteins, "proteins_only"),
    extract_one(model_object$fit_combined, "combined")
  )
}

model_coefficients <- purrr::map_dfr(all_model_objects, extract_model_coefficients)

readr::write_csv(
  model_coefficients,
  file.path(dir_tables, "Supplementary_Table_2_panel_sex_specific_model_coefficients.csv")
)


# ==============================================================================
# 9. Supplementary Figure 3: individual protein incremental discrimination
# ==============================================================================

fit_individual_protein_model <- function(panel_name, sex_level, horizon, protein) {
  derivation <- datasets_by_panel[[panel_name]]$derivation

  train_cov <- prepare_analysis_data(
    data = derivation,
    proteins = character(0),
    include_sex = FALSE,
    horizon = horizon,
    sex_level = sex_level
  )

  training_levels <- capture_factor_levels(train_cov$data, train_cov$covariates)

  train_one <- prepare_analysis_data(
    data = derivation,
    proteins = protein,
    include_sex = FALSE,
    horizon = horizon,
    training_levels = training_levels,
    sex_level = sex_level
  )

  fit_one <- fit_cox_prediction_model(train_one$data, c(train_one$covariates, protein))
  fit_cov <- fit_cox_prediction_model(train_cov$data, train_cov$covariates)

  list(
    panel = panel_name,
    sex = sex_level,
    horizon = horizon,
    protein = protein,
    training_levels = training_levels,
    covariates = train_one$covariates,
    fit_covariates = fit_cov,
    fit_one = fit_one
  )
}

predict_individual_protein_model <- function(model_object, data, dataset_name) {
  dat <- prepare_analysis_data(
    data = data,
    proteins = model_object$protein,
    include_sex = FALSE,
    horizon = model_object$horizon,
    training_levels = model_object$training_levels,
    sex_level = model_object$sex
  )$data

  tibble::tibble(
    panel = model_object$panel,
    sex = model_object$sex,
    dataset = dataset_name,
    horizon = model_object$horizon,
    protein = model_object$protein,
    event = dat$horizon_event,
    risk_covariates = predict_absolute_risk(model_object$fit_covariates, dat, model_object$horizon),
    risk_single_protein = predict_absolute_risk(model_object$fit_one, dat, model_object$horizon)
  )
}

individual_results <- list()

for (panel_name in panel_names) {
  for (sex_level in c("Men", "Women")) {
    proteins <- get_selected_proteins(selected_proteins, panel_name, sex_level)
    for (protein in proteins) {
      one_model <- fit_individual_protein_model(panel_name, sex_level, horizon = 10, protein = protein)

      for (dataset_name in dataset_order) {
        pred_one <- predict_individual_protein_model(
          one_model,
          data = datasets_by_panel[[panel_name]][[dataset_name]],
          dataset_name = dataset_name
        )

        auc_cov <- bootstrap_auc_ci(pred_one$event, pred_one$risk_covariates, n_boot = 200)
        auc_one <- bootstrap_auc_ci(pred_one$event, pred_one$risk_single_protein, n_boot = 200)

        individual_results[[paste(panel_name, sex_level, protein, dataset_name, sep = "_")]] <-
          tibble::tibble(
            panel = panel_name,
            sex = sex_level,
            dataset = dataset_name,
            protein = protein,
            C_covariates = auc_cov$C_statistic,
            C_covariates_plus_protein = auc_one$C_statistic,
            Delta_C = auc_one$C_statistic - auc_cov$C_statistic
          )
      }
    }
  }
}

individual_results <- dplyr::bind_rows(individual_results)
readr::write_csv(individual_results, file.path(dir_tables, "Supplementary_Figure_3_individual_protein_incremental_discrimination.csv"))

plot_individual_incremental <- function(data) {
  data %>%
    dplyr::mutate(
      panel_label = panel_labels[panel],
      dataset_label = dataset_labels[dataset],
      protein = forcats::fct_reorder(protein, Delta_C)
    ) %>%
    ggplot2::ggplot(ggplot2::aes(x = Delta_C, y = protein)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    ggplot2::geom_point(size = 1.8) +
    ggplot2::facet_grid(panel_label + sex ~ dataset_label, scales = "free_y", space = "free_y") +
    ggplot2::labs(
      x = "Incremental C-statistic after adding one protein",
      y = NULL
    ) +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "grey95", color = NA),
      strip.text = ggplot2::element_text(face = "bold"),
      axis.text.y = ggplot2::element_text(size = 7)
    )
}

p_individual <- plot_individual_incremental(individual_results)

ggplot2::ggsave(
  file.path(dir_figures, "Supplementary_Figure_3_individual_protein_incremental_discrimination.pdf"),
  p_individual,
  width = 12,
  height = 9
)

ggplot2::ggsave(
  file.path(dir_figures, "Supplementary_Figure_3_individual_protein_incremental_discrimination.png"),
  p_individual,
  width = 12,
  height = 9,
  dpi = 600
)


# ==============================================================================
# 10. Supplementary Figures 4–5: calibration plots
# ==============================================================================

make_calibration_data <- function(prediction_df, risk_col = "risk_combined", n_groups = 10) {
  prediction_df %>%
    dplyr::mutate(
      predicted_risk = .data[[risk_col]],
      risk_decile = dplyr::ntile(predicted_risk, n_groups)
    ) %>%
    dplyr::group_by(panel, sex, dataset, horizon, risk_decile) %>%
    dplyr::summarise(
      mean_predicted = mean(predicted_risk, na.rm = TRUE),
      observed = mean(event == 1, na.rm = TRUE),
      n = dplyr::n(),
      .groups = "drop"
    )
}

plot_calibration_panel <- function(calibration_data, panel_name) {
  calibration_data %>%
    dplyr::filter(panel == panel_name, horizon == 10) %>%
    dplyr::mutate(
      dataset_label = dataset_labels[dataset],
      sex = factor(sex, levels = c("Men", "Women"))
    ) %>%
    ggplot2::ggplot(ggplot2::aes(x = mean_predicted, y = observed)) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey60") +
    ggplot2::geom_point(size = 1.8) +
    ggplot2::geom_smooth(method = "loess", se = FALSE, linewidth = 0.5, formula = y ~ x) +
    ggplot2::facet_grid(sex ~ dataset_label) +
    ggplot2::scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    ggplot2::labs(
      x = "Predicted 10-year mortality risk",
      y = "Observed 10-year mortality proportion",
      title = panel_labels[[panel_name]]
    ) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      strip.background = ggplot2::element_rect(fill = "grey95", color = NA),
      strip.text = ggplot2::element_text(face = "bold")
    )
}

calibration_data <- make_calibration_data(prediction_long, risk_col = "risk_combined")
readr::write_csv(calibration_data, file.path(dir_tables, "Supplementary_Figures_4_5_calibration_data.csv"))

p_cal_ti <- plot_calibration_panel(calibration_data, "target_inflammation")
p_cal_ex <- plot_calibration_panel(calibration_data, "explore")

ggplot2::ggsave(file.path(dir_figures, "Supplementary_Figure_4_calibration_target_inflammation.pdf"),
                p_cal_ti, width = 9, height = 5.5)
ggplot2::ggsave(file.path(dir_figures, "Supplementary_Figure_4_calibration_target_inflammation.png"),
                p_cal_ti, width = 9, height = 5.5, dpi = 600)

ggplot2::ggsave(file.path(dir_figures, "Supplementary_Figure_5_calibration_explore.pdf"),
                p_cal_ex, width = 9, height = 5.5)
ggplot2::ggsave(file.path(dir_figures, "Supplementary_Figure_5_calibration_explore.png"),
                p_cal_ex, width = 9, height = 5.5, dpi = 600)


# ==============================================================================
# 11. Supplementary Tables 6–9: categorical NRI cross-tabulations
# ==============================================================================

make_categorical_nri_tables <- function(pred_df, population_label) {
  pred_df %>%
    dplyr::group_by(panel, dataset, horizon) %>%
    dplyr::group_modify(~ make_category_crosstab(
      event = .x$event,
      risk_old = .x$risk_covariates,
      risk_new = .x$risk_combined,
      categories = risk_categories
    )) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(population = population_label)
}

# Supplementary Table 6: overall 10-year mortality
supp_table6 <- overall_prediction %>%
  dplyr::filter(horizon == 10) %>%
  make_categorical_nri_tables(population_label = "Overall")

# Supplementary Table 7: overall 5-year mortality
supp_table7 <- overall_prediction %>%
  dplyr::filter(horizon == 5) %>%
  make_categorical_nri_tables(population_label = "Overall")

# Supplementary Table 8: sex-specific 10-year mortality in men
supp_table8 <- prediction_long %>%
  dplyr::filter(horizon == 10, sex == "Men") %>%
  make_categorical_nri_tables(population_label = "Men")

# Supplementary Table 9: sex-specific 10-year mortality in women
supp_table9 <- prediction_long %>%
  dplyr::filter(horizon == 10, sex == "Women") %>%
  make_categorical_nri_tables(population_label = "Women")

readr::write_csv(supp_table6, file.path(dir_tables, "Supplementary_Table_6_categorical_NRI_overall_10_year.csv"))
readr::write_csv(supp_table7, file.path(dir_tables, "Supplementary_Table_7_categorical_NRI_overall_5_year.csv"))
readr::write_csv(supp_table8, file.path(dir_tables, "Supplementary_Table_8_categorical_NRI_men_10_year.csv"))
readr::write_csv(supp_table9, file.path(dir_tables, "Supplementary_Table_9_categorical_NRI_women_10_year.csv"))


# ==============================================================================
# 12. Console summary
# ==============================================================================

cat("\nPrediction performance analysis completed successfully.\n")
cat("\nGenerated outputs:\n")
cat("  - Supplementary Figure 2 correlation heatmaps: ", dir_figures, "\n", sep = "")
cat("  - Supplementary Table 2 model coefficients: ", file.path(dir_tables, "Supplementary_Table_2_panel_sex_specific_model_coefficients.csv"), "\n", sep = "")
cat("  - Table 2 10-year prediction performance: ", file.path(dir_tables, "Table_2_prediction_performance_10_year.csv"), "\n", sep = "")
cat("  - Supplementary Table 3 5-year prediction performance: ", file.path(dir_tables, "Supplementary_Table_3_prediction_performance_5_year.csv"), "\n", sep = "")
cat("  - Supplementary Figure 3 individual protein analyses: ", dir_figures, "\n", sep = "")
cat("  - Supplementary Tables 4–5 sex-specific prediction performance: ", dir_tables, "\n", sep = "")
cat("  - Supplementary Figures 4–5 calibration plots: ", dir_figures, "\n", sep = "")
cat("  - Supplementary Tables 6–9 categorical NRI cross-tabulations: ", dir_tables, "\n", sep = "")
