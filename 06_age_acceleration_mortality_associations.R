################################################################################
# 06_age_acceleration_mortality_associations.R
#
# Biological age acceleration is associated with all-cause and cause-specific mortality
#
# Outputs:
#   results/tables/Table_3_age_acceleration_all_cause_mortality.csv
#   results/tables/Table_3_age_acceleration_all_cause_mortality_formatted.csv
#   results/tables/Supplementary_Table_11_age_acceleration_cause_specific_mortality.csv
#   results/tables/Supplementary_Table_11_age_acceleration_cause_specific_mortality_formatted.csv
################################################################################

library(dplyr)
library(readr)
library(tidyr)
library(purrr)
library(stringr)
library(survival)
library(broom)

dir_tables <- "results/tables"
dir_data <- "results/analysis_datasets"
dir.create(dir_tables, showWarnings = FALSE, recursive = TRUE)

file_ukb_derivation_mrage <- file.path(dir_tables, "mr_ages_derivation.csv")
file_ukb_internal_mrage <- file.path(dir_tables, "mr_ages_internal.csv")
file_esther_ti_mrage <- file.path(dir_tables, "mr_ages_esther_TI.csv")
file_esther_ht_mrage <- file.path(dir_tables, "mr_ages_esther_HT.csv")
file_cause_death <- "data/ukb_cause_specific_death.csv"

clock_map <- tibble::tibble(
  clock = c("Inflammation-MR-Clock1", "Inflammation-MR-Clock2", "Proteomics-MR-Clock1", "Proteomics-MR-Clock2"),
  aa = c("InflamAA1", "InflamAA2", "ProtAA1", "ProtAA2"),
  panel = c("Inflammation", "Inflammation", "Proteomics", "Proteomics"),
  clock_type = c("Clock1", "Clock2", "Clock1", "Clock2")
)

standardize_sex <- function(x){
  x_chr <- as.character(x)
  dplyr::case_when(
    x_chr %in% c("Men", "Male", "M", "1") ~ "Men",
    x_chr %in% c("Women", "Female", "F", "0", "2") ~ "Women",
    TRUE ~ x_chr
  )
}

detect_time_var <- function(data){
  candidates <- c("followup_time", "followup_year", "followup")
  hit <- candidates[candidates %in% names(data)]
  if(length(hit) == 0) stop("No follow-up time variable found.")
  hit[1]
}

detect_event_var <- function(data){
  candidates <- c("death_10y", "mortality10y")
  hit <- candidates[candidates %in% names(data)]
  if(length(hit) == 0) stop("No 10-year mortality variable found.")
  hit[1]
}

prepare_mrage_data <- function(data, dataset_name){
  time_var <- detect_time_var(data)
  event_var <- detect_event_var(data)
  data %>%
    mutate(
      dataset = dataset_name,
      followup_time = as.numeric(.data[[time_var]]),
      mortality10y = as.integer(as.numeric(as.character(.data[[event_var]])) == 1),
      sex_standardized = standardize_sex(sex_share)
    )
}

add_age_acceleration <- function(data){
  out <- data
  for(i in seq_len(nrow(clock_map))){
    clock_i <- clock_map$clock[i]
    aa_i <- clock_map$aa[i]
    if(clock_i %in% names(out)){
      out[[aa_i]] <- as.numeric(out[[clock_i]]) - as.numeric(out$age_share)
    }
  }
  out
}

prepare_covariates <- function(data, include_sex = TRUE){
  covars <- c("age_share", "sex_share")
  if(!include_sex) covars <- setdiff(covars, "sex_share")
  covars <- intersect(covars, names(data))
  for(v in covars){
    if(v == "age_share"){
      data[[v]] <- as.numeric(data[[v]])
    } else {
      data[[v]] <- as.factor(data[[v]])
    }
  }
  list(data = droplevels(data), covariates = covars)
}

run_cox_for_aa <- function(data, aa_var, outcome_var, group_name, dataset_name, panel_name, clock_type, include_sex = TRUE){
  prep <- prepare_covariates(data, include_sex = include_sex)
  dat <- prep$data
  covars <- prep$covariates
  analysis_vars <- unique(c("followup_time", outcome_var, aa_var, covars))
  dat <- dat[, analysis_vars, drop = FALSE] %>%
    mutate(
      followup_time = as.numeric(followup_time),
      outcome = as.integer(as.numeric(as.character(.data[[outcome_var]])) == 1),
      aa_value = as.numeric(.data[[aa_var]])
    ) %>%
    filter(!is.na(followup_time), followup_time > 0, !is.na(outcome), !is.na(aa_value)) %>%
    tidyr::drop_na() %>%
    droplevels()
  if(nrow(dat) == 0 || length(unique(dat$outcome)) < 2){
    return(tibble(
      Dataset = dataset_name, Group = group_name, Outcome = outcome_var, Panel = panel_name,
      Clock = clock_type, AA = aa_var, N = nrow(dat), Events = sum(dat$outcome == 1, na.rm = TRUE),
      HR = NA_real_, CI_lower = NA_real_, CI_upper = NA_real_, P_value = NA_real_
    ))
  }
  covars_final <- covars[purrr::map_lgl(covars, function(v){
    if(!v %in% names(dat)) return(FALSE)
    if(is.factor(dat[[v]])) nlevels(dat[[v]]) >= 2 else sd(dat[[v]], na.rm = TRUE) > 0
  })]
  rhs <- paste(c("aa_value", covars_final), collapse = " + ")
  fml <- as.formula(paste0("survival::Surv(followup_time, outcome) ~ ", rhs))
  fit <- survival::coxph(fml, data = dat, x = TRUE, y = TRUE, model = TRUE)
  res <- broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(term == "aa_value") %>%
    transmute(
      Dataset = dataset_name,
      Group = group_name,
      Outcome = outcome_var,
      Panel = panel_name,
      Clock = clock_type,
      AA = aa_var,
      N = nrow(dat),
      Events = sum(dat$outcome == 1, na.rm = TRUE),
      HR = estimate,
      CI_lower = conf.low,
      CI_upper = conf.high,
      P_value = p.value
    )
  res
}

run_all_groups <- function(data, aa_var, outcome_var, dataset_name, panel_name, clock_type){
  overall <- run_cox_for_aa(data, aa_var, outcome_var, "Overall", dataset_name, panel_name, clock_type, include_sex = TRUE)
  men <- run_cox_for_aa(filter(data, sex_standardized == "Men"), aa_var, outcome_var, "Men", dataset_name, panel_name, clock_type, include_sex = FALSE)
  women <- run_cox_for_aa(filter(data, sex_standardized == "Women"), aa_var, outcome_var, "Women", dataset_name, panel_name, clock_type, include_sex = FALSE)
  bind_rows(overall, men, women)
}

format_hr_ci <- function(hr, lower, upper){
  sprintf("%.2f (%.2f–%.2f)", hr, lower, upper)
}

format_results <- function(data){
  data %>%
    mutate(
      HR_95CI = ifelse(is.na(HR), NA_character_, format_hr_ci(HR, CI_lower, CI_upper)),
      P_value = signif(P_value, 3),
      FDR = signif(FDR, 3)
    )
}

ukb_internal <- read_csv(file_ukb_internal_mrage, show_col_types = FALSE) %>%
  prepare_mrage_data("UKB internal validation") %>%
  add_age_acceleration()

esther_ti <- read_csv(file_esther_ti_mrage, show_col_types = FALSE) %>%
  prepare_mrage_data("ESTHER Target Inflammation") %>%
  add_age_acceleration()

esther_ht <- read_csv(file_esther_ht_mrage, show_col_types = FALSE) %>%
  prepare_mrage_data("ESTHER Explore HT") %>%
  add_age_acceleration()

table3_internal <- clock_map %>%
  filter(aa %in% names(ukb_internal)) %>%
  pmap_dfr(function(clock, aa, panel, clock_type){
    run_all_groups(ukb_internal, aa, "mortality10y", "UKB internal validation", panel, clock_type)
  })

table3_esther_ti <- clock_map %>%
  filter(panel == "Inflammation", aa %in% names(esther_ti)) %>%
  pmap_dfr(function(clock, aa, panel, clock_type){
    run_all_groups(esther_ti, aa, "mortality10y", "ESTHER external validation", panel, clock_type)
  })

table3_esther_ht <- clock_map %>%
  filter(panel == "Proteomics", aa %in% names(esther_ht)) %>%
  pmap_dfr(function(clock, aa, panel, clock_type){
    run_all_groups(esther_ht, aa, "mortality10y", "ESTHER external validation", panel, clock_type)
  })

table3 <- bind_rows(table3_internal, table3_esther_ti, table3_esther_ht) %>%
  group_by(Dataset, Outcome, Group) %>%
  mutate(FDR = p.adjust(P_value, method = "BH")) %>%
  ungroup() %>%
  arrange(Dataset, Group, Panel, Clock)

write_csv(table3, file.path(dir_tables, "Table_3_age_acceleration_all_cause_mortality.csv"))
write_csv(format_results(table3), file.path(dir_tables, "Table_3_age_acceleration_all_cause_mortality_formatted.csv"))

ukb_derivation <- read_csv(file_ukb_derivation_mrage, show_col_types = FALSE) %>%
  prepare_mrage_data("UKB derivation")
ukb_internal_for_cause <- read_csv(file_ukb_internal_mrage, show_col_types = FALSE) %>%
  prepare_mrage_data("UKB internal validation")

ukb_all <- bind_rows(ukb_derivation, ukb_internal_for_cause) %>%
  add_age_acceleration()

if(file.exists(file_cause_death)){
  cause_death <- read_csv(file_cause_death, show_col_types = FALSE)
  ukb_all <- ukb_all %>%
    left_join(cause_death %>% select(id, CVD_mort, CA_mort), by = "id")
} else {
  warning("Cause-specific mortality file not found: ", file_cause_death)
  ukb_all$CVD_mort <- NA_real_
  ukb_all$CA_mort <- NA_real_
}

ukb_all <- ukb_all %>%
  mutate(
    CVD_mortality10y = ifelse(mortality10y == 1 & as.numeric(as.character(CVD_mort)) == 1, 1, 0),
    cancer_mortality10y = ifelse(mortality10y == 1 & as.numeric(as.character(CA_mort)) == 1, 1, 0)
  )

cause_summary <- ukb_all %>%
  summarise(
    N = n(),
    all_cause_deaths = sum(mortality10y == 1, na.rm = TRUE),
    CVD_deaths = sum(CVD_mortality10y == 1, na.rm = TRUE),
    cancer_deaths = sum(cancer_mortality10y == 1, na.rm = TRUE)
  )

write_csv(cause_summary, file.path(dir_tables, "Supplementary_Table_11_cause_specific_event_counts.csv"))

supp_table11 <- expand_grid(
  aa = clock_map$aa,
  outcome_var = c("CVD_mortality10y", "cancer_mortality10y")
) %>%
  left_join(clock_map, by = "aa") %>%
  filter(aa %in% names(ukb_all)) %>%
  pmap_dfr(function(aa, outcome_var, clock, panel, clock_type){
    run_all_groups(ukb_all, aa, outcome_var, "UK Biobank", panel, clock_type)
  }) %>%
  group_by(Outcome, Group) %>%
  mutate(FDR = p.adjust(P_value, method = "BH")) %>%
  ungroup() %>%
  arrange(Outcome, Group, Panel, Clock)

write_csv(supp_table11, file.path(dir_tables, "Supplementary_Table_11_age_acceleration_cause_specific_mortality.csv"))
write_csv(format_results(supp_table11), file.path(dir_tables, "Supplementary_Table_11_age_acceleration_cause_specific_mortality_formatted.csv"))

cat("\nAge acceleration mortality association analyses completed.\n")
cat("Generated Table 3 and Supplementary Table 11.\n")
