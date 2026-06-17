# ==============================================================================
# 04_effect_modification_models.R
# ==============================================================================
# Title: Climate, Air Pollution, Effect Modification and Child ARI
# Purpose:
#   Fit stratified nonlinear INLA models assessing whether ambient air pollutants
#   modify climate anomaly–ARI associations.
#
# Outputs:
#   outputs/effect_modification/results_unadj.rds
#   outputs/effect_modification/results_adj.rds
#   outputs/effect_modification/fit_stats.xlsx
#   outputs/effect_modification/sessionInfo_models.txt
# ==============================================================================

# 0. Setup ---------------------------------------------------------------------

set.seed(123)

library(dplyr)
library(tibble)
library(purrr)
library(INLA)
library(here)
library(openxlsx)

dir.create(here("outputs"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("outputs", "effect_modification"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("outputs", "effect_modification", "models"), recursive = TRUE, showWarnings = FALSE)

# 1. Load Data -----------------------------------------------------------------

DF_ARI <- readRDS(here("data", "DHS_ARI.rds"))

# 2. Data Validation -----------------------------------------------------------

required_vars <- c(
  "ch_ari",
  "IU_Tmp_z", "PU_Tmp_z",
  "IU_Tmx_z", "PU_Tmx_z",
  "IU_Preci_z", "PU_Preci_z",
  "IU_PM", "PU_PM",
  "IU_NO", "PU_NO",
  "IU_CO", "PU_CO",
  "IU_SO", "PU_SO",
  "IU_O3", "PU_O3",
  "b4", "b8", "ch_allvac_either",
  "rc_edu", "rc_empl", "mat_age_cat",
  "v190", "Surv_year"
)

missing_vars <- setdiff(required_vars, names(DF_ARI))

if (length(missing_vars) > 0) {
  stop(
    "The following required variables are missing from DF_ARI: ",
    paste(missing_vars, collapse = ", ")
  )
}

# 3. Prior Specification -------------------------------------------------------

# Penalized complexity prior for RW2 smooth precision.
# This mirrors the prior structure used in the manuscript analyses.
rw2_prior <- list(
  prec = list(
    prior = "pc.prec",
    param = c(0.5, 0.01)
  )
)

# 4. Pollutant Categorization --------------------------------------------------

pollutants <- c("PM", "NO", "CO", "SO", "O3")

for (pollutant in pollutants) {
  
  iu_var <- paste0("IU_", pollutant)
  pu_var <- paste0("PU_", pollutant)
  
  iu_median <- median(DF_ARI[[iu_var]], na.rm = TRUE)
  pu_median <- median(DF_ARI[[pu_var]], na.rm = TRUE)
  
  DF_ARI[[paste0(iu_var, "_cat2")]] <- ifelse(
    DF_ARI[[iu_var]] > iu_median,
    2,
    1
  )
  
  DF_ARI[[paste0(pu_var, "_cat2")]] <- ifelse(
    DF_ARI[[pu_var]] > pu_median,
    2,
    1
  )
}

# 5. Climate Variable Definitions ---------------------------------------------

climates <- list(
  TMP = list(
    prefix = "TMP",
    iu = "IU_Tmp_z",
    pu = "PU_Tmp_z",
    label = "Mean Temperature Anomaly"
  ),
  TMX = list(
    prefix = "TMX",
    iu = "IU_Tmx_z",
    pu = "PU_Tmx_z",
    label = "Maximum Temperature Anomaly"
  ),
  PRE = list(
    prefix = "PRE",
    iu = "IU_Preci_z",
    pu = "PU_Preci_z",
    label = "Precipitation Anomaly"
  )
)

# 6. Covariates ----------------------------------------------------------------

covariates <- paste(
  "factor(b4)",
  "factor(b8)",
  "factor(ch_allvac_either)",
  "factor(rc_edu)",
  "factor(rc_empl)",
  "factor(mat_age_cat)",
  "factor(v190)",
  "Surv_year",
  sep = " + "
)

# 7. Helper Functions ----------------------------------------------------------

build_formula <- function(climate_var, group_var, adjusted = FALSE) {
  
  smooth_term <- paste0(
    "f(inla.group(",
    climate_var,
    "), model = 'rw2', group = ",
    group_var,
    ", hyper = rw2_prior)"
  )
  
  formula_text <- paste(
    "ch_ari ~ 1 +",
    smooth_term
  )
  
  if (adjusted) {
    formula_text <- paste(
      formula_text,
      "+",
      covariates
    )
  }
  
  as.formula(formula_text, env = parent.frame())
}

run_inla_model <- function(formula_object, data) {
  
  inla(
    formula_object,
    data = data,
    family = "binomial",
    Ntrials = 1,
    control.family = list(link = "cloglog"),
    control.compute = list(
      dic = TRUE,
      waic = TRUE
    ),
    verbose = FALSE
  )
}

extract_smooth <- function(model, climate_var, pollutant) {
  
  random_names <- names(model$summary.random)
  
  smooth_name <- random_names[
    grepl(climate_var, random_names, fixed = TRUE)
  ][1]
  
  if (is.na(smooth_name)) {
    stop(
      "Could not identify the smooth term for climate variable: ",
      climate_var
    )
  }
  
  smooth_df <- as.data.frame(model$summary.random[[smooth_name]])
  
  smooth_df <- smooth_df %>%
    rename(
      lci = `0.025quant`,
      uci = `0.975quant`
    ) %>%
    mutate(
      mean = exp(mean),
      lci = exp(lci),
      uci = exp(uci)
    )
  
  if (nrow(smooth_df) %% 2 != 0) {
    stop(
      "The extracted smooth for ",
      climate_var,
      " does not have an even number of rows. ",
      "Check INLA group structure before assigning low/high pollutant levels."
    )
  }
  
  n_per_group <- nrow(smooth_df) / 2
  
  smooth_df <- smooth_df %>%
    mutate(
      Pollutant = pollutant,
      Pollutant_Level = factor(
        c(
          rep(paste("Low", pollutant), n_per_group),
          rep(paste("High", pollutant), n_per_group)
        ),
        levels = c(
          paste("Low", pollutant),
          paste("High", pollutant)
        )
      )
    )
  
  smooth_df
}

extract_fit_stats <- function(model, model_name, adjusted) {
  
  data.frame(
    Model = model_name,
    Adjusted = adjusted,
    DIC = model$dic$dic,
    WAIC = model$waic$waic
  )
}

# 8. Run Models ----------------------------------------------------------------

results_unadj <- list()
results_adj <- list()
fit_stats <- list()

# Full model objects are not saved by default to reduce repository/output size.
# Set save_model_objects <- TRUE if you want to archive all fitted INLA objects.
save_model_objects <- FALSE

models_unadj <- list()
models_adj <- list()

cat("Running effect-modification models...\n")

for (pollutant in pollutants) {
  
  cat("\nPollutant:", pollutant, "\n")
  
  for (clim_name in names(climates)) {
    
    clim <- climates[[clim_name]]
    
    cat("  Climate:", clim$prefix, "\n")
    
    # Prenatal pollutant modifier with prenatal climate exposure
    iu_group <- paste0("IU_", pollutant, "_cat2")
    iu_id <- paste0(clim$prefix, "_IU_", pollutant)
    
    formula_iu_unadj <- build_formula(
      climate_var = clim$iu,
      group_var = iu_group,
      adjusted = FALSE
    )
    
    formula_iu_adj <- build_formula(
      climate_var = clim$iu,
      group_var = iu_group,
      adjusted = TRUE
    )
    
    model_iu_unadj <- run_inla_model(
      formula_object = formula_iu_unadj,
      data = DF_ARI
    )
    
    model_iu_adj <- run_inla_model(
      formula_object = formula_iu_adj,
      data = DF_ARI
    )
    
    results_unadj[[iu_id]] <- extract_smooth(
      model = model_iu_unadj,
      climate_var = clim$iu,
      pollutant = pollutant
    )
    
    results_adj[[iu_id]] <- extract_smooth(
      model = model_iu_adj,
      climate_var = clim$iu,
      pollutant = pollutant
    )
    
    fit_stats[[paste0(iu_id, "_unadj")]] <- extract_fit_stats(
      model = model_iu_unadj,
      model_name = iu_id,
      adjusted = FALSE
    )
    
    fit_stats[[paste0(iu_id, "_adj")]] <- extract_fit_stats(
      model = model_iu_adj,
      model_name = iu_id,
      adjusted = TRUE
    )
    
    if (save_model_objects) {
      models_unadj[[iu_id]] <- model_iu_unadj
      models_adj[[iu_id]] <- model_iu_adj
    }
    
    # Postnatal pollutant modifier with postnatal climate exposure
    pu_group <- paste0("PU_", pollutant, "_cat2")
    pu_id <- paste0(clim$prefix, "_PU_", pollutant)
    
    formula_pu_unadj <- build_formula(
      climate_var = clim$pu,
      group_var = pu_group,
      adjusted = FALSE
    )
    
    formula_pu_adj <- build_formula(
      climate_var = clim$pu,
      group_var = pu_group,
      adjusted = TRUE
    )
    
    model_pu_unadj <- run_inla_model(
      formula_object = formula_pu_unadj,
      data = DF_ARI
    )
    
    model_pu_adj <- run_inla_model(
      formula_object = formula_pu_adj,
      data = DF_ARI
    )
    
    results_unadj[[pu_id]] <- extract_smooth(
      model = model_pu_unadj,
      climate_var = clim$pu,
      pollutant = pollutant
    )
    
    results_adj[[pu_id]] <- extract_smooth(
      model = model_pu_adj,
      climate_var = clim$pu,
      pollutant = pollutant
    )
    
    fit_stats[[paste0(pu_id, "_unadj")]] <- extract_fit_stats(
      model = model_pu_unadj,
      model_name = pu_id,
      adjusted = FALSE
    )
    
    fit_stats[[paste0(pu_id, "_adj")]] <- extract_fit_stats(
      model = model_pu_adj,
      model_name = pu_id,
      adjusted = TRUE
    )
    
    if (save_model_objects) {
      models_unadj[[pu_id]] <- model_pu_unadj
      models_adj[[pu_id]] <- model_pu_adj
    }
    
    rm(
      model_iu_unadj,
      model_iu_adj,
      model_pu_unadj,
      model_pu_adj,
      formula_iu_unadj,
      formula_iu_adj,
      formula_pu_unadj,
      formula_pu_adj
    )
    
    gc()
  }
}

# 9. Export Results ------------------------------------------------------------

saveRDS(
  results_unadj,
  here("outputs", "effect_modification", "results_unadj.rds")
)

saveRDS(
  results_adj,
  here("outputs", "effect_modification", "results_adj.rds")
)

fit_stats_table <- bind_rows(fit_stats)

write.xlsx(
  fit_stats_table,
  here("outputs", "effect_modification", "Pollutant_Effect_Modification_Fit_Stats.xlsx"),
  overwrite = TRUE
)

if (save_model_objects) {
  
  saveRDS(
    models_unadj,
    here("outputs", "effect_modification", "models", "models_unadj.rds")
  )
  
  saveRDS(
    models_adj,
    here("outputs", "effect_modification", "models", "models_adj.rds")
  )
}

writeLines(
  capture.output(sessionInfo()),
  here("outputs", "effect_modification", "sessionInfo_models.txt")
)

if (FALSE) {
  
  # Optional diagnostic check:
  # Run this after fitting one model if you want to verify INLA's grouped smooth
  # ordering before using first-half/second-half pollutant labels.
  
  # example_model <- readRDS(
  #   here("outputs", "effect_modification", "models", "models_adj.rds")
  # )[[1]]
  #
  # names(example_model$summary.random)
  # head(example_model$summary.random[[1]])
  # tail(example_model$summary.random[[1]])
}

cat("\nEffect-modification models completed successfully.\n")
# ==============================================================================
# END
# ==============================================================================