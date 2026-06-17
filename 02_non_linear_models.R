# ==============================================================================
# 02_non_linear_models.R
# ==============================================================================
# Title:
# Climate Anomalies, Ambient Air Pollution, and Childhood Acute Respiratory
# Infections in Sub-Saharan Africa
#
# Script:
# Nonlinear Bayesian Hierarchical INLA/GAM Models for Prenatal and Postnatal
# Climate Exposures
#
# Author:
# Prince M. Amegbor
#
# Purpose:
# Fit nonlinear Bayesian hierarchical models using R-INLA to estimate smooth
# exposure-response associations between prenatal/postnatal climate anomalies
# and childhood acute respiratory infections (ARI), with progressive adjustment
# for ambient air pollution and child, maternal, and household covariates.
#
# Input:
#   data/DHS_ARI.rds
#
# Outputs:
#   outputs/tables/Prenatal_GAM.xlsx
#   outputs/tables/Postnatal_GAM.xlsx
#   outputs/tables/Nonlinear_Model_Fit_Statistics.xlsx
#   outputs/nonlinear_models/prenatal_gam_results.rds
#   outputs/nonlinear_models/postnatal_gam_results.rds
#   outputs/nonlinear_models/sessionInfo_non_linear_models.txt
#
# Notes:
# - The nonlinear climate smooths use RW2 terms.
# - The survey/sample weight term uses an RW1 prior consistently across all models.
# - Regional and cluster-level random effects are modeled as IID terms.
# - The cloglog link and INLA computational settings are preserved from the
#   original analysis script.
# - Data are not included in this repository. Place the analytic dataset in:
#     data/DHS_ARI.rds
# ==============================================================================


# ==============================================================================
# 0. Setup
# ==============================================================================

set.seed(123)

# Install INLA if needed:
# install.packages(
#   "INLA",
#   repos = c(
#     getOption("repos"),
#     INLA = "https://inla.r-inla-download.org/R/stable"
#   )
# )

library(dplyr)
library(purrr)
library(tibble)
library(INLA)
library(openxlsx)
library(here)

dir.create(here("outputs"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("outputs", "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("outputs", "nonlinear_models"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("outputs", "nonlinear_models", "models"), recursive = TRUE, showWarnings = FALSE)


# ==============================================================================
# 1. Load Data
# ==============================================================================

DF_ARI <- readRDS(
  here("data", "DHS_ARI.rds")
)


# ==============================================================================
# 2. Data Validation
# ==============================================================================

required_vars <- c(
  # Outcome and design variables
  "ch_ari",
  "Surv_year",
  "wt",
  "Reg_ID",
  "admclust",
  
  # Prenatal climate exposures
  "IU_Tmp_z",
  "IU_Tmx_z",
  "IU_Preci_z",
  
  # Postnatal climate exposures
  "PU_Tmp_z",
  "PU_Tmx_z",
  "PU_Preci_z",
  
  # Prenatal air pollution and greenness exposures
  "IU_PM",
  "IU_CO",
  "IU_NO",
  "IU_SO",
  "IU_SU",
  "IU_O3",
  "IU_BC",
  "IU_OC",
  "IU_EVI",
  
  # Postnatal air pollution and greenness exposures
  "PU_PM",
  "PU_CO",
  "PU_NO",
  "PU_SO",
  "PU_SU",
  "PU_O3",
  "PU_BC",
  "PU_OC",
  "PU_EVI",
  
  # Child, maternal, and household covariates
  "b4",
  "b8",
  "ch_allvac_either",
  "rc_edu",
  "rc_empl",
  "mat_age_cat",
  "v190"
)

missing_vars <- setdiff(
  required_vars,
  names(DF_ARI)
)

if (length(missing_vars) > 0) {
  stop(
    "The following required variables are missing from DF_ARI: ",
    paste(missing_vars, collapse = ", ")
  )
}


# ==============================================================================
# 3. Prior Specification
# ==============================================================================

# Prior for precision of the RW1 term used for survey/sample weight smoothing.
prec.prior1 <- list(
  theta = list(
    prior = "pc.prec",
    param = c(0.5, 0.01)
  )
)

# Prior for precision of IID regional and cluster random effects.
prec.prior3 <- list(
  prec = list(
    prior = "pc.prec",
    param = c(1, 0.01)
  )
)


# ==============================================================================
# 4. Variable Naming Notes
# ==============================================================================

# Exposure naming convention:
#
# IU  = prenatal exposure window
# PU  = postnatal exposure window
#
# Tmp   = mean temperature anomaly
# Tmx   = maximum temperature anomaly
# Preci = precipitation anomaly
#
# PM  = particulate matter
# CO  = carbon monoxide
# NO  = nitrogen dioxide
# SO  = sulfur dioxide
# SU  = sulfate aerosol
# O3  = ozone
# BC  = black carbon
# OC  = organic carbon
# EVI = enhanced vegetation index


# ==============================================================================
# 5. Model Components
# ==============================================================================

make_climate_smooth <- function(variable_name) {
  
  paste0(
    "f(inla.group(",
    variable_name,
    "), model = 'rw2', scale.model = TRUE)"
  )
}

make_weight_smooth <- function() {
  
  paste(
    "f(inla.group(wt), model = 'rw1',",
    "scale.model = TRUE, hyper = prec.prior1)"
  )
}

make_random_effects <- function() {
  
  paste(
    make_weight_smooth(),
    "+ f(Reg_ID, model = 'iid', hyper = prec.prior3)",
    "+ f(admclust, model = 'iid', hyper = prec.prior3)"
  )
}

prenatal_pollution_terms_mod4 <- c(
  "IU_PM",
  "IU_CO",
  "IU_NO",
  "IU_SO",
  "IU_SU",
  "IU_O3",
  "IU_BC",
  "IU_OC",
  "IU_EVI"
)

# NOTE:
# The uploaded original nonlinear prenatal Mod5 omitted IU_EVI, although IU_EVI
# was included in prenatal Mod4 and analogous EVI terms were used elsewhere.
# To reproduce the uploaded original script exactly, IU_EVI is not included below.
# If this omission was unintended, add "IU_EVI" to this vector.
prenatal_pollution_terms_mod5 <- c(
  "IU_PM",
  "IU_CO",
  "IU_NO",
  "IU_SO",
  "IU_SU",
  "IU_O3",
  "IU_BC",
  "IU_OC"
)

postnatal_pollution_terms <- c(
  "PU_PM",
  "PU_CO",
  "PU_NO",
  "PU_SO",
  "PU_SU",
  "PU_O3",
  "PU_BC",
  "PU_OC",
  "PU_EVI"
)

socio_demo_terms <- c(
  "factor(b4)",
  "factor(b8)",
  "factor(ch_allvac_either)",
  "factor(rc_edu)",
  "factor(rc_empl)",
  "factor(mat_age_cat)",
  "factor(v190)"
)


# ==============================================================================
# 6. Helper Functions
# ==============================================================================

build_gam_formula <- function(
    smooth_terms,
    linear_terms = NULL,
    include_socio_demo = FALSE
) {
  
  rhs_terms <- smooth_terms
  
  if (!is.null(linear_terms)) {
    rhs_terms <- c(
      rhs_terms,
      linear_terms
    )
  }
  
  if (include_socio_demo) {
    rhs_terms <- c(
      rhs_terms,
      socio_demo_terms
    )
  }
  
  rhs_terms <- c(
    rhs_terms,
    "Surv_year",
    make_random_effects()
  )
  
  formula_text <- paste(
    "ch_ari ~",
    paste(rhs_terms, collapse = " + ")
  )
  
  as.formula(
    formula_text,
    env = parent.frame()
  )
}

run_inla_model <- function(
    formula_object,
    data
) {
  
  inla(
    formula_object,
    data = data,
    family = "binomial",
    Ntrials = 1,
    inla.mode = "experimental",
    
    control.predictor = list(
      compute = TRUE
    ),
    
    control.family = list(
      link = "cloglog"
    ),
    
    control.inla = list(
      strategy = "adaptive",
      int.strategy = "eb"
    ),
    
    num.threads = "8:1",
    
    control.compute = list(
      openmp.strategy = "huge",
      smtp = "default",
      dic = TRUE,
      waic = TRUE
    ),
    
    verbose = FALSE
  )
}

find_smooth_name <- function(
    model,
    variable_name
) {
  
  random_names <- names(model$summary.random)
  
  smooth_name <- random_names[
    grepl(
      variable_name,
      random_names,
      fixed = TRUE
    )
  ][1]
  
  if (is.na(smooth_name)) {
    stop(
      "Could not find smooth term for variable: ",
      variable_name,
      ". Available random-effect terms are: ",
      paste(random_names, collapse = ", ")
    )
  }
  
  smooth_name
}

extract_smooth_result <- function(
    model,
    variable_name
) {
  
  smooth_name <- find_smooth_name(
    model = model,
    variable_name = variable_name
  )
  
  as.data.frame(
    model$summary.random[[smooth_name]]
  ) %>%
    rename(
      lci = `0.025quant`,
      uci = `0.975quant`
    ) %>%
    mutate(
      across(
        c(
          mean,
          lci,
          uci
        ),
        exp
      ),
      DIC = model$dic$dic,
      WAIC = model$waic$waic
    )
}

extract_model_fit <- function(
    model,
    model_name,
    exposure_window
) {
  
  data.frame(
    Exposure_Window = exposure_window,
    Model = model_name,
    DIC = model$dic$dic,
    WAIC = model$waic$waic
  )
}

run_model_set <- function(
    formula_list,
    output_smooth_map,
    exposure_window,
    data,
    save_model_objects = FALSE
) {
  
  model_results <- list()
  smooth_results <- list()
  fit_results <- list()
  
  for (model_name in names(formula_list)) {
    
    cat(
      "Running",
      exposure_window,
      "nonlinear model:",
      model_name,
      "\n"
    )
    
    fitted_model <- run_inla_model(
      formula_object = formula_list[[model_name]],
      data = data
    )
    
    fit_results[[model_name]] <- extract_model_fit(
      model = fitted_model,
      model_name = model_name,
      exposure_window = exposure_window
    )
    
    for (sheet_name in names(output_smooth_map[[model_name]])) {
      
      variable_name <- output_smooth_map[[model_name]][[sheet_name]]
      
      smooth_results[[sheet_name]] <- extract_smooth_result(
        model = fitted_model,
        variable_name = variable_name
      )
    }
    
    if (save_model_objects) {
      model_results[[model_name]] <- fitted_model
    }
    
    if (!save_model_objects) {
      rm(fitted_model)
      gc()
    }
  }
  
  list(
    smooths = smooth_results,
    fit = fit_results,
    models = model_results
  )
}


# ==============================================================================
# 7. Prenatal Nonlinear Model Specifications
# ==============================================================================

prenatal_formulas <- list(
  
  Mod1a = build_gam_formula(
    smooth_terms = c(
      make_climate_smooth("IU_Tmp_z")
    )
  ),
  
  Mod1b = build_gam_formula(
    smooth_terms = c(
      make_climate_smooth("IU_Tmx_z")
    )
  ),
  
  Mod2 = build_gam_formula(
    smooth_terms = c(
      make_climate_smooth("IU_Preci_z")
    )
  ),
  
  Mod3 = build_gam_formula(
    smooth_terms = c(
      make_climate_smooth("IU_Tmp_z"),
      make_climate_smooth("IU_Tmx_z"),
      make_climate_smooth("IU_Preci_z")
    )
  ),
  
  Mod4 = build_gam_formula(
    smooth_terms = c(
      make_climate_smooth("IU_Tmp_z"),
      make_climate_smooth("IU_Tmx_z"),
      make_climate_smooth("IU_Preci_z")
    ),
    linear_terms = prenatal_pollution_terms_mod4
  ),
  
  Mod5 = build_gam_formula(
    smooth_terms = c(
      make_climate_smooth("IU_Tmp_z"),
      make_climate_smooth("IU_Tmx_z"),
      make_climate_smooth("IU_Preci_z")
    ),
    linear_terms = prenatal_pollution_terms_mod5,
    include_socio_demo = TRUE
  )
)

prenatal_output_smooths <- list(
  
  Mod1a = list(
    IN_GAM1_TEMP = "IU_Tmp_z"
  ),
  
  Mod1b = list(
    IN_GAM1_TMX = "IU_Tmx_z"
  ),
  
  Mod2 = list(
    IN_GAM2_RAIN = "IU_Preci_z"
  ),
  
  Mod3 = list(
    IN_GAM3_TEMP = "IU_Tmp_z",
    IN_GAM3_TMX = "IU_Tmx_z",
    IN_GAM3_RAIN = "IU_Preci_z"
  ),
  
  Mod4 = list(
    IN_GAM4_TEMP = "IU_Tmp_z",
    IN_GAM4_TMX = "IU_Tmx_z",
    IN_GAM4_RAIN = "IU_Preci_z"
  ),
  
  Mod5 = list(
    IN_GAM5_TEMP = "IU_Tmp_z",
    IN_GAM5_TMX = "IU_Tmx_z",
    IN_GAM5_RAIN = "IU_Preci_z"
  )
)


# ==============================================================================
# 8. Postnatal Nonlinear Model Specifications
# ==============================================================================

# Correction:
# In the uploaded original script, postnatal Mod1a and Mod1b omitted
# hyper = prec.prior1 for the RW1 weight term. This corrected version applies
# the same RW1 prior consistently across all postnatal models.

postnatal_formulas <- list(
  
  Mod1a = build_gam_formula(
    smooth_terms = c(
      make_climate_smooth("PU_Tmp_z")
    )
  ),
  
  Mod1b = build_gam_formula(
    smooth_terms = c(
      make_climate_smooth("PU_Tmx_z")
    )
  ),
  
  Mod2 = build_gam_formula(
    smooth_terms = c(
      make_climate_smooth("PU_Preci_z")
    )
  ),
  
  Mod3 = build_gam_formula(
    smooth_terms = c(
      make_climate_smooth("PU_Tmp_z"),
      make_climate_smooth("PU_Tmx_z"),
      make_climate_smooth("PU_Preci_z")
    )
  ),
  
  Mod4 = build_gam_formula(
    smooth_terms = c(
      make_climate_smooth("PU_Tmp_z"),
      make_climate_smooth("PU_Tmx_z"),
      make_climate_smooth("PU_Preci_z")
    ),
    linear_terms = postnatal_pollution_terms
  ),
  
  Mod5 = build_gam_formula(
    smooth_terms = c(
      make_climate_smooth("PU_Tmp_z"),
      make_climate_smooth("PU_Tmx_z"),
      make_climate_smooth("PU_Preci_z")
    ),
    linear_terms = postnatal_pollution_terms,
    include_socio_demo = TRUE
  )
)

postnatal_output_smooths <- list(
  
  Mod1a = list(
    PU_GAM1_TEMP = "PU_Tmp_z"
  ),
  
  Mod1b = list(
    PU_GAM1_TMX = "PU_Tmx_z"
  ),
  
  Mod2 = list(
    PU_GAM2_RAIN = "PU_Preci_z"
  ),
  
  Mod3 = list(
    PU_GAM3_TEMP = "PU_Tmp_z",
    PU_GAM3_TMX = "PU_Tmx_z",
    PU_GAM3_RAIN = "PU_Preci_z"
  ),
  
  Mod4 = list(
    PU_GAM4_TEMP = "PU_Tmp_z",
    PU_GAM4_TMX = "PU_Tmx_z",
    PU_GAM4_RAIN = "PU_Preci_z"
  ),
  
  Mod5 = list(
    PU_GAM5_TEMP = "PU_Tmp_z",
    PU_GAM5_TMX = "PU_Tmx_z",
    PU_GAM5_RAIN = "PU_Preci_z"
  )
)


# ==============================================================================
# 9. Run Models
# ==============================================================================

# Full INLA model objects can be large. Keep this FALSE for routine GitHub use.
# Set TRUE if you want to save all fitted model objects for supplementary checks.
save_model_objects <- FALSE

cat(
  "\nStarting prenatal nonlinear models...\n"
)

prenatal_output <- run_model_set(
  formula_list = prenatal_formulas,
  output_smooth_map = prenatal_output_smooths,
  exposure_window = "Prenatal",
  data = DF_ARI,
  save_model_objects = save_model_objects
)

prenatal_gam_results <- prenatal_output$smooths

cat(
  "\nStarting postnatal nonlinear models...\n"
)

postnatal_output <- run_model_set(
  formula_list = postnatal_formulas,
  output_smooth_map = postnatal_output_smooths,
  exposure_window = "Postnatal",
  data = DF_ARI,
  save_model_objects = save_model_objects
)

postnatal_gam_results <- postnatal_output$smooths


# ==============================================================================
# 10. Export Results
# ==============================================================================

write.xlsx(
  prenatal_gam_results,
  file = here(
    "outputs",
    "tables",
    "Prenatal_GAM.xlsx"
  ),
  overwrite = TRUE
)

write.xlsx(
  postnatal_gam_results,
  file = here(
    "outputs",
    "tables",
    "Postnatal_GAM.xlsx"
  ),
  overwrite = TRUE
)

saveRDS(
  prenatal_gam_results,
  file = here(
    "outputs",
    "nonlinear_models",
    "prenatal_gam_results.rds"
  )
)

saveRDS(
  postnatal_gam_results,
  file = here(
    "outputs",
    "nonlinear_models",
    "postnatal_gam_results.rds"
  )
)

fit_statistics <- bind_rows(
  c(
    prenatal_output$fit,
    postnatal_output$fit
  )
)

write.xlsx(
  fit_statistics,
  file = here(
    "outputs",
    "tables",
    "Nonlinear_Model_Fit_Statistics.xlsx"
  ),
  overwrite = TRUE
)

if (save_model_objects) {
  
  saveRDS(
    prenatal_output$models,
    file = here(
      "outputs",
      "nonlinear_models",
      "models",
      "prenatal_gam_models.rds"
    )
  )
  
  saveRDS(
    postnatal_output$models,
    file = here(
      "outputs",
      "nonlinear_models",
      "models",
      "postnatal_gam_models.rds"
    )
  )
}


# ==============================================================================
# 11. Session Information
# ==============================================================================

writeLines(
  capture.output(
    sessionInfo()
  ),
  con = here(
    "outputs",
    "nonlinear_models",
    "sessionInfo_non_linear_models.txt"
  )
)



cat(
  "\nNonlinear model analysis completed successfully.\n"
)

# ==============================================================================
# END
# ==============================================================================