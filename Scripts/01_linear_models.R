# ==============================================================================
# 01_linear_models.R
# ==============================================================================
# Title:
# Climate Anomalies, Ambient Air Pollution, and Childhood Acute Respiratory
# Infections in Sub-Saharan Africa
#
# Script:
# Linear Bayesian Hierarchical INLA Models for Prenatal and Postnatal Exposures
#
# Author:
# Prince M. Amegbor
#
# Purpose:
# Fit Bayesian hierarchical binomial models using R-INLA to estimate associations
# between prenatal/postnatal climate anomalies, ambient air pollution, and
# childhood acute respiratory infections (ARI).
#
# Input:
#   data/DHS_ARI.rds
#
# Outputs:
#   outputs/tables/Prenatal_ARI.xlsx
#   outputs/tables/Postnatal_ARI.xlsx
#   outputs/linear_models/sessionInfo_linear_models.txt
#
# Notes:
# - The model specifications, priors, covariates, random effects, cloglog link,
#   DIC/WAIC computation, and INLA computational settings are preserved from
#   the original analysis script.
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
dir.create(here("outputs", "linear_models"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("outputs", "linear_models", "models"), recursive = TRUE, showWarnings = FALSE)


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

# Prior for precision of IID random effects.
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

random_effects <- paste(
  "f(inla.group(wt), model = 'rw1', scale.model = TRUE, hyper = prec.prior1)",
  "+ f(Reg_ID, model = 'iid', hyper = prec.prior3)",
  "+ f(admclust, model = 'iid', hyper = prec.prior3)"
)

prenatal_pollution_terms <- c(
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

build_formula <- function(
    exposure_terms,
    pollution_terms = NULL,
    include_socio_demo = FALSE
) {
  
  rhs_terms <- exposure_terms
  
  if (!is.null(pollution_terms)) {
    rhs_terms <- c(
      rhs_terms,
      pollution_terms
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
    "Surv_year"
  )
  
  rhs <- paste(
    rhs_terms,
    collapse = " + "
  )
  
  formula_text <- paste(
    "ch_ari ~",
    rhs,
    "+",
    random_effects
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

extract_linear_results <- function(
    model
) {
  
  output <- as.data.frame(
    model$summary.fixed[
      ,
      c(
        "mean",
        "sd",
        "0.025quant",
        "0.975quant"
      )
    ]
  ) %>%
    rownames_to_column(
      "variables"
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
      )
    )
  
  output$DIC <- model$dic$dic
  output$WAIC <- model$waic$waic
  
  output
}

run_model_set <- function(
    formula_list,
    data,
    save_model_objects = FALSE
) {
  
  result_list <- list()
  model_list <- list()
  
  for (model_name in names(formula_list)) {
    
    cat(
      "Running model:",
      model_name,
      "\n"
    )
    
    fitted_model <- run_inla_model(
      formula_object = formula_list[[model_name]],
      data = data
    )
    
    result_list[[model_name]] <- extract_linear_results(
      model = fitted_model
    )
    
    if (save_model_objects) {
      model_list[[model_name]] <- fitted_model
    }
    
    if (!save_model_objects) {
      rm(fitted_model)
      gc()
    }
  }
  
  list(
    results = result_list,
    models = model_list
  )
}


# ==============================================================================
# 7. Prenatal Model Specifications
# ==============================================================================

prenatal_formulas <- list(
  
  Mod1a = build_formula(
    exposure_terms = c(
      "IU_Tmp_z"
    )
  ),
  
  Mod1b = build_formula(
    exposure_terms = c(
      "IU_Tmx_z"
    )
  ),
  
  Mod2 = build_formula(
    exposure_terms = c(
      "IU_Preci_z"
    )
  ),
  
  Mod3 = build_formula(
    exposure_terms = c(
      "IU_Tmp_z",
      "IU_Tmx_z",
      "IU_Preci_z"
    )
  ),
  
  Mod4 = build_formula(
    exposure_terms = c(
      "IU_Tmp_z",
      "IU_Tmx_z",
      "IU_Preci_z"
    ),
    pollution_terms = prenatal_pollution_terms
  ),
  
  Mod5 = build_formula(
    exposure_terms = c(
      "IU_Tmp_z",
      "IU_Tmx_z",
      "IU_Preci_z"
    ),
    pollution_terms = prenatal_pollution_terms,
    include_socio_demo = TRUE
  )
)


# ==============================================================================
# 8. Postnatal Model Specifications
# ==============================================================================

postnatal_formulas <- list(
  
  Mod1a = build_formula(
    exposure_terms = c(
      "PU_Tmp_z"
    )
  ),
  
  Mod1b = build_formula(
    exposure_terms = c(
      "PU_Tmx_z"
    )
  ),
  
  Mod2 = build_formula(
    exposure_terms = c(
      "PU_Preci_z"
    )
  ),
  
  Mod3 = build_formula(
    exposure_terms = c(
      "PU_Tmp_z",
      "PU_Tmx_z",
      "PU_Preci_z"
    )
  ),
  
  Mod4 = build_formula(
    exposure_terms = c(
      "PU_Tmp_z",
      "PU_Tmx_z",
      "PU_Preci_z"
    ),
    pollution_terms = postnatal_pollution_terms
  ),
  
  Mod5 = build_formula(
    exposure_terms = c(
      "PU_Tmp_z",
      "PU_Tmx_z",
      "PU_Preci_z"
    ),
    pollution_terms = postnatal_pollution_terms,
    include_socio_demo = TRUE
  )
)


# ==============================================================================
# 9. Run Models
# ==============================================================================

# Full INLA model objects can be large. Keep this FALSE for routine GitHub use.
# Set TRUE if you want to save all fitted model objects for supplementary checks.
save_model_objects <- FALSE

cat(
  "\nStarting prenatal linear models...\n"
)

prenatal_output <- run_model_set(
  formula_list = prenatal_formulas,
  data = DF_ARI,
  save_model_objects = save_model_objects
)

prenatal_results <- prenatal_output$results

cat(
  "\nStarting postnatal linear models...\n"
)

postnatal_output <- run_model_set(
  formula_list = postnatal_formulas,
  data = DF_ARI,
  save_model_objects = save_model_objects
)

postnatal_results <- postnatal_output$results


# ==============================================================================
# 10. Export Results
# ==============================================================================

write.xlsx(
  prenatal_results,
  file = here(
    "outputs",
    "tables",
    "Prenatal_ARI.xlsx"
  ),
  overwrite = TRUE
)

write.xlsx(
  postnatal_results,
  file = here(
    "outputs",
    "tables",
    "Postnatal_ARI.xlsx"
  ),
  overwrite = TRUE
)

saveRDS(
  prenatal_results,
  file = here(
    "outputs",
    "linear_models",
    "prenatal_linear_results.rds"
  )
)

saveRDS(
  postnatal_results,
  file = here(
    "outputs",
    "linear_models",
    "postnatal_linear_results.rds"
  )
)

if (save_model_objects) {
  
  saveRDS(
    prenatal_output$models,
    file = here(
      "outputs",
      "linear_models",
      "models",
      "prenatal_linear_models.rds"
    )
  )
  
  saveRDS(
    postnatal_output$models,
    file = here(
      "outputs",
      "linear_models",
      "models",
      "postnatal_linear_models.rds"
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
    "linear_models",
    "sessionInfo_linear_models.txt"
  )
)


# ==============================================================================
# 12. Optional
# ==============================================================================

if (FALSE) {
  
  # ---------------------------------------------------------------------------
  # Optional: export a compact model comparison table
  # ---------------------------------------------------------------------------
  
  extract_model_fit_table <- function(result_list, exposure_window) {
    
    bind_rows(
      lapply(
        names(result_list),
        function(model_name) {
          
          data.frame(
            Exposure_Window = exposure_window,
            Model = model_name,
            DIC = unique(result_list[[model_name]]$DIC),
            WAIC = unique(result_list[[model_name]]$WAIC)
          )
        }
      )
    )
  }
  
  prenatal_fit_table <- extract_model_fit_table(
    prenatal_results,
    exposure_window = "Prenatal"
  )
  
  postnatal_fit_table <- extract_model_fit_table(
    postnatal_results,
    exposure_window = "Postnatal"
  )
  
  linear_fit_table <- bind_rows(
    prenatal_fit_table,
    postnatal_fit_table
  )
  
  write.xlsx(
    linear_fit_table,
    file = here(
      "outputs",
      "tables",
      "Linear_Model_Fit_Statistics.xlsx"
    ),
    overwrite = TRUE
  )
  
}

cat(
  "\nLinear model analysis completed successfully.\n"
)

# ==============================================================================
# END
# ==============================================================================