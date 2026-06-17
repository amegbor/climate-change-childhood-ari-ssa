# ==============================================================================
# 03_effect_modification_evi.R
# ==============================================================================
# Title:
# Climate Anomalies, Vegetation Cover, and Childhood Acute Respiratory Infections
#
# Script:
# Effect Modification by Enhanced Vegetation Index (EVI)
#
# Author:
# Prince M. Amegbor
#
# Purpose:
# Fit stratified nonlinear INLA models to assess whether prenatal and postnatal
# vegetation cover, defined using median-split EVI, modifies associations between
# climate anomalies and childhood acute respiratory infections (ARI).
#
# Input:
#   data/DHS_ARI.rds
#
# Outputs:
#   outputs/effect_modification_evi/evi_smooths_unadjusted.rds
#   outputs/effect_modification_evi/evi_smooths_adjusted.rds
#   outputs/effect_modification_evi/EVI_Effect_Modification_Fit_Stats.xlsx
#   outputs/figures/Figure_2.tif
#   outputs/figures/Figure_C1.tif
#   outputs/figures/Figure_C2.tif
#   outputs/figures/Figure_C3.tif
#   outputs/effect_modification_evi/sessionInfo_evi_effect_modification.txt
#
# Notes:
# - IU = prenatal exposure window.
# - PU = postnatal exposure window.
# - EVI categories are created using window-specific medians:
#     1 = Low EVI
#     2 = High EVI
# - Prenatal climate anomalies are paired with prenatal EVI.
# - Postnatal climate anomalies are paired with postnatal EVI.
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
library(ggplot2)
library(INLA)
library(here)
library(cowplot)
library(openxlsx)

dir.create(here("outputs"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("outputs", "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("outputs", "effect_modification_evi"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("outputs", "effect_modification_evi", "models"), recursive = TRUE, showWarnings = FALSE)


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
  # Outcome
  "ch_ari",
  
  # Prenatal climate anomalies
  "IU_Tmp_z",
  "IU_Tmx_z",
  "IU_Preci_z",
  
  # Postnatal climate anomalies
  "PU_Tmp_z",
  "PU_Tmx_z",
  "PU_Preci_z",
  
  # Effect modifier
  "IU_EVI",
  "PU_EVI",
  
  # Adjustment variables
  "b4",
  "b8",
  "ch_allvac_either",
  "rc_edu",
  "rc_empl",
  "mat_age_cat",
  "v190",
  "Surv_year"
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
# 3. Define EVI Categories
# ==============================================================================

iu_evi_median <- median(
  DF_ARI$IU_EVI,
  na.rm = TRUE
)

pu_evi_median <- median(
  DF_ARI$PU_EVI,
  na.rm = TRUE
)

DF_ARI <- DF_ARI %>%
  mutate(
    IU_EVI_cat = ifelse(
      IU_EVI > iu_evi_median,
      "high",
      "low"
    ),
    PU_EVI_cat = ifelse(
      PU_EVI > pu_evi_median,
      "high",
      "low"
    ),
    IU_EVI_cat2 = ifelse(
      IU_EVI_cat == "high",
      2,
      1
    ),
    PU_EVI_cat2 = ifelse(
      PU_EVI_cat == "high",
      2,
      1
    )
  )


# ==============================================================================
# 4. Model Configuration
# ==============================================================================

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

configs <- list(
  TMP_IN = list(
    id = "TMP_IN",
    climate_var = "IU_Tmp_z",
    evi_group = "IU_EVI_cat2",
    window = "Prenatal",
    xlab = "Mean Temperature Anomaly"
  ),
  
  TMP_OUT = list(
    id = "TMP_OUT",
    climate_var = "PU_Tmp_z",
    evi_group = "PU_EVI_cat2",
    window = "Postnatal",
    xlab = "Mean Temperature Anomaly"
  ),
  
  TMX_IN = list(
    id = "TMX_IN",
    climate_var = "IU_Tmx_z",
    evi_group = "IU_EVI_cat2",
    window = "Prenatal",
    xlab = "Maximum Temperature Anomaly"
  ),
  
  TMX_OUT = list(
    id = "TMX_OUT",
    climate_var = "PU_Tmx_z",
    evi_group = "PU_EVI_cat2",
    window = "Postnatal",
    xlab = "Maximum Temperature Anomaly"
  ),
  
  PRE_IN = list(
    id = "PRE_IN",
    climate_var = "IU_Preci_z",
    evi_group = "IU_EVI_cat2",
    window = "Prenatal",
    xlab = "Precipitation Anomaly"
  ),
  
  PRE_OUT = list(
    id = "PRE_OUT",
    climate_var = "PU_Preci_z",
    evi_group = "PU_EVI_cat2",
    window = "Postnatal",
    xlab = "Precipitation Anomaly"
  )
)


# ==============================================================================
# 5. Helper Functions
# ==============================================================================

build_formula <- function(
    climate_var,
    evi_group,
    adjusted = FALSE
) {
  
  smooth_term <- paste0(
    "f(inla.group(",
    climate_var,
    "), model = 'rw2', group = ",
    evi_group,
    ")"
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
    control.family = list(
      link = "cloglog"
    ),
    control.compute = list(
      dic = TRUE,
      waic = TRUE
    ),
    verbose = FALSE
  )
}

find_smooth_name <- function(
    model,
    climate_var
) {
  
  random_names <- names(
    model$summary.random
  )
  
  smooth_name <- random_names[
    grepl(
      climate_var,
      random_names,
      fixed = TRUE
    )
  ][1]
  
  if (is.na(smooth_name)) {
    stop(
      "Could not find smooth term for climate variable: ",
      climate_var,
      ". Available random-effect terms are: ",
      paste(random_names, collapse = ", ")
    )
  }
  
  smooth_name
}

extract_smooth <- function(
    model,
    climate_var
) {
  
  smooth_name <- find_smooth_name(
    model = model,
    climate_var = climate_var
  )
  
  smooth_df <- as.data.frame(
    model$summary.random[[smooth_name]]
  ) %>%
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
      "Check INLA grouped smooth structure before assigning low/high EVI levels."
    )
  }
  
  n_per_group <- nrow(smooth_df) / 2
  
  smooth_df <- smooth_df %>%
    mutate(
      EVI = factor(
        c(
          rep("Low EVI", n_per_group),
          rep("High EVI", n_per_group)
        ),
        levels = c(
          "Low EVI",
          "High EVI"
        )
      ),
      DIC = model$dic$dic,
      WAIC = model$waic$waic
    )
  
  smooth_df
}

extract_fit_stats <- function(
    model,
    model_id,
    adjusted
) {
  
  data.frame(
    Model = model_id,
    Adjusted = adjusted,
    DIC = model$dic$dic,
    WAIC = model$waic$waic
  )
}


# ==============================================================================
# 6. Fit Models and Extract Smooths
# ==============================================================================

smooths_unadj <- list()
smooths_adj <- list()

fit_stats <- list()

save_model_objects <- FALSE
models_unadj <- list()
models_adj <- list()

cat("Fitting EVI effect-modification models...\n")

for (cfg_name in names(configs)) {
  
  cfg <- configs[[cfg_name]]
  
  cat("Running:", cfg$id, "\n")
  
  # Unadjusted model
  formula_unadj <- build_formula(
    climate_var = cfg$climate_var,
    evi_group = cfg$evi_group,
    adjusted = FALSE
  )
  
  model_unadj <- run_inla_model(
    formula_object = formula_unadj,
    data = DF_ARI
  )
  
  smooths_unadj[[cfg$id]] <- extract_smooth(
    model = model_unadj,
    climate_var = cfg$climate_var
  )
  
  fit_stats[[paste0(cfg$id, "_unadj")]] <- extract_fit_stats(
    model = model_unadj,
    model_id = cfg$id,
    adjusted = FALSE
  )
  
  if (save_model_objects) {
    models_unadj[[cfg$id]] <- model_unadj
  }
  
  # Adjusted model
  formula_adj <- build_formula(
    climate_var = cfg$climate_var,
    evi_group = cfg$evi_group,
    adjusted = TRUE
  )
  
  model_adj <- run_inla_model(
    formula_object = formula_adj,
    data = DF_ARI
  )
  
  smooths_adj[[cfg$id]] <- extract_smooth(
    model = model_adj,
    climate_var = cfg$climate_var
  )
  
  fit_stats[[paste0(cfg$id, "_adj")]] <- extract_fit_stats(
    model = model_adj,
    model_id = cfg$id,
    adjusted = TRUE
  )
  
  if (save_model_objects) {
    models_adj[[cfg$id]] <- model_adj
  }
  
  rm(
    model_unadj,
    model_adj,
    formula_unadj,
    formula_adj
  )
  
  gc()
}


# ==============================================================================
# 7. Export Model Outputs
# ==============================================================================

saveRDS(
  smooths_unadj,
  here(
    "outputs",
    "effect_modification_evi",
    "evi_smooths_unadjusted.rds"
  )
)

saveRDS(
  smooths_adj,
  here(
    "outputs",
    "effect_modification_evi",
    "evi_smooths_adjusted.rds"
  )
)

fit_stats_table <- bind_rows(
  fit_stats
)

write.xlsx(
  fit_stats_table,
  here(
    "outputs",
    "effect_modification_evi",
    "EVI_Effect_Modification_Fit_Stats.xlsx"
  ),
  overwrite = TRUE
)

if (save_model_objects) {
  
  saveRDS(
    models_unadj,
    here(
      "outputs",
      "effect_modification_evi",
      "models",
      "evi_models_unadjusted.rds"
    )
  )
  
  saveRDS(
    models_adj,
    here(
      "outputs",
      "effect_modification_evi",
      "models",
      "evi_models_adjusted.rds"
    )
  )
}


# ==============================================================================
# 8. Plot Function
# ==============================================================================

create_plot <- function(
    df,
    title,
    xlab,
    tag = NULL,
    show_legend = FALSE
) {
  
  ggplot(
    df,
    aes(
      x = ID,
      y = mean,
      color = EVI,
      fill = EVI
    )
  ) +
    geom_line(
      linewidth = 1
    ) +
    geom_ribbon(
      aes(
        ymin = lci,
        ymax = uci
      ),
      alpha = 0.20,
      color = NA
    ) +
    labs(
      title = title,
      x = xlab,
      y = "Hazard Ratio",
      tag = tag
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1
      ),
      legend.position = ifelse(
        show_legend,
        "right",
        "none"
      ),
      plot.title = element_text(
        hjust = 0.5
      )
    )
}

make_plot <- function(
    result_list,
    config_id,
    tag = NULL,
    show_legend = FALSE
) {
  
  if (!config_id %in% names(result_list)) {
    stop(
      "Missing result object: ",
      config_id
    )
  }
  
  cfg <- configs[[config_id]]
  
  create_plot(
    df = result_list[[config_id]],
    title = cfg$window,
    xlab = cfg$xlab,
    tag = tag,
    show_legend = show_legend
  )
}

save_figure <- function(
    filename,
    plot_list,
    legend,
    width = 10,
    height = 12
) {
  
  tiff(
    here(
      "outputs",
      "figures",
      filename
    ),
    units = "in",
    width = width,
    height = height,
    res = 300,
    compression = "lzw"
  )
  
  print(
    plot_grid(
      plot_grid(
        plotlist = plot_list,
        ncol = 2
      ),
      legend,
      ncol = 2,
      rel_widths = c(3, 0.4)
    )
  )
  
  dev.off()
}


# ==============================================================================
# 9. Shared Legend
# ==============================================================================

shared_legend <- get_legend(
  make_plot(
    result_list = smooths_adj,
    config_id = "TMP_IN",
    show_legend = TRUE
  ) +
    theme(
      legend.title = element_blank()
    )
)


# ==============================================================================
# 10. Main Figure 2: Adjusted Models Only
# ==============================================================================

figure_2_plots <- list(
  make_plot(
    result_list = smooths_adj,
    config_id = "TMP_IN",
    tag = "A)"
  ),
  make_plot(
    result_list = smooths_adj,
    config_id = "TMP_OUT"
  ),
  make_plot(
    result_list = smooths_adj,
    config_id = "TMX_IN",
    tag = "B)"
  ),
  make_plot(
    result_list = smooths_adj,
    config_id = "TMX_OUT"
  ),
  make_plot(
    result_list = smooths_adj,
    config_id = "PRE_IN",
    tag = "C)"
  ),
  make_plot(
    result_list = smooths_adj,
    config_id = "PRE_OUT"
  )
)

save_figure(
  filename = "Figure_2.tif",
  plot_list = figure_2_plots,
  legend = shared_legend,
  width = 10,
  height = 12
)


# ==============================================================================
# 11. Supplementary Figures C1-C3
# ==============================================================================

build_supplementary_figure <- function(
    id_in,
    id_out,
    filename
) {
  
  supplementary_plots <- list(
    make_plot(
      result_list = smooths_unadj,
      config_id = id_in,
      tag = "A)"
    ),
    make_plot(
      result_list = smooths_unadj,
      config_id = id_out
    ),
    make_plot(
      result_list = smooths_adj,
      config_id = id_in,
      tag = "B)"
    ),
    make_plot(
      result_list = smooths_adj,
      config_id = id_out
    )
  )
  
  save_figure(
    filename = filename,
    plot_list = supplementary_plots,
    legend = shared_legend,
    width = 10,
    height = 8
  )
}

build_supplementary_figure(
  id_in = "TMP_IN",
  id_out = "TMP_OUT",
  filename = "Figure_C1.tif"
)

build_supplementary_figure(
  id_in = "TMX_IN",
  id_out = "TMX_OUT",
  filename = "Figure_C2.tif"
)

build_supplementary_figure(
  id_in = "PRE_IN",
  id_out = "PRE_OUT",
  filename = "Figure_C3.tif"
)


# ==============================================================================
# 12. Session Information
# ==============================================================================

writeLines(
  capture.output(
    sessionInfo()
  ),
  here(
    "outputs",
    "effect_modification_evi",
    "sessionInfo_evi_effect_modification.txt"
  )
)


cat(
  "\nEVI effect-modification models evaluated and figures saved successfully.\n"
)

# ==============================================================================
# END
# ==============================================================================