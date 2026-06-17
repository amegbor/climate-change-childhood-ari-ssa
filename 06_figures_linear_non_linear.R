# ==============================================================================
# 06_figures_linear_non_linear.R
# ==============================================================================
# Title:
# Climate Anomalies, Ambient Air Pollution, and Childhood Acute Respiratory
# Infections in Sub-Saharan Africa
#
# Script:
# Figures for Linear and Nonlinear Main Models
#
# Author:
# Prince M. Amegbor
#
# Purpose:
# Generate manuscript and supplementary figures from the outputs of:
#
#   01_linear_models.R
#   02_non_linear_models.R
#
# Inputs:
#   outputs/tables/Prenatal_ARI.xlsx
#   outputs/tables/Postnatal_ARI.xlsx
#   outputs/tables/Prenatal_GAM.xlsx
#   outputs/tables/Postnatal_GAM.xlsx
#
# Optional Input:
#   outputs/tables/Visualize_Fixed_Effect.xlsx
#
# Outputs:
#   outputs/figures/Pre_Post_Fixed_effect.tiff
#   outputs/figures/Pre_Post_Fixed_effect.jpeg
#   outputs/figures/Figure_B1.tiff
#   outputs/figures/Figure_B2.tiff
#   outputs/figures/Figure_B3.tiff
#   outputs/figures/Figure_B4.tiff
#   outputs/figures/Figure_B5.tiff
#   outputs/figures/Figure_B6.tiff
#   outputs/figures/sessionInfo_figures_linear_non_linear.txt
#
# Notes:
# - The fixed-effect coefficient plot uses exponentiated estimates from the
#   linear model outputs.
# - The nonlinear plots use exponentiated smooth estimates from the GAM outputs.
# - The script will use Visualize_Fixed_Effect.xlsx if available. Otherwise, it
#   builds the plotting data directly from Prenatal_ARI.xlsx and Postnatal_ARI.xlsx.
# ==============================================================================


# ==============================================================================
# 0. Setup
# ==============================================================================

library(dplyr)
library(purrr)
library(readxl)
library(ggplot2)
library(ggpubr)
library(openxlsx)
library(here)

theme_set(theme_pubr())

dir.create(here("outputs"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("outputs", "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("outputs", "tables"), recursive = TRUE, showWarnings = FALSE)


# ==============================================================================
# 1. Helper Functions
# ==============================================================================

read_excel_workbook <- function(file_path) {
  
  if (!file.exists(file_path)) {
    stop("File not found: ", file_path)
  }
  
  sheet_names <- readxl::excel_sheets(file_path)
  
  sheets <- lapply(
    sheet_names,
    function(sheet_name) {
      readxl::read_excel(
        path = file_path,
        sheet = sheet_name
      ) %>%
        as.data.frame()
    }
  )
  
  names(sheets) <- sheet_names
  
  sheets
}

standardize_result_columns <- function(df) {
  
  names(df) <- gsub(
    "^.*\\.{2}",
    "",
    names(df)
  )
  
  if ("0.025quant" %in% names(df)) {
    df <- df %>%
      rename(
        lci = `0.025quant`
      )
  }
  
  if ("0.975quant" %in% names(df)) {
    df <- df %>%
      rename(
        uci = `0.975quant`
      )
  }
  
  df
}

format_estimate <- function(mean, lci, uci, digits = 2) {
  
  sprintf(
    paste0(
      "%.",
      digits,
      "f (%.",
      digits,
      "f, %.",
      digits,
      "f)"
    ),
    mean,
    lci,
    uci
  )
}

save_tiff <- function(filename, plot, width, height, res = 300) {
  
  tiff(
    filename = here(
      "outputs",
      "figures",
      filename
    ),
    units = "in",
    width = width,
    height = height,
    res = res,
    compression = "lzw"
  )
  
  print(plot)
  
  dev.off()
}

save_jpeg <- function(filename, plot, width, height, res = 300) {
  
  jpeg(
    filename = here(
      "outputs",
      "figures",
      filename
    ),
    units = "in",
    width = width,
    height = height,
    res = res
  )
  
  print(plot)
  
  dev.off()
}


# ==============================================================================
# 2. Linear Fixed-Effect Plotting Data
# ==============================================================================

build_fixed_effect_plot_data <- function(
    linear_results_file,
    exposure_window = c("Prenatal", "Postnatal")
) {
  
  exposure_window <- match.arg(exposure_window)
  
  model_sheets <- read_excel_workbook(
    linear_results_file
  )
  
  model_sheets <- lapply(
    model_sheets,
    standardize_result_columns
  )
  
  if (exposure_window == "Prenatal") {
    
    variable_map <- data.frame(
      Models = c(
        "Model 1", "Model 1", "Model 1",
        "Model 2", "Model 2", "Model 2",
        "Model 3", "Model 3", "Model 3",
        "Model 4", "Model 4", "Model 4"
      ),
      Exposure = c(
        "Mean Temperature anomaly",
        "Maximum Temperature anomaly",
        "Precipitation anomaly",
        "Mean Temperature anomaly",
        "Maximum Temperature anomaly",
        "Precipitation anomaly",
        "Mean Temperature anomaly",
        "Maximum Temperature anomaly",
        "Precipitation anomaly",
        "Mean Temperature anomaly",
        "Maximum Temperature anomaly",
        "Precipitation anomaly"
      ),
      Sheet = c(
        "Mod1a", "Mod1b", "Mod2",
        "Mod3", "Mod3", "Mod3",
        "Mod4", "Mod4", "Mod4",
        "Mod5", "Mod5", "Mod5"
      ),
      Variable = c(
        "IU_Tmp_z",
        "IU_Tmx_z",
        "IU_Preci_z",
        "IU_Tmp_z",
        "IU_Tmx_z",
        "IU_Preci_z",
        "IU_Tmp_z",
        "IU_Tmx_z",
        "IU_Preci_z",
        "IU_Tmp_z",
        "IU_Tmx_z",
        "IU_Preci_z"
      ),
      stringsAsFactors = FALSE
    )
    
  } else {
    
    variable_map <- data.frame(
      Models = c(
        "Model 1", "Model 1", "Model 1",
        "Model 2", "Model 2", "Model 2",
        "Model 3", "Model 3", "Model 3",
        "Model 4", "Model 4", "Model 4"
      ),
      Exposure = c(
        "Mean Temperature anomaly",
        "Maximum Temperature anomaly",
        "Precipitation anomaly",
        "Mean Temperature anomaly",
        "Maximum Temperature anomaly",
        "Precipitation anomaly",
        "Mean Temperature anomaly",
        "Maximum Temperature anomaly",
        "Precipitation anomaly",
        "Mean Temperature anomaly",
        "Maximum Temperature anomaly",
        "Precipitation anomaly"
      ),
      Sheet = c(
        "Mod1a", "Mod1b", "Mod2",
        "Mod3", "Mod3", "Mod3",
        "Mod4", "Mod4", "Mod4",
        "Mod5", "Mod5", "Mod5"
      ),
      Variable = c(
        "PU_Tmp_z",
        "PU_Tmx_z",
        "PU_Preci_z",
        "PU_Tmp_z",
        "PU_Tmx_z",
        "PU_Preci_z",
        "PU_Tmp_z",
        "PU_Tmx_z",
        "PU_Preci_z",
        "PU_Tmp_z",
        "PU_Tmx_z",
        "PU_Preci_z"
      ),
      stringsAsFactors = FALSE
    )
  }
  
  plot_data <- lapply(
    seq_len(nrow(variable_map)),
    function(i) {
      
      sheet_name <- variable_map$Sheet[i]
      variable_name <- variable_map$Variable[i]
      
      if (!sheet_name %in% names(model_sheets)) {
        stop(
          "Expected sheet not found in ",
          linear_results_file,
          ": ",
          sheet_name
        )
      }
      
      model_df <- model_sheets[[sheet_name]]
      
      if (!"variables" %in% names(model_df)) {
        stop(
          "The sheet ",
          sheet_name,
          " must contain a column named 'variables'."
        )
      }
      
      variable_row <- model_df %>%
        filter(
          variables == variable_name
        )
      
      if (nrow(variable_row) != 1) {
        stop(
          "Expected exactly one row for variable ",
          variable_name,
          " in sheet ",
          sheet_name,
          ". Found ",
          nrow(variable_row),
          "."
        )
      }
      
      data.frame(
        Models = variable_map$Models[i],
        Exposure = variable_map$Exposure[i],
        variables = variable_name,
        mean = variable_row$mean,
        lci = variable_row$lci,
        uci = variable_row$uci,
        DIC = variable_row$DIC,
        WAIC = variable_row$WAIC,
        stringsAsFactors = FALSE
      )
    }
  ) %>%
    bind_rows() %>%
    mutate(
      estimates = format_estimate(
        mean = mean,
        lci = lci,
        uci = uci
      )
    )
  
  plot_data
}

load_or_build_fixed_effect_data <- function() {
  
  curated_file <- here(
    "outputs",
    "tables",
    "Visualize_Fixed_Effect.xlsx"
  )
  
  prenatal_linear_file <- here(
    "outputs",
    "tables",
    "Prenatal_ARI.xlsx"
  )
  
  postnatal_linear_file <- here(
    "outputs",
    "tables",
    "Postnatal_ARI.xlsx"
  )
  
  if (file.exists(curated_file)) {
    
    fixed_sheets <- read_excel_workbook(
      curated_file
    )
    
    if (!all(c("Pre_natal_ARI", "Post_natal_ARI") %in% names(fixed_sheets))) {
      stop(
        "Visualize_Fixed_Effect.xlsx must contain sheets named ",
        "'Pre_natal_ARI' and 'Post_natal_ARI'."
      )
    }
    
    Pre_natal_ARI <- fixed_sheets[["Pre_natal_ARI"]]
    Post_natal_ARI <- fixed_sheets[["Post_natal_ARI"]]
    
  } else {
    
    Pre_natal_ARI <- build_fixed_effect_plot_data(
      linear_results_file = prenatal_linear_file,
      exposure_window = "Prenatal"
    )
    
    Post_natal_ARI <- build_fixed_effect_plot_data(
      linear_results_file = postnatal_linear_file,
      exposure_window = "Postnatal"
    )
    
    write.xlsx(
      list(
        Pre_natal_ARI = Pre_natal_ARI,
        Post_natal_ARI = Post_natal_ARI
      ),
      file = curated_file,
      overwrite = TRUE
    )
  }
  
  list(
    Pre_natal_ARI = Pre_natal_ARI,
    Post_natal_ARI = Post_natal_ARI
  )
}


# ==============================================================================
# 3. Linear Fixed-Effect Figure
# ==============================================================================

fixed_effect_data <- load_or_build_fixed_effect_data()

Pre_natal_ARI <- fixed_effect_data$Pre_natal_ARI
Post_natal_ARI <- fixed_effect_data$Post_natal_ARI

Pre_natal_ARI <- Pre_natal_ARI %>%
  mutate(
    Models = factor(
      Models,
      levels = c(
        "Model 4",
        "Model 3",
        "Model 2",
        "Model 1"
      )
    ),
    Exposure = factor(
      Exposure,
      levels = c(
        "Precipitation anomaly",
        "Maximum Temperature anomaly",
        "Mean Temperature anomaly"
      )
    )
  )

Post_natal_ARI <- Post_natal_ARI %>%
  mutate(
    Models = factor(
      Models,
      levels = c(
        "Model 4",
        "Model 3",
        "Model 2",
        "Model 1"
      )
    ),
    Exposure = factor(
      Exposure,
      levels = c(
        "Precipitation anomaly",
        "Maximum Temperature anomaly",
        "Mean Temperature anomaly"
      )
    )
  )

make_fixed_effect_plot <- function(
    df,
    ylab,
    panel_title
) {
  
  dodge <- position_dodge(
    width = 0.75
  )
  
  ggplot(
    df,
    aes(
      x = mean,
      y = Exposure,
      colour = Models
    )
  ) +
    geom_vline(
      xintercept = 1,
      linewidth = 0.25,
      linetype = "dashed"
    ) +
    geom_point(
      position = dodge,
      size = 2.25
    ) +
    geom_errorbarh(
      aes(
        xmin = lci,
        xmax = uci
      ),
      position = dodge,
      height = 0
    ) +
    scale_x_log10() +
    scale_color_brewer(
      palette = "Set1",
      breaks = c(
        "Model 1",
        "Model 2",
        "Model 3",
        "Model 4"
      )
    ) +
    labs(
      x = "Hazard Ratio",
      y = ylab,
      colour = "Models",
      title = panel_title
    ) +
    theme_classic() +
    theme(
      axis.text.y.left = element_text(
        margin = margin(
          r = 5,
          l = 5
        )
      ),
      plot.title = element_text(
        hjust = 0
      )
    )
}

Pre_FE <- make_fixed_effect_plot(
  df = Pre_natal_ARI,
  ylab = "Prenatal Climate Anomalies",
  panel_title = "A"
)

Post_FE <- make_fixed_effect_plot(
  df = Post_natal_ARI,
  ylab = "Postnatal Climate Anomalies",
  panel_title = "B"
)

fixed_effect_figure <- ggarrange(
  Pre_FE,
  Post_FE,
  ncol = 2,
  common.legend = TRUE,
  legend = "bottom"
)

save_tiff(
  filename = "Figure_1.tiff",
  plot = fixed_effect_figure,
  width = 16,
  height = 10
)


# ==============================================================================
# 4. Load Nonlinear/GAM Outputs
# ==============================================================================

prenatal_gam_file <- here(
  "outputs",
  "tables",
  "Prenatal_GAM.xlsx"
)

postnatal_gam_file <- here(
  "outputs",
  "tables",
  "Postnatal_GAM.xlsx"
)

prenatal_gam <- read_excel_workbook(
  prenatal_gam_file
)

postnatal_gam <- read_excel_workbook(
  postnatal_gam_file
)

prenatal_gam <- lapply(
  prenatal_gam,
  standardize_result_columns
)

postnatal_gam <- lapply(
  postnatal_gam,
  standardize_result_columns
)


# ==============================================================================
# 5. Nonlinear/GAM Plot Function
# ==============================================================================

make_gam_plot <- function(
    df,
    xlab,
    title
) {
  
  required_cols <- c(
    "ID",
    "mean",
    "lci",
    "uci"
  )
  
  missing_cols <- setdiff(
    required_cols,
    names(df)
  )
  
  if (length(missing_cols) > 0) {
    stop(
      "GAM result data frame is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  ggplot(
    data = df,
    aes(
      x = ID,
      y = mean,
      ymin = lci,
      ymax = uci
    )
  ) +
    geom_line(
      linewidth = 0.8
    ) +
    geom_ribbon(
      alpha = 0.12
    ) +
    xlab(
      xlab
    ) +
    ylab(
      "Estimated effect (HR)"
    ) +
    labs(
      title = title
    ) +
    theme_pubr()
}

get_sheet <- function(
    workbook_list,
    sheet_name
) {
  
  if (!sheet_name %in% names(workbook_list)) {
    stop(
      "Expected sheet not found: ",
      sheet_name,
      ". Available sheets are: ",
      paste(names(workbook_list), collapse = ", ")
    )
  }
  
  workbook_list[[sheet_name]]
}

temperature_xlab <- as.expression(
  expression(
    paste(
      "Temperature anomaly (",
      sd,
      ")"
    )
  )
)

mean_temperature_xlab <- as.expression(
  expression(
    paste(
      "Mean temperature anomaly (",
      sd,
      ")"
    )
  )
)

max_temperature_xlab <- as.expression(
  expression(
    paste(
      "Maximum temperature anomaly (",
      sd,
      ")"
    )
  )
)

precipitation_xlab <- as.expression(
  expression(
    paste(
      "Precipitation anomaly (",
      sd,
      ")"
    )
  )
)


# ==============================================================================
# 6. Prenatal Nonlinear Figures
# ==============================================================================

T1 <- make_gam_plot(
  df = get_sheet(prenatal_gam, "IN_GAM1_TEMP"),
  xlab = mean_temperature_xlab,
  title = "Model 1"
)

T2 <- make_gam_plot(
  df = get_sheet(prenatal_gam, "IN_GAM3_TEMP"),
  xlab = mean_temperature_xlab,
  title = "Model 2"
)

T3 <- make_gam_plot(
  df = get_sheet(prenatal_gam, "IN_GAM4_TEMP"),
  xlab = mean_temperature_xlab,
  title = "Model 3"
)

T4 <- make_gam_plot(
  df = get_sheet(prenatal_gam, "IN_GAM5_TEMP"),
  xlab = mean_temperature_xlab,
  title = "Model 4"
)

TM1 <- make_gam_plot(
  df = get_sheet(prenatal_gam, "IN_GAM1_TMX"),
  xlab = max_temperature_xlab,
  title = "Model 1"
)

TM2 <- make_gam_plot(
  df = get_sheet(prenatal_gam, "IN_GAM3_TMX"),
  xlab = max_temperature_xlab,
  title = "Model 2"
)

TM3 <- make_gam_plot(
  df = get_sheet(prenatal_gam, "IN_GAM4_TMX"),
  xlab = max_temperature_xlab,
  title = "Model 3"
)

TM4 <- make_gam_plot(
  df = get_sheet(prenatal_gam, "IN_GAM5_TMX"),
  xlab = max_temperature_xlab,
  title = "Model 4"
)

P1 <- make_gam_plot(
  df = get_sheet(prenatal_gam, "IN_GAM2_RAIN"),
  xlab = precipitation_xlab,
  title = "Model 1"
)

P2 <- make_gam_plot(
  df = get_sheet(prenatal_gam, "IN_GAM3_RAIN"),
  xlab = precipitation_xlab,
  title = "Model 2"
)

P3 <- make_gam_plot(
  df = get_sheet(prenatal_gam, "IN_GAM4_RAIN"),
  xlab = precipitation_xlab,
  title = "Model 3"
)

P4 <- make_gam_plot(
  df = get_sheet(prenatal_gam, "IN_GAM5_RAIN"),
  xlab = precipitation_xlab,
  title = "Model 4"
)


# ==============================================================================
# 7. Postnatal Nonlinear Figures
# ==============================================================================

T_1a <- make_gam_plot(
  df = get_sheet(postnatal_gam, "PU_GAM1_TEMP"),
  xlab = mean_temperature_xlab,
  title = "Model 1"
)

T_2a <- make_gam_plot(
  df = get_sheet(postnatal_gam, "PU_GAM3_TEMP"),
  xlab = mean_temperature_xlab,
  title = "Model 2"
)

T_3a <- make_gam_plot(
  df = get_sheet(postnatal_gam, "PU_GAM4_TEMP"),
  xlab = mean_temperature_xlab,
  title = "Model 3"
)

T_4a <- make_gam_plot(
  df = get_sheet(postnatal_gam, "PU_GAM5_TEMP"),
  xlab = mean_temperature_xlab,
  title = "Model 4"
)

TM_1a <- make_gam_plot(
  df = get_sheet(postnatal_gam, "PU_GAM1_TMX"),
  xlab = max_temperature_xlab,
  title = "Model 1"
)

TM_2a <- make_gam_plot(
  df = get_sheet(postnatal_gam, "PU_GAM3_TMX"),
  xlab = max_temperature_xlab,
  title = "Model 2"
)

TM_3a <- make_gam_plot(
  df = get_sheet(postnatal_gam, "PU_GAM4_TMX"),
  xlab = max_temperature_xlab,
  title = "Model 3"
)

TM_4a <- make_gam_plot(
  df = get_sheet(postnatal_gam, "PU_GAM5_TMX"),
  xlab = max_temperature_xlab,
  title = "Model 4"
)

P_1a <- make_gam_plot(
  df = get_sheet(postnatal_gam, "PU_GAM2_RAIN"),
  xlab = precipitation_xlab,
  title = "Model 1"
)

P_2a <- make_gam_plot(
  df = get_sheet(postnatal_gam, "PU_GAM3_RAIN"),
  xlab = precipitation_xlab,
  title = "Model 2"
)

P_3a <- make_gam_plot(
  df = get_sheet(postnatal_gam, "PU_GAM4_RAIN"),
  xlab = precipitation_xlab,
  title = "Model 3"
)

P_4a <- make_gam_plot(
  df = get_sheet(postnatal_gam, "PU_GAM5_RAIN"),
  xlab = precipitation_xlab,
  title = "Model 4"
)


# ==============================================================================
# 8. Export Nonlinear/GAM Figures
# ==============================================================================

Figure_B1 <- ggarrange(
  T1,
  T2,
  T3,
  T4,
  ncol = 2,
  nrow = 2
)

Figure_B2 <- ggarrange(
  T_1a,
  T_2a,
  T_3a,
  T_4a,
  ncol = 2,
  nrow = 2
)

Figure_B3 <- ggarrange(
  TM1,
  TM2,
  TM3,
  TM4,
  ncol = 2,
  nrow = 2
)

Figure_B4 <- ggarrange(
  TM_1a,
  TM_2a,
  TM_3a,
  TM_4a,
  ncol = 2,
  nrow = 2
)

Figure_B5 <- ggarrange(
  P1,
  P2,
  P3,
  P4,
  ncol = 2,
  nrow = 2
)

Figure_B6 <- ggarrange(
  P_1a,
  P_2a,
  P_3a,
  P_4a,
  ncol = 2,
  nrow = 2
)

save_tiff(
  filename = "Figure_B1.tiff",
  plot = Figure_B1,
  width = 10,
  height = 8
)

save_tiff(
  filename = "Figure_B2.tiff",
  plot = Figure_B2,
  width = 10,
  height = 8
)

save_tiff(
  filename = "Figure_B3.tiff",
  plot = Figure_B3,
  width = 10,
  height = 8
)

save_tiff(
  filename = "Figure_B4.tiff",
  plot = Figure_B4,
  width = 10,
  height = 8
)

save_tiff(
  filename = "Figure_B5.tiff",
  plot = Figure_B5,
  width = 10,
  height = 8
)

save_tiff(
  filename = "Figure_B6.tiff",
  plot = Figure_B6,
  width = 10,
  height = 8
)




# ==============================================================================
# 9. Session Information
# ==============================================================================

writeLines(
  capture.output(
    sessionInfo()
  ),
  con = here(
    "outputs",
    "figures",
    "sessionInfo_figures_linear_non_linear.txt"
  )
)

cat(
  "\nLinear and nonlinear figures generated successfully.\n"
)

# ==============================================================================
# END
# ==============================================================================