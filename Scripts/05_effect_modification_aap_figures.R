# ==============================================================================
# 05_effect_modification_figures.R
# ==============================================================================
# Title: Climate, Air Pollution, Effect Modification and Child ARI
# Purpose:
#   Generate manuscript Figures 3-5 and Supplementary Figures D1-D5 from saved
#   effect-modification model outputs.
#
# Inputs:
#   outputs/effect_modification/results_unadj.rds
#   outputs/effect_modification/results_adj.rds
#
# Outputs:
#   outputs/figures/Figure_3.tif
#   outputs/figures/Figure_4.tif
#   outputs/figures/Figure_5.tif
#   outputs/figures/Figure_D1.tif
#   outputs/figures/Figure_D2.tif
#   outputs/figures/Figure_D3.tif
#   outputs/figures/Figure_D4.tif
#   outputs/figures/Figure_D5.tif
#   outputs/effect_modification/sessionInfo_figures.txt
# ==============================================================================

# 0. Setup ---------------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(cowplot)
library(here)

dir.create(here("outputs"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("outputs", "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("outputs", "effect_modification"), recursive = TRUE, showWarnings = FALSE)

# 1. Load Saved Results --------------------------------------------------------

results_unadj <- readRDS(
  here("outputs", "effect_modification", "results_unadj.rds")
)

results_adj <- readRDS(
  here("outputs", "effect_modification", "results_adj.rds")
)

# 2. Definitions ---------------------------------------------------------------

pollutants <- c("PM", "NO", "SO", "CO", "O3")

figure_pollutants <- list(
  D1 = "PM",
  D2 = "CO",
  D3 = "NO",
  D4 = "SO",
  D5 = "O3"
)

climate_labels <- list(
  TMP = "Mean Temperature Anomaly",
  TMX = "Maximum Temperature Anomaly",
  PRE = "Precipitation Anomaly"
)

# 3. Plot Function -------------------------------------------------------------

create_plot <- function(df, title, xlab, show_legend = FALSE, tag = NULL) {
  
  ggplot(
    df,
    aes(
      x = ID,
      y = mean,
      color = Pollutant_Level,
      fill = Pollutant_Level
    )
  ) +
    geom_line(linewidth = 1) +
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
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = ifelse(show_legend, "bottom", "none"),
      plot.title = element_text(hjust = 0.5)
    )
}

make_plot <- function(result_list, prefix, period, pollutant, tag = NULL) {
  
  id <- paste0(prefix, "_", period, "_", pollutant)
  
  if (!id %in% names(result_list)) {
    stop("Missing result object: ", id)
  }
  
  create_plot(
    df = result_list[[id]],
    title = ifelse(period == "IU", "Prenatal", "Postnatal"),
    xlab = climate_labels[[prefix]],
    show_legend = FALSE,
    tag = tag
  )
}

save_plot_grid <- function(filename, plot_list, legend, width = 8, height = 12) {
  
  tiff(
    here("outputs", "figures", filename),
    units = "in",
    width = width,
    height = height,
    res = 300,
    compression = "lzw"
  )
  
  print(
    plot_grid(
      plot_grid(
        plotlist = c(plot_list, list(legend)),
        ncol = 2
      ),
      rel_widths = c(3, 0.4)
    )
  )
  
  dev.off()
}

# 4. Shared Legend -------------------------------------------------------------

legend_source <- results_adj[["PRE_IU_O3"]]

if (is.null(legend_source)) {
  stop("Legend source PRE_IU_O3 not found in adjusted results.")
}

legend <- get_legend(
  create_plot(
    df = legend_source,
    title = "",
    xlab = "",
    show_legend = TRUE
  ) +
    theme(
      legend.title = element_blank()
    )
)

# 5. Main Figures 3-5 ----------------------------------------------------------

# Figure 3: Adjusted mean temperature anomaly models
figure_3_plots <- list(
  make_plot(results_adj, "TMP", "IU", "PM", "A)"),
  make_plot(results_adj, "TMP", "PU", "PM"),
  make_plot(results_adj, "TMP", "IU", "NO"),
  make_plot(results_adj, "TMP", "PU", "NO"),
  make_plot(results_adj, "TMP", "IU", "SO"),
  make_plot(results_adj, "TMP", "PU", "SO"),
  make_plot(results_adj, "TMP", "IU", "CO"),
  make_plot(results_adj, "TMP", "PU", "CO"),
  make_plot(results_adj, "TMP", "IU", "O3"),
  make_plot(results_adj, "TMP", "PU", "O3")
)

save_plot_grid(
  filename = "Figure_3.tif",
  plot_list = figure_3_plots,
  legend = legend,
  width = 8,
  height = 12
)

# Figure 4: Adjusted maximum temperature anomaly models
figure_4_plots <- list(
  make_plot(results_adj, "TMX", "IU", "PM", "C)"),
  make_plot(results_adj, "TMX", "PU", "PM"),
  make_plot(results_adj, "TMX", "IU", "NO"),
  make_plot(results_adj, "TMX", "PU", "NO"),
  make_plot(results_adj, "TMX", "IU", "SO"),
  make_plot(results_adj, "TMX", "PU", "SO"),
  make_plot(results_adj, "TMX", "IU", "CO"),
  make_plot(results_adj, "TMX", "PU", "CO"),
  make_plot(results_adj, "TMX", "IU", "O3"),
  make_plot(results_adj, "TMX", "PU", "O3")
)

save_plot_grid(
  filename = "Figure_4.tif",
  plot_list = figure_4_plots,
  legend = legend,
  width = 8,
  height = 12
)

# Figure 5: Adjusted precipitation anomaly models
figure_5_plots <- list(
  make_plot(results_adj, "PRE", "IU", "PM", "E)"),
  make_plot(results_adj, "PRE", "PU", "PM"),
  make_plot(results_adj, "PRE", "IU", "NO"),
  make_plot(results_adj, "PRE", "PU", "NO"),
  make_plot(results_adj, "PRE", "IU", "SO"),
  make_plot(results_adj, "PRE", "PU", "SO"),
  make_plot(results_adj, "PRE", "IU", "CO"),
  make_plot(results_adj, "PRE", "PU", "CO"),
  make_plot(results_adj, "PRE", "IU", "O3"),
  make_plot(results_adj, "PRE", "PU", "O3")
)

save_plot_grid(
  filename = "Figure_5.tif",
  plot_list = figure_5_plots,
  legend = legend,
  width = 8,
  height = 12
)

# 6. Supplementary Figures D1-D5 ----------------------------------------------

make_supplementary_plot_list <- function(pollutant) {
  
  list(
    make_plot(results_unadj, "TMP", "IU", pollutant, "A)"),
    make_plot(results_unadj, "TMP", "PU", pollutant),
    make_plot(results_adj, "TMP", "IU", pollutant, "B)"),
    make_plot(results_adj, "TMP", "PU", pollutant),
    
    make_plot(results_unadj, "TMX", "IU", pollutant, "C)"),
    make_plot(results_unadj, "TMX", "PU", pollutant),
    make_plot(results_adj, "TMX", "IU", pollutant, "D)"),
    make_plot(results_adj, "TMX", "PU", pollutant),
    
    make_plot(results_unadj, "PRE", "IU", pollutant, "E)"),
    make_plot(results_unadj, "PRE", "PU", pollutant),
    make_plot(results_adj, "PRE", "IU", pollutant, "F)"),
    make_plot(results_adj, "PRE", "PU", pollutant)
  )
}

for (fig_id in names(figure_pollutants)) {
  
  pollutant <- figure_pollutants[[fig_id]]
  
  supplementary_plots <- make_supplementary_plot_list(
    pollutant = pollutant
  )
  
  save_plot_grid(
    filename = paste0("Figure_", fig_id, ".tif"),
    plot_list = supplementary_plots,
    legend = legend,
    width = 8,
    height = 12
  )
}

# 7. Session Information -------------------------------------------------------

writeLines(
  capture.output(sessionInfo()),
  here("outputs", "effect_modification", "sessionInfo_figures.txt")
)

cat("\nFigures 3-5 and D1-D5 generated successfully.\n")
# ==============================================================================
# END
# ==============================================================================