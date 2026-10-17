# =============================================================================
# 00_setup.R  -  packages, paths and shared settings
# Project: Kidney function and survival in heart failure
# Author:  Afolasewa Duduyemi
# =============================================================================
# Run from the project root (open the folder in RStudio, or setwd() to it).

required_packages <- c("survival", "dplyr", "ggplot2", "broom", "tidyr", "scales")
missing <- required_packages[!required_packages %in% rownames(installed.packages())]
if (length(missing) > 0) install.packages(missing)
invisible(lapply(required_packages, library, character.only = TRUE))

data_file  <- file.path("data", "heart_failure_clinical_records_dataset.csv")
figure_dir <- file.path("outputs", "figures")
table_dir  <- file.path("outputs", "tables")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir,  recursive = TRUE, showWarnings = FALSE)

# One fixed colour per eGFR group, ordered from best to worst kidney function
egfr_colours <- c(">=60" = "#2a78d6", "45-59" = "#1baf7a", "30-44" = "#eda100", "<30" = "#e34948")

theme_report <- function() {
  theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"),
          panel.grid.minor = element_blank(),
          legend.position = "bottom")
}

save_plot <- function(plot, name, width = 8, height = 5) {
  ggsave(file.path(figure_dir, paste0(name, ".png")), plot,
         width = width, height = height, dpi = 300, bg = "white")
}
