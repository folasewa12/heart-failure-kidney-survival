# Run the full analysis. Tables and charts are written to outputs/.
source("R/00_setup.R")
source("R/01_prepare_data.R")
source("R/02_survival_descriptives.R")
source("R/03_cox_models.R")
message("Done - see the outputs/ folder.")
