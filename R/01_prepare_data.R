# =============================================================================
# 01_prepare_data.R  -  data checks, eGFR and analysis variables
# =============================================================================
# Output: `hf` - 299 patients, 96 deaths, follow-up 4-285 days

hf <- read.csv(data_file)

# ---- 1. Basic checks ------------------------------------------------------------
stopifnot(nrow(hf) == 299, !anyNA(hf), !any(duplicated(hf)))
summary(hf)

# ---- 2. Values that look filled in rather than measured ---------------------
# The exact column mean appears repeatedly, which suggests missing values were
# replaced by the mean before the data was published. These are flagged, kept in
# the main analysis, and dropped in a sensitivity analysis (script 03).
data_quality <- tibble::tibble(
  check = c("Platelets equal to the column mean (263,358.03)",
            "Creatinine phosphokinase equal to 582 (rounded mean)",
            "Serum creatinine equal to 1.18 or 1.83 mg/dL",
            "Non-integer age (60.667)"),
  rows  = c(sum(hf$platelets == 263358.03),
            sum(hf$creatinine_phosphokinase == 582),
            sum(hf$serum_creatinine %in% c(1.18, 1.83)),
            sum(hf$age %% 1 != 0))
)
print(data_quality)
write.csv(data_quality, file.path(table_dir, "data_quality_checks.csv"), row.names = FALSE)

# ---- 3. Kidney function: eGFR (CKD-EPI 2021, race-free equation) -------------
# Serum creatinine is in mg/dL. sex: 1 = male, 0 = female.
ckd_epi_2021 <- function(scr, age, female) {
  kappa <- ifelse(female, 0.7, 0.9)
  alpha <- ifelse(female, -0.241, -0.302)
  142 * pmin(scr / kappa, 1)^alpha * pmax(scr / kappa, 1)^-1.200 *
    0.9938^age * ifelse(female, 1.012, 1)
}

hf <- hf %>%
  mutate(
    female   = sex == 0,
    egfr     = ckd_epi_2021(serum_creatinine, age, female),
    egfr_grp = cut(egfr, breaks = c(-Inf, 30, 45, 60, Inf), right = FALSE,
                   labels = c("<30", "30-44", "45-59", ">=60")),
    egfr_grp = relevel(egfr_grp, ref = ">=60"),
    # Scaled so hazard ratios read naturally
    egfr_per10_lower   = -egfr / 10,               # per 10 mL/min/1.73m2 LOWER eGFR
    age_per10          = age / 10,                  # per 10 years older
    ef_per5_lower      = -ejection_fraction / 5,    # per 5 points LOWER ejection fraction
    sodium_per5_lower  = -serum_sodium / 5,         # per 5 mmol/L LOWER sodium
    suspected_imputed  = serum_creatinine %in% c(1.18, 1.83)
  )

message("Patients: ", nrow(hf), " | deaths: ", sum(hf$DEATH_EVENT),
        " | median eGFR: ", round(median(hf$egfr), 1))
