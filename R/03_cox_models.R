# =============================================================================
# 03_cox_models.R  -  Cox proportional hazards models, checks and sensitivity analysis
# =============================================================================

adjusters <- "age_per10 + sex + ef_per5_lower + sodium_per5_lower +
              anaemia + diabetes + high_blood_pressure + smoking"

fit_cox <- function(rhs, data = hf) {
  coxph(as.formula(paste("Surv(time, DEATH_EVENT) ~", rhs)), data = data)
}

models <- list(
  unadjusted_continuous = fit_cox("egfr_per10_lower"),
  unadjusted_groups     = fit_cox("egfr_grp"),
  adjusted_continuous   = fit_cox(paste("egfr_per10_lower +", adjusters)),
  adjusted_groups       = fit_cox(paste("egfr_grp +", adjusters)),
  # Sensitivity: drop the 19 patients whose creatinine looks filled in
  sensitivity_no_imputed = fit_cox(paste("egfr_per10_lower +", adjusters),
                                   data = filter(hf, !suspected_imputed)),
  # Sensitivity: creatinine itself, since age is also inside the eGFR equation
  sensitivity_log_creatinine = fit_cox(paste("log(serum_creatinine) +", adjusters))
)
lapply(models, summary)

# ---- Hazard ratio table ---------------------------------------------------------
hr_table <- bind_rows(lapply(names(models), function(m) {
  tidy(models[[m]], exponentiate = TRUE, conf.int = TRUE) %>% mutate(model = m)
})) %>%
  select(model, term, hazard_ratio = estimate, conf.low, conf.high, p.value)
print(hr_table, n = Inf)
write.csv(hr_table, file.path(table_dir, "cox_hazard_ratios.csv"), row.names = FALSE)

# ---- Checks -------------------------------------------------------------------
# Proportional hazards: a significant global test would mean an effect changes over time
cox.zph(models$adjusted_continuous)

# Linearity of eGFR: separate slopes below and above 60
hf <- hf %>%
  mutate(egfr_below60_per10 = pmax(60 - egfr, 0) / 10,
         egfr_above60_per10 = pmin(60 - egfr, 0) / 10)
piecewise <- fit_cox(paste("egfr_below60_per10 + egfr_above60_per10 +", adjusters))
summary(piecewise)

# ---- Forest plot of the adjusted model ----------------------------------------
labels <- c(egfr_per10_lower = "eGFR, per 10 lower", age_per10 = "Age, per 10 years",
            ef_per5_lower = "Ejection fraction, per 5 lower", high_blood_pressure = "High blood pressure",
            anaemia = "Anaemia", sodium_per5_lower = "Sodium, per 5 mmol/L lower",
            smoking = "Smoking", diabetes = "Diabetes", sex = "Male sex")

forest <- tidy(models$adjusted_continuous, exponentiate = TRUE, conf.int = TRUE) %>%
  mutate(label = factor(labels[term], levels = rev(labels)))

fig_forest <- ggplot(forest, aes(estimate, label)) +
  geom_vline(xintercept = 1, colour = "grey60") +
  geom_errorbar(aes(xmin = conf.low, xmax = conf.high), width = 0.2, orientation = "y") +
  geom_point(size = 2.5, colour = "#2a78d6") +
  scale_x_log10() +
  labs(title = "Lower eGFR predicts death independently of age and heart function",
       x = "Adjusted hazard ratio (95% CI, log scale)", y = NULL) +
  theme_report()
save_plot(fig_forest, "adjusted_hazard_ratios")
