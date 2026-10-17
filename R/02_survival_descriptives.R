# =============================================================================
# 02_survival_descriptives.R  -  mortality by kidney function, Kaplan-Meier, log-rank
# =============================================================================

# ---- 1. Deaths and mortality rate by eGFR group ------------------------------
mortality_table <- hf %>%
  group_by(egfr_grp) %>%
  summarise(patients      = n(),
            deaths        = sum(DEATH_EVENT),
            pct_died      = 100 * deaths / patients,
            person_years  = sum(time) / 365.25,
            deaths_per_100_person_years = 100 * deaths / person_years,
            .groups = "drop")
print(mortality_table)
write.csv(mortality_table, file.path(table_dir, "mortality_by_egfr.csv"), row.names = FALSE)

# ---- 2. Kaplan-Meier survival --------------------------------------------------
km_fit <- survfit(Surv(time, DEATH_EVENT) ~ egfr_grp, data = hf)
summary(km_fit, times = c(30, 90, 180))

km_data <- tidy(km_fit) %>%
  mutate(egfr_grp = sub("egfr_grp=", "", strata))

fig_km <- ggplot(km_data, aes(time, estimate, colour = egfr_grp)) +
  geom_step(linewidth = 0.9) +
  scale_colour_manual(values = egfr_colours, name = "eGFR (mL/min/1.73m²)") +
  scale_y_continuous(labels = percent, limits = c(0, 1)) +
  labs(title = "Survival falls sharply once eGFR is below 45",
       x = "Days of follow-up", y = "Surviving") +
  theme_report()
save_plot(fig_km, "kaplan_meier_by_egfr")

# ---- 3. Log-rank test across the four groups ----------------------------------
survdiff(Surv(time, DEATH_EVENT) ~ egfr_grp, data = hf)
