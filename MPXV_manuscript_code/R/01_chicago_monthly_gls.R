# ============================================================================
# Chicago monthly epidemiological model (first wave, June 2022 - March 2023)
#
# Predictors: JYNNEOS first and second doses (CDPH), total diagnostic tests in
# the MSM population (Howard Brown Health), and self-reported STI positivity
# rate and unprotected/protected sex acts ratio (RADAR cohort).
#
# GLS AR(1), OLS and negative binomial models are compared by AIC/BIC across
# combinations of variables and 1- and 2-month vaccine lags. The final model
# is GLS AR(1) with all variables; the 1- and 2-month lags are reported as
# complementary analyses.
# ============================================================================

source(file.path("R", "00_functions.R"))

monthly <- read.csv(file.path(DATA_DIR, "Chi_EpiData_2022_23.csv"))
monthly$Month <- as.Date(monthly$Month)
monthly <- monthly %>%
  filter(Month < as.Date("2023-04-01")) %>%
  arrange(Month) %>%
  prepare_series()

doses <- c("Dose1_lag", "Dose2_lag")
variable_sets <- list(
  doses_only      = doses,
  doses_sti       = c(doses, "STIPositivityRate"),
  doses_tests     = c(doses, "Total_Tests"),
  all_variables   = c(doses, "STIPositivityRate", "unprotected_protected", "Total_Tests")
)
lags <- c(1, 2)

comparison <- list()
final_fits <- list()

for (k in lags) {
  d  <- add_dose_lags(monthly, k)
  dw <- dw_test(doses, d)
  for (v in names(variable_sets)) {
    models <- fit_model_set(variable_sets[[v]], d)
    comparison[[paste(k, v)]] <- cbind(lag_months = k, variables = v,
                                       compare_model_set(models),
                                       DW_p_value = dw$p)
    if (v == "all_variables") final_fits[[as.character(k)]] <- models[["GLS AR(1)"]]
  }
}

comparison_table <- bind_rows(comparison)
coef_table <- bind_rows(lapply(lags, function(k) {
  cbind(lag_months = k, tidy_coefs(final_fits[[as.character(k)]]))
}))

write_table(comparison_table, "chicago_monthly_model_comparison.csv")
write_table(coef_table,       "chicago_monthly_final_model_coefficients.csv")

cat("\n==== Chicago monthly: model comparison ====\n")
print(comparison_table[, c("lag_months", "variables", "model_type", "n",
                           "AIC", "BIC", "AR1_phi", "DW_p_value")],
      digits = 4, row.names = FALSE)

for (k in lags) {
  cat("\n---- Final model: GLS AR(1), all variables,", k, "month lag ----\n")
  print(summary(final_fits[[as.character(k)]]))
}
