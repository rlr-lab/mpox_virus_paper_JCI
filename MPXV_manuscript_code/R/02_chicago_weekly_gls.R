# ============================================================================
# Chicago weekly model (CDPH weekly cases and vaccine doses, first wave)
#
# Predictors: JYNNEOS first and second doses.
# GLS AR(1), OLS and negative binomial models are compared by AIC/BIC across
# vaccine lags of 1-4 weeks; the model with the lowest AIC is retained
# (GLS AR(1), 4-week lag).
# ============================================================================

source(file.path("R", "00_functions.R"))

weekly <- read.csv(file.path(DATA_DIR, "AllCDPHData_FirstWave_Granular.csv"))
weekly$date <- as.Date(weekly$date)
weekly <- weekly %>% arrange(date) %>% prepare_series()

predictors <- c("Dose1_lag", "Dose2_lag")
lag_weeks  <- 1:4

comparison <- list()
fits       <- list()

for (k in lag_weeks) {
  d  <- add_dose_lags(weekly, k)
  dw <- dw_test(predictors, d)
  fits[[k]] <- fit_model_set(predictors, d)
  comparison[[k]] <- cbind(lag_weeks = k, compare_model_set(fits[[k]]),
                           DW_p_value = dw$p)
}
comparison_table <- bind_rows(comparison)

# ---- Final model: lowest AIC ------------------------------------------------
best        <- comparison_table[which.min(comparison_table$AIC), ]
final_model <- fits[[best$lag_weeks]][[best$model_type]]

write_table(comparison_table, "chicago_weekly_model_comparison.csv")
write_table(cbind(lag_weeks = best$lag_weeks, model_type = best$model_type,
                  tidy_coefs(final_model)),
            "chicago_weekly_final_model_coefficients.csv")

cat("\n==== Chicago weekly: model comparison ====\n")
print(comparison_table[, c("lag_weeks", "model_type", "n", "AIC", "BIC",
                           "AR1_phi", "Dose1_p_value", "Dose2_p_value", "DW_p_value")],
      digits = 4, row.names = FALSE)
cat("\nSelected:", best$model_type, "with", best$lag_weeks, "week lag\n")
print(summary(final_model))
