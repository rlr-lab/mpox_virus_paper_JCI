# ============================================================================
# Three-city model: Chicago, Los Angeles and New York (weekly, first wave)
#
# Predictors: JYNNEOS first and second doses.
# Data: CDPH (Chicago), NYC Department of Health GitHub repository (New York)
# and Los Angeles County Department of Public Health (Los Angeles).
#
# For each city, GLS AR(1), OLS and negative binomial models are compared by
# AIC/BIC across vaccine lags of 1-4 weeks; the model with the lowest AIC is
# retained (GLS AR(1), 4-week lag, in all three cities).
# ============================================================================

source(file.path("R", "00_functions.R"))

cities_data <- read.csv(file.path(DATA_DIR, "AllCities_FirstWave.csv"))
cities_data$Week <- as.Date(cities_data$Week)

cities     <- c("Chicago", "LA", "NY")
predictors <- c("Dose1_lag", "Dose2_lag")
lag_weeks  <- 1:4

comparison <- list()
coefs      <- list()

for (city in cities) {
  series <- cities_data %>% filter(City == city) %>% arrange(Week) %>% prepare_series()

  fits <- list()
  city_rows <- list()
  for (k in lag_weeks) {
    d  <- add_dose_lags(series, k)
    dw <- dw_test(predictors, d)
    fits[[k]] <- fit_model_set(predictors, d)
    city_rows[[k]] <- cbind(City = city, lag_weeks = k,
                            compare_model_set(fits[[k]]), DW_p_value = dw$p)
  }
  city_table <- bind_rows(city_rows)
  comparison[[city]] <- city_table

  # ---- Final model: lowest AIC ----
  best        <- city_table[which.min(city_table$AIC), ]
  final_model <- fits[[best$lag_weeks]][[best$model_type]]
  coefs[[city]] <- cbind(City = city, lag_weeks = best$lag_weeks,
                         model_type = best$model_type, tidy_coefs(final_model))

  cat("\n==== ", city, ": ", best$model_type, ", ", best$lag_weeks,
      "-week lag ====\n", sep = "")
  print(summary(final_model))
}

comparison_table <- bind_rows(comparison)
coef_table       <- bind_rows(coefs)

write_table(comparison_table, "three_cities_model_comparison.csv")
write_table(coef_table,       "three_cities_final_model_coefficients.csv")

cat("\n==== Three cities: model comparison ====\n")
print(comparison_table[, c("City", "lag_weeks", "model_type", "n", "AIC", "BIC",
                           "AR1_phi", "Dose1_p_value", "Dose2_p_value", "DW_p_value")],
      digits = 4, row.names = FALSE)
