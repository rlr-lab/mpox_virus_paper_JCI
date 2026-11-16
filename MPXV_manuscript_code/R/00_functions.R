# ============================================================================
# Shared helper functions for the MPXV vaccination GLS analyses
# ============================================================================

suppressPackageStartupMessages({
  library(nlme)    # gls(), corAR1()
  library(MASS)    # glm.nb()
  library(dplyr)
  library(car)     # durbinWatsonTest()
  library(ggplot2)
})

DATA_DIR    <- "data"
TABLE_DIR   <- file.path("results", "tables")
FIGURE_DIR  <- file.path("results", "figures")
dir.create(TABLE_DIR,  recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)

# Add a sequential time index and the log-transformed outcome.
# Must be called on the full, time-ordered series *before* lagging so that
# the AR(1) structure uses the true time spacing.
prepare_series <- function(df) {
  df$time_index <- seq_len(nrow(df))
  df$log_cases  <- log(df$MPXVcases + 1)
  df
}

# Add lagged vaccine dose columns (Dose1_lag, Dose2_lag) and drop the
# leading rows that become NA.
add_dose_lags <- function(df, k) {
  df$Dose1_lag <- dplyr::lag(df$Dose1, k)
  df$Dose2_lag <- dplyr::lag(df$Dose2, k)
  df[complete.cases(df[, c("log_cases", "Dose1_lag", "Dose2_lag")]), ]
}

# GLS with AR(1) residual correlation, fit by maximum likelihood so that
# AIC is comparable across models with different fixed effects.
fit_gls_ar1 <- function(formula, data) {
  gls(formula,
      data        = data,
      correlation = corAR1(form = ~ time_index),
      method      = "ML")
}

# Fit the three model types compared in every analysis:
#   GLS with AR(1) correlation and OLS on log(cases + 1), and negative
#   binomial regression on the raw counts. A model that fails to fit is NULL.
fit_model_set <- function(predictors, data) {
  rhs <- paste(predictors, collapse = " + ")
  log_formula   <- as.formula(paste("log_cases ~", rhs))
  count_formula <- as.formula(paste("MPXVcases ~", rhs))
  try_fit <- function(expr) tryCatch(expr, error = function(e) NULL)
  list(
    "GLS AR(1)"         = try_fit(fit_gls_ar1(log_formula, data)),
    "OLS"               = try_fit(lm(log_formula, data = data)),
    "Negative binomial" = try_fit(glm.nb(count_formula, data = data))
  )
}

# Durbin-Watson test on the OLS fit of `predictors` (bootstrap p-value).
dw_test <- function(predictors, data) {
  set.seed(2022)
  durbinWatsonTest(lm(as.formula(paste("log_cases ~",
                                       paste(predictors, collapse = " + "))),
                      data = data))
}

# One row per fitted model: fit statistics plus vaccine dose estimates.
compare_model_set <- function(models) {
  bind_rows(lapply(names(models), function(type) {
    m <- models[[type]]
    if (is.null(m)) return(data.frame(model_type = type, converged = FALSE))
    cf <- tidy_coefs(m)
    cbind(model_type = type, converged = TRUE, glance_model(m),
          Dose1_estimate = cf$estimate[cf$term == "Dose1_lag"],
          Dose1_p_value  = cf$p_value[cf$term == "Dose1_lag"],
          Dose2_estimate = cf$estimate[cf$term == "Dose2_lag"],
          Dose2_p_value  = cf$p_value[cf$term == "Dose2_lag"])
  }))
}

# Estimated AR(1) autocorrelation parameter (phi).
ar1_phi <- function(model) {
  unname(coef(model$modelStruct$corStruct, unconstrained = FALSE))
}

# Coefficient table as a data frame for any of gls / lm / glm.nb.
tidy_coefs <- function(model) {
  tab <- if (inherits(model, "gls")) summary(model)$tTable else coef(summary(model))
  out <- data.frame(term = rownames(tab), tab[, 1:4], row.names = NULL)
  names(out) <- c("term", "estimate", "std_error", "statistic", "p_value")
  out
}

# One-row model fit summary.
glance_model <- function(model) {
  data.frame(
    n       = nobs(model),
    AIC     = AIC(model),
    BIC     = BIC(model),
    logLik  = as.numeric(logLik(model)),
    AR1_phi = if (inherits(model, "gls")) ar1_phi(model) else NA_real_
  )
}

write_table <- function(df, name) {
  write.csv(df, file.path(TABLE_DIR, name), row.names = FALSE)
}
