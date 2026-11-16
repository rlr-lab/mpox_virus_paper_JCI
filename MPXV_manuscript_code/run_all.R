# ============================================================================
# Reproduce all final models. Run from the repository root:
#   Rscript run_all.R
# ============================================================================

scripts <- c(
  "R/01_chicago_monthly_gls.R",
  "R/02_chicago_weekly_gls.R",
  "R/03_three_cities_gls.R",
  "R/04_sexual_partners.R"
)

for (s in scripts) {
  message("\n>>> Running ", s)
  source(s, echo = FALSE, local = new.env())
}

writeLines(capture.output(sessionInfo()), file.path("results", "sessionInfo.txt"))
message("\nDone. Tables in results/tables, figures in results/figures.")
