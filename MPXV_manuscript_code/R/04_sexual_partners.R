# ============================================================================
# Self-reported sexual partner numbers in the RADAR cohort (2022)
#
# Outcome: number of vaginal/anal sex partners in the past 6 months
# (vasexpartners_p6m), from participants with more than one visit in 2022.
#
# Longitudinal change: negative binomial generalized linear mixed-effects
# models with participant-level random intercepts, with time as
#   (a) months since the participant's first 2022 visit and
#   (b) pre vs post mpox outbreak onset (4 June 2022).
# ============================================================================

source(file.path("R", "00_functions.R"))
suppressPackageStartupMessages(library(lme4))

OUTBREAK_START <- as.Date("2022-06-04")

partners <- read.csv(file.path(DATA_DIR, "RADAR_STI_2022_NonUnique_Partners.csv"))
partners <- partners %>%
  mutate(id = factor(id), date = as.Date(date)) %>%
  filter(date < as.Date("2023-01-01")) %>%
  group_by(id) %>%
  filter(n() > 1) %>%                                  # repeated measures only
  mutate(months_since_first = as.numeric(date - min(date)) / (365.25 / 12)) %>%
  ungroup() %>%
  mutate(period = factor(ifelse(date < OUTBREAK_START, "Pre", "Post"),
                         levels = c("Pre", "Post"))) %>%
  filter(!is.na(vasexpartners_p6m)) %>%
  droplevels()

# ---- Descriptives -------------------------------------------------------------
by_month <- partners %>%
  group_by(date) %>%
  summarise(n = n(), median = median(vasexpartners_p6m),
            q1 = quantile(vasexpartners_p6m, 0.25),
            q3 = quantile(vasexpartners_p6m, 0.75), .groups = "drop")

# ---- Negative binomial GLMMs ------------------------------------------------
glmm_time   <- glmer.nb(vasexpartners_p6m ~ months_since_first + (1 | id), data = partners)
glmm_period <- glmer.nb(vasexpartners_p6m ~ period + (1 | id), data = partners)

tidy_glmm <- function(model, label) {
  tab <- coef(summary(model))
  ci  <- confint(model, parm = "beta_", method = "Wald")
  data.frame(model = label, term = rownames(tab),
             estimate = tab[, "Estimate"], std_error = tab[, "Std. Error"],
             z = tab[, "z value"], p_value = tab[, "Pr(>|z|)"],
             rate_ratio = exp(tab[, "Estimate"]),
             rr_lower = exp(ci[, 1]), rr_upper = exp(ci[, 2]),
             row.names = NULL)
}
glmm_table <- bind_rows(tidy_glmm(glmm_time,   "months_since_first"),
                        tidy_glmm(glmm_period, "pre_post_outbreak"))

write_table(by_month, "partners_by_month.csv")
write_table(glmm_table, "partners_glmm_nb.csv")

cat("\n==== Sexual partners: observations by survey month ====\n")
print(by_month, n = Inf)
cat("\nParticipants:", nlevels(partners$id), " Observations:", nrow(partners), "\n")
cat("\n==== Negative binomial GLMM: months since first visit ====\n")
print(summary(glmm_time))
cat("\n==== Negative binomial GLMM: pre vs post outbreak ====\n")
print(summary(glmm_period))
cat("\n==== Fixed effects (rate ratios with Wald 95% CI) ====\n")
print(glmm_table, digits = 3, row.names = FALSE)

# ---- Figures ----------------------------------------------------------------
suppressPackageStartupMessages(library(ggpubr))

log1p_partners_scale <- scale_y_continuous(
  trans  = "log1p",
  breaks = c(0, 1, 2, 5, 10, 25, 50, 100, 200, 500),
  labels = scales::label_number(accuracy = 1)
)
period_labels <- c("Pre" = "Pre-Mpox", "Post" = "Post-Mpox")
period_colours <- c("Pre" = "#4E79A7", "Post" = "#F28E2B")

# Figure 2D: partner numbers before vs after outbreak onset
n_labels <- partners %>%
  group_by(period) %>%
  summarise(n = n(), .groups = "drop")

fig2d <- ggplot(partners, aes(x = period, y = vasexpartners_p6m)) +
  geom_violin(aes(fill = period), alpha = 0.4, color = NA,
              trim = TRUE, width = 0.8) +
  geom_line(aes(group = id), color = "grey60", alpha = 0.4, linewidth = 0.4) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white",
               color = "black", alpha = 0.7) +
  geom_point(aes(color = period), size = 1.8, alpha = 0.6) +
  geom_text(data = n_labels,
            aes(x = period, label = paste0("n = ", n)),
            y = Inf, vjust = 1.5, size = 3.5, fontface = "bold",
            color = "grey30", inherit.aes = FALSE) +
  log1p_partners_scale +
  scale_x_discrete(labels = period_labels) +
  scale_fill_manual(values = period_colours) +
  scale_color_manual(values = period_colours) +
  labs(x = NULL, y = "Sexual partners (past 6 months)") +
  theme_pubr() +
  theme(legend.position = "none")

ggsave(file.path(FIGURE_DIR, "Fig2D_partners_pre_post.png"),
       fig2d, width = 5, height = 6, dpi = 300)

# Supplementary Figure 3: partner numbers over time with monthly median
fig_s3 <- ggplot(partners, aes(x = date, y = vasexpartners_p6m)) +
  geom_line(aes(group = id, color = period), alpha = 0.2, linewidth = 0.35) +
  geom_point(aes(color = period), size = 1.5, alpha = 0.5) +
  geom_line(data = by_month, aes(x = date, y = median),
            color = "black", linewidth = 1.1, inherit.aes = FALSE) +
  geom_point(data = by_month, aes(x = date, y = median, size = n),
             shape = 21, fill = "white", color = "black", stroke = 1.2,
             inherit.aes = FALSE) +
  geom_text(data = by_month, aes(x = date, y = median, label = paste0("n=", n)),
            size = 3, fontface = "bold", color = "black", inherit.aes = FALSE) +
  log1p_partners_scale +
  scale_size_continuous(range = c(10, 20), guide = "none") +
  scale_color_manual(values = period_colours, labels = period_labels, name = NULL) +
  labs(x = NULL, y = "Sexual partners (past 6 months)") +
  scale_x_date(date_breaks = "1 month", date_labels = "%b %Y") +
  theme_pubr() +
  theme(legend.position = "top",
        legend.text = element_text(size = 11),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(FIGURE_DIR, "FigS3_partners_over_time.png"),
       fig_s3, width = 10, height = 6, dpi = 300)
