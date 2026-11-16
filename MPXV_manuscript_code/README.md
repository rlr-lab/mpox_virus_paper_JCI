# Analyses of the 2022 mpox outbreak

Code for the analyses in the manuscript. The code covers three sets of
analyses:

1. **Sexual behavior in the RADAR cohort.** Changes in self-reported sexual
   partner numbers during 2022 (R).
2. **Epidemiological statistical modeling.** Associations between mpox case
   counts, vaccination and behavioral metrics in Chicago, Los Angeles and
   New York (R). Modeling summaries are given in Supplemental Figure 4.
3. **Phylogenetic analysis of hMPXV genomes.** Genome assembly, maximum
   likelihood and time-scaled phylogenies, and ancestral reconstruction of
   geographic location (shell commands). See
   [`phylogenetics/README.md`](phylogenetics/README.md).

For the statistical analyses (sections 1 and 2), statistical significance was
defined as a two-sided p-value < 0.05.

```
R/                 statistical analyses (sections 1 and 2)
phylogenetics/     genome assembly and phylogenetic analyses (section 3)
data/              input data for the R analyses (not distributed)
```

## 1. Sexual partner numbers

`R/04_sexual_partners.R`

- **Data:** RADAR cohort participants with more than one visit in 2022.
- **Outcome:** self-reported number of sexual partners in the past 6 months.
- **Longitudinal change:** negative binomial generalized linear mixed-effects
  models with participant-level random intercepts (`lme4::glmer.nb`). Time
  enters each model in one of two ways:
  - months since the participant's first 2022 visit;
  - before vs after the start of the mpox outbreak (4 June 2022).

## 2. Epidemiological statistical modeling

For every dataset:

1. The outcome is `log(MPXVcases + 1)`, or the raw case counts for negative
   binomial models. First and second JYNNEOS vaccine doses are entered at a
   time lag.
2. Residual autocorrelation is tested with the Durbin–Watson test on the OLS
   fit. It was significant in every dataset, so correlation structures are
   added to the models.
3. Three model types are fit:
   - generalized least squares with an autoregressive AR(1) correlation
     structure (GLS AR(1));
   - ordinary least squares (OLS);
   - negative binomial regression (NB).

   Each type is fit at several vaccine lags (temporal structures).
4. Models are compared with the Akaike (AIC) and Bayesian (BIC) information
   criteria, and the model with the lowest AIC is retained.

GLS AR(1) had the best performance in all analyses.

| Script | Analysis | Predictors | Lags compared | Final model |
|---|---|---|---|---|
| `R/01_chicago_monthly_gls.R` | Chicago monthly epidemiology model, Jun 2022 – Mar 2023 | Dose 1, Dose 2, total diagnostic tests, STI positivity rate, unprotected/protected sex acts ratio (combinations of variables compared) | 1 and 2 months | GLS AR(1) with all variables; 1- and 2-month lags reported as complementary analyses |
| `R/02_chicago_weekly_gls.R` | Chicago weekly (granular) model | Dose 1, Dose 2 | 1–4 weeks | GLS AR(1), 4-week lag |
| `R/03_three_cities_gls.R` | Chicago, Los Angeles and New York weekly model, fit for each city | Dose 1, Dose 2 | 1–4 weeks | GLS AR(1), 4-week lag (all three cities) |

`R/00_functions.R` contains the shared helper functions.

## 3. Phylogenetic analysis of hMPXV genomes

Shell scripts in `phylogenetics/`, run in order. Software versions, the
command for each step and the Gblocks masking (done interactively in SeaView)
are described in [`phylogenetics/README.md`](phylogenetics/README.md).

| Script | Steps |
|---|---|
| `01_consensus.sh` | Two-round reference mapping of Nanopore reads (minimap2) and consensus calling (samtools) |
| `02_alignment_ml_trees.sh` | MAFFT alignment to NC_063383.1; IQ-TREE maximum likelihood trees (trimmed and untrimmed alignments) |
| `03_contextual_dataset.sh` | Nextclade QC of 831 public US genomes; joint alignment and IQ-TREE tree |
| `04_timetree_mugration.sh` | TreeTime time-scaled phylogeny and mugration reconstruction of US Census region |

## Data (R analyses)

Data are not included in this repository. To run the code, place the
following files in `data/` with the column names shown.

| File | Contents | Sources |
|---|---|---|
| `data/RADAR_STI_2022_NonUnique_Partners.csv` | One row per participant visit: `id`, `date`, `vasexpartners_p6m` | RADAR cohort |
| `data/Chi_EpiData_2022_23.csv` | Chicago monthly: `Month`, `MPXVcases`, `Dose1`, `Dose2`, `Total_Tests`, `STIPositivityRate`, `unprotected_protected` | Cases and JYNNEOS doses: Chicago Department of Public Health (CDPH). Total diagnostic tests in the MSM population (chlamydia, gonorrhea, syphilis, HIV and mpox): Howard Brown Health. STI positivity rate and sex acts ratio: RADAR cohort (summarized per month after summing all events for each person and month) |
| `data/AllCDPHData_FirstWave_Granular.csv` | Chicago weekly: `date`, `MPXVcases`, `Dose1`, `Dose2` | CDPH |
| `data/AllCities_FirstWave.csv` | Weekly: `Week`, `MPXVcases`, `Dose1`, `Dose2`, `City` | Chicago: CDPH. New York: NYC Department of Health and Mental Hygiene GitHub data repository. Los Angeles: Los Angeles County Department of Public Health |

## Running the R analyses

The manuscript analyses used R v4.4.1 with nlme v3.1-168, MASS v7.3-65,
ggplot2 v4.0.0, dplyr v1.1.4 and car v3.1-3, plus lme4 for the mixed models
and ggpubr for the sexual partner figures.

```r
install.packages(c("nlme", "MASS", "lme4", "dplyr", "car", "ggplot2", "ggpubr"))
```

From the repository root:

```sh
Rscript run_all.R
```

Each script can also be run on its own, e.g. `Rscript R/04_sexual_partners.R`.
The package versions of each run are written to `results/sessionInfo.txt`.

## Outputs of the R analyses

`results/tables/`

| File | Contents |
|---|---|
| `partners_by_month.csv` | Number of responses, median and IQR of partner numbers per survey month |
| `partners_glmm_nb.csv` | Negative binomial GLMM fixed effects, with rate ratios and Wald 95% CIs |
| `chicago_monthly_model_comparison.csv` | AIC, BIC, AR(1) φ, dose p-values and Durbin–Watson p-value for each model type × variable combination × lag |
| `chicago_monthly_final_model_coefficients.csv` | GLS AR(1) all-variable model coefficients at 1- and 2-month lags |
| `chicago_weekly_model_comparison.csv` | Model type × lag comparison |
| `chicago_weekly_final_model_coefficients.csv` | Final model coefficients |
| `three_cities_model_comparison.csv` | City × model type × lag comparison |
| `three_cities_final_model_coefficients.csv` | Final model coefficients for each city |

`results/figures/` (sexual partner analysis only)

- `Fig2D_partners_pre_post.png`: Figure 2D, partner numbers before vs after
  outbreak onset.
- `FigS3_partners_over_time.png`: Supplementary Figure 3, partner numbers per
  participant over time, with the monthly median and number of responses.

## Notes on the R analyses

- Lagging drops the first *k* observations of each series. The AR(1) time
  index is assigned before lagging.
- Durbin–Watson p-values from `car::durbinWatsonTest` are bootstrap-based
  (1,000 replicates), so a value of 0 means p < 0.001. A seed is set for
  reproducibility.
- Negative binomial models that fail to converge (some monthly models with
  n = 8–9) are reported with `converged = FALSE`.
- Negative binomial AIC is on the count scale, while GLS and OLS AIC are on
  the log scale.

## Data availability

Whole-genome sequences generated in this study are available through NCBI
GenBank under accession numbers PZ534571–PZ534651.
