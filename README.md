# mrplew: Mister P Local Equivalent Weights

**This library is currently experimental and not for public use.**

**Please don't use it unless you know and are in contact with the authors!**

## What is mrplew?

Multilevel Regression with Poststratification (MrP) is a widely-used method
for estimating population quantities from non-representative survey data. A MrP
estimate is typically computed by fitting a regression model to survey responses,
using that model to predict outcomes for every cell of the population, and then
averaging those predictions weighted by the population cell sizes.

`mrplew` computes **local equivalent weights** (LEWs) for MrP estimates. A LEW
re-expresses the MrP estimate as a simple weighted average of the individual
survey outcomes — one weight per respondent. This is useful because:

- **Interpretability**: MrP weights behave like traditional survey weights, and
  can be compared directly to raking or post-stratification weights.
- **Diagnostics**: Covariate balance can be checked using the same tools as for
  weighted surveys.
- **Uncertainty**: The weighted-average representation makes variance estimation
  straightforward.

`mrplew` supports OLS (`lm`), logistic regression (`glm`), and full Bayesian
MCMC models fitted with `brms`.

## Installation

```r
library(devtools)
install_github("https://github.com/rgiordan/mrplew", upgrade = "never", force = TRUE)
```

## Quick start: MrPlew weights from a brms logistic model

Here is a minimal end-to-end example. It uses the built-in simulation helper,
fits a logistic regression in `brms`, and then computes MrPlew weights.

```r
library(tidyverse)
library(brms)
library(mrplew)

# ------------------------------------------------------------------
# 1. Simulate survey and population data
# ------------------------------------------------------------------
set.seed(42)
sim_data <- mrplew::simulate_survey_data(
  n_groups   = 4,
  n_obs      = 2000,
  n_obs_pop  = 5000,
  degree     = 2
)

survey_df <- sim_data$survey_df   # individual survey respondents
pop_df    <- sim_data$pop_df      # population frame

# Aggregate population data (speeds up posterior prediction)
agg_list    <- aggregate_simuilation_data(sim_data)
pop_agg_df  <- agg_list$pop_agg_df

# ------------------------------------------------------------------
# 2. Fit a logistic regression with brms
# ------------------------------------------------------------------
g_sum    <- paste(sim_data$group_effects$g_cols, collapse = " + ")
reg_form <- formula(sprintf("y ~ 1 + (%s)^2", g_sum))

brms_fit <- brm(
  reg_form,
  data    = survey_df,
  family  = bernoulli(link = "logit"),
  chains  = 4,
  cores   = 4,
  warmup  = 500,
  iter    = 2000,
  seed    = 1234,
  file    = "logit_posterior"   # cache the fit
)

# ------------------------------------------------------------------
# 3. Compute MrPlew weights
# ------------------------------------------------------------------
mrplew_result <- get_mrplew_logistic_brms(
  brms_post  = brms_fit,
  survey_df  = survey_df,
  pop_df     = pop_agg_df,
  pop_w      = pop_agg_df$w
)

# The MrP estimate (posterior mean)
cat("MrP estimate:", mean(mrplew_result$mrp_draws), "\n")

# One weight per survey respondent — use these just like survey weights
weights <- mrplew_result$mrplew_w
cat("Weighted survey mean:", weighted.mean(survey_df$y, weights), "\n")
```

## Checking covariate balance

Once you have weights, you can check whether they rebalance the survey with
respect to covariates that were not included in the regression model. A
well-specified model should produce near-zero imbalance.

```r
# Use a richer formula to check balance on higher-order terms
balance_form <- formula(sprintf("y ~ 1 + (%s)^4", g_sum))

balance_result <- check_covariate_balance(
  mrplew_list = mrplew_result,
  survey_df   = survey_df,
  pop_df      = pop_agg_df,
  reg_form    = balance_form,
  pop_w       = pop_agg_df$w
)

# Tidy balance summary
balance_result$balance_df
```

The `balance_df` table reports the weighted population mean, the MrPlew-weighted
survey mean, and the difference for each covariate term. Differences close to
zero indicate good balance.

## OLS and frequentist logistic weights

`mrplew` also provides closed-form weights for frequentist models:

```r
# OLS
lm_fit     <- lm(reg_form, survey_df)
ols_result <- get_mrplew_lm(lm_fit, survey_df, pop_df)
ols_result$mrp          # MrP point estimate
ols_result$mrplew_w     # one weight per respondent

# Logistic regression (glm)
glm_fit      <- glm(reg_form, survey_df, family = binomial(link = "logit"))
glm_result   <- get_mrplew_logistic_glm(glm_fit, survey_df, pop_df)
glm_result$mrp
glm_result$mrplew_w
```

## More

A fully worked example, including balance plots and weight diagnostics, is in
`examples/example.R`.
