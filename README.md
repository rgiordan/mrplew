# mrplew: Mister P Local Equivalent Weights

This library is currently a public port of research code!  Please
use with care, and contact the authors (or post a Github issue)
if you have doubts or questions.


## What is mrplew?

This software package accompanies our paper, [Locally Equivalent Weights for
Multilevel Regression and
Poststratification](https://rgiordan.github.io/assets/mrplew_paper.pdf).

MrPlew stands for "multilevel poststratification and regression local
equivalent weights."  Please see our paper for more details and motivation. 

The package `mrplew` supports OLS (`lm`), logistic regression (`glm`), and
Gaussian or logistic Bayesian MCMC models fitted with `brms`.

## Installation

```r
library(devtools)
install_github("https://github.com/rgiordan/mrplew", force = TRUE)
```

## Quick start: MrPlew weights from a brms logistic model

An example can be found in `examples/example.R`.  Here are some excerpts.

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
  n_groups  = 4,
  n_obs     = 2000,
  n_obs_pop = 5000,
  degree    = 2
)

survey_df <- sim_data$survey_df   # individual survey respondents
pop_df    <- sim_data$pop_df      # population frame

# Aggregate population data by group (speeds up posterior prediction)
agg_list   <- mrplew::aggregate_simulation_data(sim_data)
pop_agg_df <- agg_list$pop_agg_df   # pop_agg_df$frac gives population fractions

# ------------------------------------------------------------------
# 2. Fit a logistic regression with brms
# ------------------------------------------------------------------
g_sum    <- paste(sim_data$group_effects$g_cols, collapse = " + ")
reg_form <- formula(sprintf("y ~ 1 + (%s)^2", g_sum))

brms_fit <- brm(
  reg_form,
  data   = survey_df,
  family = bernoulli(link = "logit"),
  chains = 4,
  cores  = 4,
  warmup = 500,
  iter   = 2000,
  seed   = 1234,
  file   = "logit_posterior"   # cache the fit
)

# ------------------------------------------------------------------
# 3. Compute MrP draws, then MrPlew weights
# ------------------------------------------------------------------

# First get posterior draws of the MrP estimate.
# pop_frac gives the fraction of the population in each aggregated cell.
mrp_draws <- mrplew::get_mrp_draws_brms(
  brms_post = brms_fit,
  pop_df    = pop_agg_df,
  pop_frac  = pop_agg_df$frac
)

cat("MrP estimate:", mean(mrp_draws), "\n")

# Then compute the local equivalent weights using those draws.
mrplew_result <- mrplew::get_mrplew_logistic_brms(
  brms_post = brms_fit,
  survey_df = survey_df,
  mrp_draws = mrp_draws
)

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

balance_result <- mrplew::check_covariate_balance(
  mrplew_list = mrplew_result,
  survey_df   = survey_df,
  pop_df      = pop_agg_df,
  reg_form    = balance_form,
  pop_frac    = pop_agg_df$frac
)

# Tidy balance summary
balance_result$balance_df
```

The `balance_df` table reports the weighted population mean, the MrPlew-weighted
survey mean, and the difference for each covariate term. Differences close to
zero indicate good balance.

## Variance estimation

MrPlew weights make frequentist variance estimation simple. Given the model
residuals at the survey observations, `compute_frequentist_sd` returns a
standard deviation for the MrP estimate:

```r
# Get posterior mean predictions at the survey points
yhat_draws <- brms::posterior_epred(brms_fit, newdata = survey_df)
yhat       <- colMeans(yhat_draws)
resid      <- survey_df$y - yhat

# OLS weights also work here
lm_fit     <- lm(reg_form, survey_df)
ols_result <- mrplew::get_mrplew_lm(lm_fit, survey_df, pop_df)

freq_sd <- mrplew::compute_frequentist_sd(w = ols_result$mrplew_w, resid = resid)
cat("Frequentist SD:", freq_sd, "\n")

# For comparison, the posterior SD from the MCMC draws:
cat("Posterior SD:  ", sd(mrp_draws), "\n")
```

## OLS and frequentist logistic weights

`mrplew` also provides closed-form weights for frequentist models:

```r
# OLS
lm_fit     <- lm(reg_form, survey_df)
ols_result <- mrplew::get_mrplew_lm(lm_fit, survey_df, pop_df)
ols_result$mrp       # MrP point estimate
ols_result$mrplew_w  # one weight per respondent

# Logistic regression (glm)
glm_fit    <- glm(reg_form, survey_df, family = binomial(link = "logit"))
glm_result <- mrplew::get_mrplew_logistic_glm(glm_fit, survey_df, pop_df)
glm_result$mrp
glm_result$mrplew_w
```

## More

A fully worked example, including balance plots and weight diagnostics, is in
`examples/example.R`.
