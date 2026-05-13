# mrplew Library Errata

Reviewed 2026-05-12. All issues are read-only findings; no code was changed.

---

## Critical Bugs

### 1. `binary_lib.R:46` — `n_obs` undefined in `draw_conditional_binary_data`

When `unif_base` is `NULL`, the code calls `runif(n_obs)` but `n_obs` is not a
parameter or local variable of the function. Should be `runif(length(y_orig))`.

---

### 2. `brms_lib.R:125–128` — `get_ols_mcmc_weights` calls `get_mrplew_mcmc` with wrong argument names

```r
result_list <- get_mrplew_mcmc(
    yhat_pop_draws=yhat_pop_draws,        # ← wrong name
    dloglikdy_survey_draws=dloglikdy_survey_draws,
    pop_frac=pop_frac)                    # ← not a parameter of get_mrplew_mcmc
```

`get_mrplew_mcmc` (in `mrplew_mcmc_lib.R`) takes `mrp_draws` and
`dloglikdy_survey_draws` only; it does not accept `yhat_pop_draws` or
`pop_frac`. The correct call should first convert `yhat_pop_draws` to `mrp_draws`
via `get_mrp_draws(yhat_pop_draws, pop_frac)`, then pass the result as
`mrp_draws`. Compare with the correct pattern in `get_mrplew_logistic_brms`.

---

### 3. `brms_lib.R:130` — `save_draws` is not a parameter of `get_ols_mcmc_weights`

The function body references `if (save_draws)` but the function signature is

```r
get_ols_mcmc_weights <- function(brms_post, survey_df, pop_df, pop_frac=NULL,
                                 re_formula=NULL, allow_new_levels=FALSE)
```

`save_draws` is absent. This would throw "object 'save_draws' not found"
whenever the `if` branch is reached. The analogous `get_mrplew_logistic_brms`
correctly declares `save_draws=FALSE`.

---

### 4. `brms_lib.R:131` — `expit(dloglikdy_survey_draws)` is wrong for OLS

Inside the `if (save_draws)` block of `get_ols_mcmc_weights`:

```r
yhat_survey_draws <- expit(dloglikdy_survey_draws)
```

For OLS, `dloglikdy_survey_draws` is `-resid / sigma^2` (i.e., a scaled
residual), not a log-odds. Applying `expit` to it produces meaningless values.
This is a copy-paste error from the logistic version, where `dloglikdy` really
is the log-odds `eta` and `expit(eta) = yhat` is meaningful.

---

### 5. `balance_lib.R:80` — `check_balance_matrices` always passes the column-name check

```r
if (any(colnames(x2) != colnames(x2))) {   # compares x2 to itself
```

Both sides are `x2`; the condition is identically `FALSE`. The check never
fires regardless of whether `x1` and `x2` have matching column names. Should be
`colnames(x1) != colnames(x2)`.

---

### 6. `mrplew_mcmc_lib.R:97` — Wrong sign in Gaussian log-likelihood derivative

```r
dloglikdy_survey_draws <- -1 * resid_draws / (sigma_draws^2)
```

The Gaussian log-likelihood is `-(y - yhat)^2 / (2 sigma^2)`. Its derivative
with respect to `y_n` is `+(y_n - yhat_n) / sigma^2 = resid_n / sigma^2`. The
code negates this. The immediately preceding comment (line 96) correctly states
the expected sign:

```
# The log likelihood derivative for the n^th datapoint is
# sigma^{-2} (y_n - \hat{y}_n)
```

but the code contradicts it. The logistic version (`posterior_linpred`, line 54)
uses the correct positive sign. This sign error propagates into every OLS weight
computation via the MCMC path.

---

## Test File Bugs

### 7. `test_opt.R:31,34` — Nonexistent function names

The test calls `get_ols_weights` and `get_logit_weights`, neither of which
exists. The correct names are `get_mrplew_lm` and `get_mrplew_logistic_glm`.

---

### 8. `test_mcmc.R:28` — Nonexistent function `get_logit_mcmc_weights`

The test references `get_logit_mcmc_weights`; the actual function is
`get_mrplew_logistic_brms`.

---

### 9. `test_opt.R:51` and `test_mcmc.R:48` — Wrong field name `$w`

Both tests access `mrp_weights$w` / `mcmc_mrp$w`, but all weight-returning
functions store the weights under `$mrplew_w`.

---

### 10. `test_mcmc.R:39` — `aggregate_simulation_data` called with two arguments

```r
agg_list <- aggregate_simulation_data(sim_data, y_col)
```

The package's exported `aggregate_simulation_data` (in `simulation_lib.R`)
accepts only `sim_data`. The two-argument version lives in `tests/testthat/helper.R`
(local only) and would shadow the exported one inside the test. However, the
`helper.R` version calls its second parameter `resp` and the call passes
`y_col`. If `y_col` is a string, this happens to work, but the variable name
mismatch is confusing and the functions are not the same.

---

## Weight Scaling Inconsistencies

### 11. `simulation_lib.R:204` vs `tests/testthat/helper.R:44` — `w_opt` differs by a factor of `nsur`

**`simulation_lib.R` (exported `aggregate_simulation_data`):**
```r
mutate(w_opt = frac_pop / frac_sur)
```
When joined to individual survey rows and summed: each row in group *s* gets
weight `frac_pop_s / frac_sur_s = frac_pop_s * nsur / count_sur_s`.
Sum over all survey rows = `nsur * sum(frac_pop) = nsur`.

**`helper.R` (local test version):**
```r
mutate(w_opt = w_pop / count_sur)
```
When joined to individual survey rows and summed: each row in group *s* gets
weight `frac_pop_s / count_sur_s`. Sum = `sum(frac_pop) = 1`.

These two formulas differ by a factor of `nsur`. The test in `test_opt.R`
compares `mrp_weights$mrplew_w` (which uses the library's convention and sums to
`nsur`) against `w_opt` from the helper version (which sums to 1), so the
assertion would fail by a factor of `nsur` even if the function names were
correct.

---

### 12. Weight-sum convention: `mrplew_w` sums to approximately `n_obs`

`get_mrplew_lm` returns `w_ols = nsur * (pop_frac^T X_pop) (X_sur^T X_sur)^{-1} X_sur^T`,
which sums to `nsur` when the model contains an intercept.

`get_mrplew_mcmc` returns `mrplew_w = n_obs * cov(mrp_draws, dloglikdy)`,
where `n_obs = ncol(dloglikdy_survey_draws)` = number of survey observations.

The `compute_frequentist_sd` docstring (line 131 of `balance_lib.R`) explicitly
states weights should "sum to the number of effective observations," consistent
with this convention.

However, `check_covariate_balance` normalises the survey weights by `nsur`
before passing them to `get_balance_df`, which expects weights summing to ~1.
For the OLS case this is mathematically correct. For the MCMC case, the sum of
`n_obs * cov(mrp, dloglik)` is `n_obs * cov(mrp, sum_n dloglik_n)`, which does
not generally equal `n_obs`, so the renormalised weights in the balance check
may not sum to 1 as expected.

---

## Documentation / Docstring Errors

### 13. `balance_lib.R:94` — Duplicate `@param w1`, missing `@param w2`

```r
#' @param w1 A vector of normalized weights for x1, summing to ~ 1
#' @param w1 A vector of normalized weights for x2, summing to ~ 1   ← should be @param w2
```

---

### 14. `utils_lib.R:17` — Comment says "ones" but code defaults to `1/N`

```r
# Use pop_frac for weights if specified, otherwise use
# a vector of ones as long as pop_df.
```
The actual default is `rep(1, nrow(pop_df)) / nrow(pop_df)` — a vector of
`1/N`, not ones. A vector of ones would fail the `abs(weight_sum - 1) > 1e-6`
check on the next line.

The same "Defaults to ones" phrasing appears in the `@param pop_frac`
docstrings of `get_mrp_draws_brms`, `get_mrplew_logistic_brms`, and
`get_ols_mcmc_weights`. All should say "Defaults to `1/N`" or "uniform weights."

`get_mrp_draws` (line 18–21 of `mrplew_mcmc_lib.R`) also defaults to `1/npop`,
not ones.

---

### 15. `brms_lib.R:21` — Wrong docstring title on `get_mrplew_logistic_brms`

The docstring reads "Get mrplew weights for the brms logistic MCMC estimator"
— identical to `get_mrp_draws_brms` immediately above it (line 3). It was
apparently copied without updating the title.

---

### 16. `utils_lib.R:83–90` — `safe_get_eta_draws` looks for `eta_draws` that is never stored

`safe_get_eta_draws` checks for `"eta_draws"` in `mrplew_list`, but no function
in the library stores a field named `eta_draws`. The relevant quantity stored by
`get_mrplew_logistic_brms` when `save_draws=TRUE` is `dloglikdy_survey_draws`.
The function will always fall through to calling `posterior_linpred`.

---

## Typos

### 17. `balance_lib.R:57` — `id_cols` should be `id_col`

```r
sprintf("ID column %d is already present, which should never happen.", id_cols)
```
The variable is named `id_col` (singular); `id_cols` is undefined. Additionally,
the format specifier `%d` (integer) is used for a string column name, and the
`sprintf` result is silently discarded (no `warning()` or `stop()`), so the
message is never emitted.

---

### 18. `binary_lib.R:32` — "expectatoins" should be "expectations"

```r
#' @param e_y A vector of estimated expectatoins of y_orig
```

---

### 19. `simulation_lib.R:138,139` — Duplicate `@param n_obs`

```r
#' @param n_obs The number of survey observations
#' @param n_obs The number of population observations
```
The second entry should be `@param n_obs_pop`.

---

### 20. `simulation_lib.R:142` — "simualted" should be "simulated"

```r
#' @return A list of simualted data
```

---

## Minor Scoping / Logic Issues

### 21. `balance_lib.R:53–58` — `get_consistent_regressors` silently swallows collision error

```r
if (id_col %in% c(names(df1), names(df2))) {
  sprintf("ID column %d is already present, which should never happen.", id_cols)
}
```
There is no `stop()` or `warning()` call; the sprintf result is discarded. If
the randomly-generated `id_col` happens to collide with an existing column name,
execution continues silently and the downstream splitting logic will produce
incorrect results.

---

### 22. `brms_lib.R:43–44` — `is.null(mrp_draws) & is.null(pop_df)` uses scalar `&` in an `if`

In R `if` conditions, the recommended operator is `&&` (short-circuit). Using
`&` works for scalar logicals but is non-idiomatic and will trigger a warning if
either side is ever a vector. This pattern appears in one place.

---

### 23. `test_mrpaw.R:12–13` — `context()` is deprecated in testthat 3

Both test files use `context("mrplew")`, which is deprecated since testthat 3.x
and does nothing. The file name serves as context automatically.
