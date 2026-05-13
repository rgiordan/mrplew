# mrplew Library Errata

Reviewed 2026-05-12 (initial); updated 2026-05-12 after user fixes.
No code was changed during this review.

Items marked **FIXED** were corrected by the user.
Items marked **NOTE: errata was wrong** were incorrect in the original review.

---

## Remaining Issues

### 1. `brms_lib.R:124` — `get_ols_mcmc_weights` crashes when `pop_frac=NULL`

The new code (after removing the `get_mrp_draws` call) is:

```r
mrp_draws <- yhat_pop_draws %*% pop_frac
```

The parameter `pop_frac` defaults to `NULL`, so this line throws an error when
`pop_frac` is not supplied.  The analogous `get_mrplew_logistic_brms` reaches
the NULL case via `get_mrp_draws_brms → get_mrp_draws`, which substitutes
`rep(1/npop, npop)` when `pop_frac=NULL`.  The OLS path needs the same
treatment, e.g. replace the raw matrix multiply with a call to
`get_mrp_draws(yhat_pop_draws, pop_frac)`.

---

### 2. `brms_lib.R:3,21` — Both docstring titles on the first two exported functions are identical

```r
#' Get mrplew weights for the brms logistic MCMC estimator.   ← get_mrp_draws_brms (line 3)
...
#' Get mrplew weights for the brms logistic MCMC estimator.   ← get_mrplew_logistic_brms (line 21)
```

The title on `get_mrp_draws_brms` should describe that function (e.g., "Get MrP
posterior draws from a brms model"), not the weight-computation step.

---

### 3. `balance_lib.R:55–58` — `sprintf` return silently discarded; no error is raised on collision

The variable name was fixed (`id_cols` → `id_col`), but the `if` block still
takes no action:

```r
if (id_col %in% c(names(df1), names(df2))) {
    sprintf("ID column %d is already present, which should never happen.", id_col)
}
```

The `sprintf` result is not passed to `stop()`, `warning()`, or `message()`.
A column-name collision would be silently ignored and the downstream row-split
would produce incorrect results.  Also note `%d` (integer format) is used for
the string `id_col`; should be `%s`.

---

### 4. `utils_lib.R:17` — Inline comment still says "ones" when default is `1/N`

The external-facing docstrings were updated to "Defaults to 1/N", but the
comment inside `get_population_frac` still reads:

```r
# Use pop_frac for weights if specified, otherwise use
# a vector of ones as long as pop_df.
```

The actual default is `rep(1, nrow(pop_df)) / nrow(pop_df)`, not ones.

---

### 5. `utils_lib.R:83–90` — `safe_get_eta_draws` looks for a field that is never stored

```r
safe_get_eta_draws <- function(mrplew_list, post, survey_df) {
    if ("eta_draws" %in% names(mrplew_list)) {
        eta_draws <- mrplew_list$eta_draws
    } else {
        eta_draws <- posterior_linpred(post, newdata=survey_df)
    }
    ...
}
```

No function in the library stores a field called `eta_draws` in its result list.
`get_mrplew_logistic_brms` saves `dloglikdy_survey_draws` (which equals eta for
logistic), not `eta_draws`. The fast-path branch will never be taken.

---

### 6. `simulation_lib.R:141` — "simualted" typo

```r
#' @return A list of simualted data
```

---

### 7. `weights_opt_lib.R` — Inconsistent `mrplew_w` scaling between OLS and logistic GLM

`get_mrplew_lm` multiplies the raw derivative by `nsur`:
```r
w_ols <- nsur * t(pop_frac) %*% x_pop %*% solve(xtx, t(x_ols))
```
so `sum(w_ols) ≈ nsur`.

`get_mrplew_logistic_glm` stores the raw derivative without this factor:
```r
w_logit <- t(mrp_chain_rule_term) %*% solve(hess, t(x_ols))
```
so `sum(w_logit) ≈ 1`.

Consequence: for a correctly specified saturated model, `w_ols ≈ w_opt` while
`w_logit ≈ w_opt / nsur`.  The MCMC function `get_mrplew_mcmc` uses the
`nsur`-multiplied convention (`mrplew_w <- n_obs * dmrp_dy`), consistent with
the OLS convention but inconsistent with the logistic GLM convention.
`compute_frequentist_sd` also assumes `sum(w) ≈ n_obs`.

The test in `test_opt.R` was updated to assert each method against its own
scale (`w_opt` for OLS, `w_opt / nsur` for logit with tolerance `1e-6`) to
reflect this inconsistency rather than mask it.

---

### 8. `brms_lib.R:115–116` — `re_formula` and `allow_new_levels` declared but never used

```r
get_ols_mcmc_weights <- function(brms_post, survey_df, pop_df, pop_frac=NULL,
                                 re_formula=NULL, allow_new_levels=FALSE)
```

Neither `re_formula` nor `allow_new_levels` appears anywhere in the function
body.  If these are intended for a future `posterior_epred` call they should be
documented; if not needed they should be removed to avoid confusing callers.

---

## Items Fixed by User

| # | Location | What was fixed |
|---|---|---|
| F1 | `binary_lib.R:46` | `runif(n_obs)` → `runif(length(y_orig))` |
| F2 | `balance_lib.R:80` | `colnames(x2) != colnames(x2)` → `colnames(x1) != colnames(x2)` |
| F3 | `balance_lib.R:94` | Duplicate `@param w1` → `@param w2` |
| F4 | `balance_lib.R:57` | `id_cols` → `id_col` |
| F5 | `binary_lib.R:32` | "expectatoins" → "expectations" |
| F6 | `brms_lib.R:125–128` | `get_mrplew_mcmc` now called with correct argument name `mrp_draws=` |
| F7 | `brms_lib.R:130` | Removed dangling `if (save_draws)` block (parameter was absent from signature) |
| F8 | `brms_lib.R:131` | Removed nonsensical `expit(dloglikdy)` for OLS (copy-paste from logistic) |
| F9 | `simulation_lib.R:138–139` | Duplicate `@param n_obs` → `@param n_obs_pop` |
| F10 | `brms_lib.R`, `mrplew_mcmc_lib.R` | Docstrings: "Defaults to ones" → "Defaults to 1/N" |
| F11 | `helper.R:44` | `w_opt=w_pop / count_sur` → `w_opt=frac_pop / frac_sur` (now matches `simulation_lib.R`) |
| F12 | `brms_lib.R:96` | Comment updated to match code: now correctly reads `-sigma^{-2}(y_n - ŷ_n)` |

---

## Correction to Original Errata

**Original errata item #4 was incorrect.**  The original review claimed the sign
in `get_gaussian_dloglikdy_draws` was wrong.  In fact, the correct derivative is:

`d/dy_i [-(y_i - ŷ_i)² / (2σ²)] = -(y_i - ŷ_i)/σ² = -resid/σ²`

The code (`-1 * resid_draws / sigma_draws^2`) was always correct.  The old
comment (before the user's fix) stated the positive sign, which was wrong.
The user's fix to the comment was appropriate.

---

## Test File Issues (Acknowledged as Not Yet Fixed)

| # | Location | Issue |
|---|---|---|
| T1 | `test_opt.R:31,34` | Calls `get_ols_weights` and `get_logit_weights`; correct names are `get_mrplew_lm` and `get_mrplew_logistic_glm` |
| T2 | `test_mcmc.R:28` | Calls `get_logit_mcmc_weights`; correct name is `get_mrplew_logistic_brms` |
| T3 | `test_opt.R:51`, `test_mcmc.R:48` | Accesses `$w`; correct field name is `$mrplew_w` |
| T4 | `test_mcmc.R:39` | `aggregate_simulation_data(sim_data, y_col)` — the exported library function takes only one argument; the two-argument version is local to `helper.R` |
