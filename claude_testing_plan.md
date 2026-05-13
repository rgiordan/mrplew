# Testing update plan

This document records the inconsistencies between the current test suite and
the current source code.  It is a guide for a future editing session; no code
changes should be made based on this file alone.

---

## Overview of what changed in the source

The key API changes that broke the tests are:

| Old name / parameter | New name / parameter | Location |
|---|---|---|
| `get_ols_weights` | `get_mrplew_lm` | `weights_opt_lib.R` |
| `get_logit_weights` | `get_mrplew_logistic_glm` | `weights_opt_lib.R` |
| `get_logit_mcmc_weights` | `get_mrplew_logistic_brms` (see below) | `brms_lib.R` |
| `mrp_weights$w` | `mrp_weights$mrplew_w` | all return lists |
| `pop_w=` | `pop_frac=` | `get_mrp_draws`, `check_covariate_balance` |
| `aggregate_simuilation_data` (typo) | `aggregate_simulation_data` | `simulation_lib.R` |
| single-step logit MCMC weights | two-step: `get_mrp_draws_brms` then `get_mrplew_logistic_brms` | `brms_lib.R` |

Additionally, the new function `compute_frequentist_sd` in `balance_lib.R` is
not yet tested.

---

## File-by-file issues

### `tests/testthat/test_opt.R`

**Issue 1 — wrong function names (lines 31, 33)**

The test calls `get_ols_weights` and `get_logit_weights`, neither of which
exists.  These should be `get_mrplew_lm` and `get_mrplew_logistic_glm`
respectively.

```r
# Current (broken)
mrp_weights <- get_ols_weights(fit, sim_data$survey_df, sim_data$pop_df)
mrp_weights <- get_logit_weights(fit, sim_data$survey_df, sim_data$pop_df)

# Should be
mrp_weights <- get_mrplew_lm(fit, sim_data$survey_df, sim_data$pop_df)
mrp_weights <- get_mrplew_logistic_glm(fit, sim_data$survey_df, sim_data$pop_df)
```

**Issue 2 — wrong return field name (line 51)**

The test checks `mrp_weights$w`, but the return lists from both functions use
the field `mrplew_w`.

```r
# Current (broken)
AssertNearlyEqual(mrp_weights$w, w_opt)

# Should be
AssertNearlyEqual(mrp_weights$mrplew_w, w_opt)
```

**Issue 3 — the `logit` branch of the loop sets `fit` to a `glm` result**

After fixing Issues 1 and 2, verify that the assertion on line 41 (`predict`
+ `AssertNearlyEqual`) is still logically correct for both the OLS and logit
iterations of the loop.  The `predict(..., type="response")` call works for
both `lm` and `glm`, so the logic is probably fine, but it should be checked.

---

### `tests/testthat/test_mcmc.R`

**Issue 1 — `get_logit_mcmc_weights` no longer exists (lines 29, 47)**

The test dispatches to `get_logit_mcmc_weights` but this function has been
removed.  The replacement is a two-step call:

```r
# Step 1: get MrP posterior draws
mrp_draws <- get_mrp_draws_brms(
  brms_post = post,
  pop_df    = agg_list$pop_agg_df,
  pop_frac  = agg_list$pop_agg_df$frac)

# Step 2: get weights from those draws
mcmc_mrp <- get_mrplew_logistic_brms(
  brms_post = post,
  survey_df = sim_data$survey_df,
  mrp_draws = mrp_draws)
```

Because the OLS and logit paths now have different call signatures, the
`mrplewFunction` dispatch pattern (single function variable called with the
same arguments for both methods) no longer works.  The test should be
refactored so that OLS and logit have separate call sites — either two
separate `test_that` blocks or explicit `if`/`else` branches that do the full
call for each family.

**Issue 2 — wrong column name in `pop_frac` argument (line 46)**

The test currently passes `pop_frac=agg_list$pop_agg_df$w`.  The helper's
`aggregate_simulation_data` (see helper.R section below) creates a column
named `frac`, not `w`, so `$w` evaluates to `NULL`.  This silently disables
population weighting rather than raising an error.

```r
# Current (broken — $w is NULL)
pop_frac=agg_list$pop_agg_df$w

# Should be
pop_frac=agg_list$pop_agg_df$frac
```

**Issue 3 — wrong return field name (line 48)**

```r
# Current (broken)
expect_true(length(mcmc_mrp$w) == nrow(sim_data$survey_df))

# Should be
expect_true(length(mcmc_mrp$mrplew_w) == nrow(sim_data$survey_df))
```

**Issue 4 — stale `agg_list` load from RDS (line 35)**

Line 35 reads `agg_list <- rds_load$agg_list`, but `RunAndCachePosterior`
(helper.R line 87–90) never saves `agg_list` in the RDS.  The result is
always `NULL`, and the variable is immediately overwritten on line 39.  This
is harmless but confusing.  The dead assignment on line 35 should be removed.

---

### `tests/testthat/helper.R`

**Issue 1 — local `aggregate_simulation_data` shadows the exported package function**

The helper defines its own `aggregate_simulation_data(sim_data, resp="y")`
on line 17.  The package now exports a function with the same name but only
one parameter.  The helper's two-argument version is still needed by
`test_mcmc.R` (which passes `y_col` as the second argument).

Options (decide before editing):
- Keep the helper's private version but rename it (e.g., `aggregate_simulation_data_with_resp`) to avoid confusion.
- Update the exported `aggregate_simulation_data` to accept an optional `resp` argument, and remove the helper's private version.

Either way, all call sites in the tests must be updated consistently.

**Issue 2 — `frac` vs `w` column naming**

The helper's `aggregate_simulation_data` returns data frames with a `frac`
column (line 29, 35).  The test_mcmc.R call at line 46 uses `$w` (see test_mcmc
Issue 2 above).  The column name in the helper is already correct (`frac`);
the consumer is wrong.

---

### Source code issues in `R/brms_lib.R` (will cause test failures even after test fixes)

These are bugs in the library source, not the tests, but they will prevent
tests from passing and must be fixed alongside the test updates.

**Source issue 1 — `get_mrp_draws_brms` passes `pop_w` to `get_mrp_draws` which now expects `pop_frac`**

`get_mrp_draws_brms` (line 14–19 of `brms_lib.R`) calls:
```r
mrp_draws <- get_mrp_draws(yhat_pop_draws=yhat_pop_draws, pop_w=pop_w)
```

But `get_mrplew_mcmc_lib.R`'s `get_mrp_draws` now uses the parameter `pop_frac`,
not `pop_w`.  The argument name must be updated in `get_mrp_draws_brms`.

Also, the `example.R` passes `pop_frac=pop_agg_df$frac` as a named argument to
`get_mrp_draws_brms`, but the function signature only has `pop_w=NULL`; the
`pop_frac` argument is silently caught by `...` and forwarded to
`posterior_epred`, which ignores it.  The `get_mrp_draws_brms` signature should
be updated to `pop_frac=NULL` and the internal call updated accordingly.

**Source issue 2 — `get_ols_mcmc_weights` calls `get_mrplew_mcmc` with wrong arguments**

`get_ols_mcmc_weights` (lines 115–138 of `brms_lib.R`) calls:
```r
result_list <- get_mrplew_mcmc(
    yhat_pop_draws=yhat_pop_draws,
    dloglikdy_survey_draws=dloglikdy_survey_draws,
    pop_w=pop_w)
```

But `get_mrplew_mcmc` now takes `(mrp_draws, dloglikdy_survey_draws)` only.
The `yhat_pop_draws` must first be converted to `mrp_draws` via `get_mrp_draws`,
and `pop_w` must be renamed to `pop_frac` throughout.

**Source issue 3 — `get_ols_mcmc_weights` references undefined `save_draws`**

The `if (save_draws)` block near the end of `get_ols_mcmc_weights` references
a variable `save_draws` that is not in the function's parameter list.  The
parameter `save_draws=FALSE` must be added to the function signature (matching
`get_mrplew_logistic_brms`).

---

### `tests/testthat/test_mrpaw.R` — no changes required

`simulate_survey_data` still has the same exported signature and return
structure.  The test builds its own aggregation without calling the package's
`aggregate_simulation_data`.  The helper functions `AssertNearlyEqual` and
`AssertNearlyZero` are unchanged.  The two `TODO` stubs (`covariate_balance`
and `influence_weighting`) remain empty and still pass trivially.

---

### `tests/testthat/test_utils.R` — no changes required

`get_block_bootstrap_covariance_draws` has the same signature and return
structure.  No issues found.

---

## New test coverage needed

`compute_frequentist_sd` (`balance_lib.R`, lines 131–143) is exported and
used in `example.R` but has no test.  A new `test_that` block should be
added to test this function — for example, verifying that it returns zero
when residuals are zero, or checking the result against a known analytical
value.

---

## Suggested edit order

1. Fix source bugs in `R/brms_lib.R` (issues 1–3 above) before touching tests,
   so that running tests after each fix gives meaningful signal.
2. Update `tests/testthat/helper.R` — resolve the `aggregate_simulation_data`
   shadowing question and pick a consistent column name strategy.
3. Fix `tests/testthat/test_opt.R` (Issues 1–2).
4. Fix `tests/testthat/test_mcmc.R` (Issues 1–4), refactoring the logit
   dispatch path to use the two-step call.
5. Add a new `test_that` block for `compute_frequentist_sd`.
6. Run `devtools::test()` and confirm all tests pass.
