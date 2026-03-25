#!/usr/bin/env Rscript
#
# Test the manual derivatives using numerical differentiation.
# Effectively, this tests GetIVSEDerivs and GetRegressionSEDerivs with
# both grouped and ungrouped standard errors.

library(mrplew)
library(testthat)
library(tidyverse)
library(brms)

context("mrplew")



test_that("ols_works", {
  #####################################
  # Run correctly specified logistic and OLS regression

  sim_data <- GetTestData()
  agg_list <- AggregateSimulationData(sim_data)
  group_effects <- sim_data$group_effects

  g_cols <- sim_data$group_effects$g_cols
  g_sum <- paste(group_effects$g_cols, collapse=" + ")
  reg_form <-sprintf( "y ~ 1 + (%s)^%d", g_sum, degree=2)

  for (method in c("ols", "logit")) {
    if (method == "ols") {
      fit <- lm(formula(reg_form), sim_data$survey_df)
      mrp_weights <- GetOLSWeights(fit, sim_data$survey_df, sim_data$pop_df)
    } else if (method == "logit") {
      fit <- glm(formula(reg_form), sim_data$survey_df, family=binomial(link="logit"))
      mrp_weights <- GetLogitWeights(fit, sim_data$survey_df, sim_data$pop_df)
    } else {
      expect_true(FALSE, "This should never happen")
    }
  }

  # Because the model is correctly specified, the fit should be exact
  yhat <- predict(fit, agg_list$survey_agg_df, type="response")
  AssertNearlyEqual(yhat, agg_list$survey_agg_df$ybar, tol=1e-9)
  
  # Because the model is correctly and dense, the weights should be
  # the optimal weights
  w_opt <-
    sim_data$survey_df %>%
    inner_join(select(agg_list$joint_df, s, w_opt), by="s") %>%
    pull(w_opt)

  AssertNearlyEqual(mrp_weights$w, w_opt)
})


