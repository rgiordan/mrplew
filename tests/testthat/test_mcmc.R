#!/usr/bin/env Rscript
#
# Test the manual derivatives using numerical differentiation.
# Effectively, this tests GetIVSEDerivs and GetRegressionSEDerivs with
# both grouped and ungrouped standard errors.

library(mrplew)
library(testthat)
library(tidyverse)
library(brms)
library(rlang)

context("mrplew")




################################################################
# Get the MrP posterior draws and weights for the normal model
test_that("mcmc_runs", {
  for (method in c("ols", "logit", "ols_response_name", "logit_response_name")) {
    print(sprintf("Testing MCMC for method %s", method))
    if (method %in% c("ols", "ols_response_name")) {
      model_type <- "ols"
    } else if (method %in% c("logit", "logit_response_name")) {
      model_type <- "logit"
    } else {
      expect_true(FALSE, sprintf("mcmc_runs: Unknown method %s", method))
    }

    rds_load <- SafeLoadPosterior(method)
    post <- rds_load$post
    sim_data <- rds_load$sim_data

    y_col <- as.character(f_lhs(as.formula(formula(post))))
    agg_list <- GetAggData(sim_data, y_col)

    # Test that this runs and produces weights of the correct length.
    if (model_type == "ols") {
      mcmc_mrp <- get_ols_mcmc_weights(
        post,
        sim_data$survey_df,
        agg_list$pop_agg_df,
        pop_frac=agg_list$pop_agg_df$frac)
    } else if (model_type == "logit") {
      mrp_draws <- get_mrp_draws_brms(
        brms_post=post,
        pop_df=agg_list$pop_agg_df,
        pop_frac=agg_list$pop_agg_df$frac)
      mcmc_mrp <- get_mrplew_logistic_brms(
        brms_post=post,
        survey_df=sim_data$survey_df,
        mrp_draws=mrp_draws)
    }

    expect_true(length(mcmc_mrp$mrplew_w) == nrow(sim_data$survey_df))

    print("Testing likelihood")
    # Test the likelihood computation
    yhat_pop <- brms::posterior_epred(post, newdata=agg_list$pop_agg_df)
    linpred_pop <- brms::posterior_linpred(post, newdata=agg_list$pop_agg_df)

    if (model_type == "ols") {
      AssertNearlyEqual(linpred_pop, yhat_pop)

      # Sanity check that I'm computing the log likelihood correctly
      # (I'm not taking into account the prior so there will be some small mismatch)
      ols_ll_draws <- get_ols_likelihood_component_draws(post, sim_data$survey_df)

      sigma_draws <- ols_ll_draws$sigma_draws
      resid_draws <- ols_ll_draws$resid_draws
      lp_draws_check <- post %>% spread_draws(lp__) %>% pull(lp__)
      lp_mat <- -0.5 * (resid_draws^2) / (sigma_draws^2) - log(sigma_draws)
      lp_draws <- apply(lp_mat, FUN=sum, MARGIN=1)
      expect_true(cor(lp_draws, lp_draws_check) > 0.99)

    } else if (model_type == "logit") {
      # posterior_epred should be yhat.
      # posterior_linpred should be theta^T x_n.
      # Draws are in rows and observations in columns.
      AssertNearlyEqual(expit(linpred_pop), yhat_pop)

      # Sanity check that I'm computing the log likelihood correctly
      # (I'm not taking into account the prior so there will be some small mismatch)
      eta_draws <- posterior_linpred(post, newdata=sim_data$survey_df)
      y <- get_response(post)
      lp_mat <- (y * t(eta_draws) - log(1 + exp(t(eta_draws)))) %>% t()
      lp_draws <- apply(lp_mat, FUN=sum, MARGIN=1)
      lp_draws_check <- post %>% spread_draws(lp__) %>% pull(lp__)
      expect_true(cor(lp_draws, lp_draws_check) > 0.99)


    } else {
      expect_true(FALSE, sprintf("mcmc_runs: Unknown model type %s", model_type))
    }
  }

})
