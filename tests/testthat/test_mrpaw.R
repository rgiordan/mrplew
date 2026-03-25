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



test_that("simulations_correct", {

  ##########################################################
  # Simualte data

  set.seed(42)

  n_groups <- 2
  degree <- 2
  n_obs <- 100000
  n_obs_pop <- 100000

  sim_data <- SimulateSurveyData(n_groups, n_obs, n_obs_pop, degree=degree)
  group_effects <- sim_data$group_effects
  survey_df <- sim_data$survey_df
  pop_df <- sim_data$pop_df

  g_cols <- group_effects$g_cols

  # Check that the averages within the categories are roughly the predicted
  # averages
  err_z <- 
    survey_df %>%
    group_by(ey) %>%
    summarize(ybar=mean(y), n=n()) %>%
    summarize(max_err=max(abs(ey - ybar) / sqrt(ey * (1 - ey) / n)))
  
  expect_true(err_z < 3.0, info="Simulated means are not expected means")


  ##########################################################
  # Aggregate across groups to check 

  # Note that the check column is not in the AggregateSimulationData function
  survey_agg_df <-
    survey_df %>%
    group_by(pick(all_of(c(g_cols, "s")))) %>%
    summarize(check=sd(ey), count=n(), ybar=mean(y), ey=mean(ey),
              .groups="drop") %>%
    mutate(w=count / sum(count))

  pop_agg_df <-
    pop_df %>%
    ungroup() %>%
    group_by(pick(all_of(c(g_cols, "s")))) %>%
    summarize(check=sd(ey), count=n(), ybar=mean(y), ey=mean(ey),
              .groups="drop") %>%
    mutate(w=count / sum(count))

  # Sanity check that the ey values are the same within every category
  AssertNearlyZero(pop_agg_df$check, tol=1e-6)
  AssertNearlyZero(survey_agg_df$check, tol=1e-6)

  # Compare the averages and weights directly from the population
  # to the sample.

  joint_df <- left_join(
    pop_agg_df,
    survey_agg_df,
    by=c(g_cols, "s"),
    suffix=c("_pop", "_sur"))
  AssertNearlyEqual(joint_df$ey_sur, joint_df$ey_pop, tol=1e-6)

})





test_that("covariate_balance", {
  # TODO!
})


test_that("influence_weighting", {
  # TODO!
})