#!/usr/bin/env Rscript
#
# Test the manual derivatives using numerical differentiation.
# Effectively, this tests GetIVSEDerivs and GetRegressionSEDerivs with
# both grouped and ungrouped standard errors.

library(mrplew)
library(testthat)
library(tidyverse)

context("mrplew")


###############
# Test
test_that("block_bootstrap_works", {
    set.seed(1423098)
    DrawSamples <- function(n_draws) {
        draws1_mat <- rnorm(2 * n_draws) %>% matrix(ncol=2)
        draws2_mat <- draws1_mat + rnorm(2 * n_draws) %>% matrix(ncol=2)
        return(list(d1=draws1_mat, d2=draws2_mat))  
    }

    n_draws <- 50300

    n_sims <- 1000
    cov_draws <- array(NA, c(n_sims, 2, 2))
    for (sim in 1:n_sims) {
    draws <- DrawSamples(n_draws)
    cov_draws[sim,,] <- cov(draws$d1, draws$d2)
    }
    true_se <- apply(cov_draws, sd, MARGIN=c(2,3))

    draws <- DrawSamples(n_draws)
    cov(draws$d1, draws$d2) # Should be roughly the identity

    num_blocks <- 500
    num_draws <- 600
    draws1_mat <- draws$d1
    draws2_mat <- draws$d2
    cov_fun <- cov
    show_progress_bar <- TRUE
    cov_se_list <- GetBlockBootstrapCovarianceDraws(
        draws1_mat, draws2_mat, num_blocks, num_draws, show_progress_bar=TRUE
    )

    prop_err <- (cov_se_list$cov_se - true_se) / true_se

    AssertNearlyEqual(prop_err, 0, tol=0.08)
})

