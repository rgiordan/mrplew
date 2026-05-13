library(tidyverse)
library(brms)




##############################################
# MCMC

#' Get mrp draws MCMC estimators.
#'
#' @param yhat_pop_draws MCMC draws x Population size matrix of yhat posterior draws
#' @param pop_frac Optional.  The weight given to each row of pop_df.  Defaults to ones.
#'
#' @return Draws from the MrP estimate.
#'
#'@export
get_mrp_draws <- function(yhat_pop_draws, pop_frac=NULL) {
    if (is.null(pop_frac)) {
        npop <- ncol(yhat_pop_draws)
        pop_frac <- rep(1 / npop, npop)
    } else {
        stopifnot(length(pop_frac) == ncol(yhat_pop_draws))
    }
    mrp_draws <- yhat_pop_draws %*% pop_frac
    return(mrp_draws)
}

#' Get mrplew weights for MCMC estimators.
#'
#' @param mrp_draws MCMC draws of the MrP estimator
#' @param dloglikdy_survey_draws MCMC draws x Survey size matrix of draws of d log p(y_i | theta) / dy_i
#'
#' @return Draws from the MrP estimate, and the weight vector
#' whose i-th entry is d E[MrP | X, Y] / d y_i where i indexes survey observations.
#'
#'@export
get_mrplew_mcmc <- function(mrp_draws, dloglikdy_survey_draws) {
    if (!is.matrix(mrp_draws)) {
      mrp_draws <- matrix(mrp_draws, ncol=1)
    }
    stopifnot(is.matrix(dloglikdy_survey_draws))
    stopifnot(nrow(mrp_draws) == nrow(dloglikdy_survey_draws)) # Number of MCMC draws should be the same
    
    dmrp_dy <- cov(mrp_draws, dloglikdy_survey_draws)[1,]
    n_obs <- ncol(dloglikdy_survey_draws)
    mrplew_w <- n_obs * dmrp_dy

    result_list <- list(
        mrp=mean(mrp_draws),
        mrp_draws=mrp_draws,
        mrplew_w=mrplew_w,
        n_obs=n_obs,
        dmrp_dy=dmrp_dy
    )
    return(result_list)
}





#' Estimate Monte Carlo standard errors of sample covariances or
#' by block bootstrapping draws from an MCMC chain.
#'
#' @param draws1_mat One set of parameter draws.
#' @param draws2_mat Another set of parameter draws.
#' @param num_blocks The number of blocks in the block bootstrap.
#' @param num_draws The number of bootstrap draws.
#' @param show_progress_par.  Optional.  If TRUE, show a progress bar.
#' By default, FALSE.
#' @return A list containing the draws of the covariance cov_samples
#' and the estimated Monte Carlo sample errors in cov_se.
#'
#' @importFrom purrr reduce
#' @export
get_block_bootstrap_covariance_draws <- function(draws1_mat, draws2_mat,
                                                 num_blocks, num_draws,
                                                 show_progress_bar=FALSE) {
  
  if (nrow(draws1_mat) != nrow(draws2_mat)) {
    stop("draws1_mat and draws2_mat must have the same number of rows.")
  }
  
  num_samples <- nrow(draws1_mat)
  
  block_size <- floor(num_samples / num_blocks)
  
  # Correction factor if the number of blocked observations is not the same
  # as the original.
  n_factor <- (block_size * num_blocks) / num_samples
  
  # The indices of each block into the MCMC samples.
  block_inds <- lapply(
    1:num_blocks,
    function(ind) { (ind - 1) * block_size + 1:block_size })
  
  base_cov <- cov(draws1_mat, draws2_mat)
  cov_samples <- array(NA, c(num_draws, ncol(draws1_mat), ncol(draws2_mat)))
  mean1_samples <- array(NA, c(num_draws, ncol(draws1_mat)))
  mean2_samples <- array(NA, c(num_draws, ncol(draws2_mat)))

  if (show_progress_bar) {
    pb <- txtProgressBar(min=1, max=num_draws, style=3)
  }

  # Pre-compute the required first and second sample moments within each block
  ComputeSums <- function(draws_mat) {
    lapply(block_inds, \(inds) colSums(draws_mat[inds, , drop=FALSE ]))
  }
  sums1 <- ComputeSums(draws1_mat)
  sums2 <- ComputeSums(draws2_mat)

  outers12 <- lapply(
    block_inds, 
    \(inds) t(draws1_mat[inds, , drop=FALSE ]) %*% draws2_mat[inds, , drop=FALSE])

  ComputeCovariance <- function(block_ind_draws) {
    n_ind_draws <- length(block_ind_draws) * block_size
    AverageOverInds <- function(sim_list) {
      reduce(sim_list[block_ind_draws], \(x, y) x + y) / n_ind_draws
    }
    d1_bar <- AverageOverInds(sums1)
    d2_bar <- AverageOverInds(sums2)
    outer_bar <- AverageOverInds(outers12)
    return(list(
      cov=outer_bar - d1_bar %*% t(d2_bar),
      mean1=d1_bar,
      mean2=d2_bar))
  }
  
  # if (FALSE) {
  #   # Fast sanity check.  Both methods should give the same answer.
  #   all_block_inds <- do.call(c, block_inds)
  #   n_samples <- nrow(draws1_mat)
  #   cov(draws1_mat[all_block_inds, , drop=FALSE], 
  #       draws2_mat[all_block_inds, , drop=FALSE]) * (n_samples - 1) / n_samples - 
  #     ComputeCovariance(1:num_blocks)
  # }
  
  for (draw in 1:num_draws) {
    if (show_progress_bar) {
      setTxtProgressBar(pb, draw)
    }
    block_ind_draws <- sample(1:num_blocks, num_blocks, replace=TRUE)
    sample_list <- ComputeCovariance(block_ind_draws)
    cov_samples[draw, , ] <- sample_list$cov
    mean1_samples[draw, ] <- sample_list$mean1
    mean2_samples[draw, ] <- sample_list$mean2
  }
  if (show_progress_bar) {
    close(pb)
  }
  
  cov_se <- sqrt(n_factor) * apply(cov_samples, MARGIN=c(2, 3), sd)
  rownames(cov_se) <- colnames(draws1_mat)
  colnames(cov_se) <- colnames(draws2_mat)
  
  return(list(cov_samples=cov_samples, cov_se=cov_se, 
              mean1_samples=mean1_samples, mean2_samples=mean2_samples))
}
