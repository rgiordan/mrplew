library(tidyverse)


# Use pop_w for weights if specfied, otherwise use 
# a vector of ones as long as pop_df.
GetPopulationWeights <- function(pop_df, pop_w=NULL) {
    if (is.null(pop_w)) {
        pop_w <- rep(1, nrow(pop_df)) / nrow(pop_df)
    }

    weight_sum <- sum(pop_w)
    if (abs(weight_sum - 1) > 1e-6) {
        warning(sprintf("The population weights do not sum to one: %f", weight_sum))
    }
    return(pop_w)
}




CheckLogitFamily <- function(logit_fit) {
    logit_family <- family(logit_fit)
    if (!(logit_family$family %in% c("binomial", "bernoulli"))) {
      warning(sprintf("Family is not binomial or bernoulli (%s)", logit_family$family))
    }
    if (logit_family$link != "logit") {
      warning(sprintf("Link is not logit (%s)", logit_family$link))
    }
}





CheckOLSFamily <- function(lin_post) {
    post_family <- family(lin_post)
    if (post_family$family != "gaussian") {
      warning(sprintf("Family is not gaussian (%s)", post_family$family))
    }
    if (post_family$link != "identity") {
      warning(sprintf("Link is not identity (%s)", post_family$link))
    }
}



#' Get the response variable (y) from the posterior.
#' I don't see this use clearly documented, so I want to factor it out
#' for testing.
#' @param post A brms posterior
#'
#' @return The numeric response variable used for the posterior fitting
#'@export
GetResponse <- function(post) {
    return(as.numeric(standata(post)$Y))
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
#' @export
GetBlockBootstrapCovarianceDraws <- function(draws1_mat, draws2_mat,
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
  if (show_progress_bar) {
    pb <- txtProgressBar(min=1, max=num_draws, style=3)
  }

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
    return(outer_bar - d1_bar %*% t(d2_bar))
  }
  
  if (FALSE) {
    # Fast sanity check.  Both methods should give the same answer.
    all_block_inds <- do.call(c, block_inds)
    n_samples <- nrow(draws1_mat)
    cov(draws1_mat[all_block_inds, , drop=FALSE], 
        draws2_mat[all_block_inds, , drop=FALSE]) * (n_samples - 1) / n_samples - 
      ComputeCovariance(1:num_blocks)
  }
  
  for (draw in 1:num_draws) {
    if (show_progress_bar) {
      setTxtProgressBar(pb, draw)
    }
    block_ind_draws <- sample(1:num_blocks, num_blocks, replace=TRUE)
    cov_samples[draw, , ] <- ComputeCovariance(block_ind_draws)
  }
  if (show_progress_bar) {
    close(pb)
  }
  
  cov_se <- sqrt(n_factor) * apply(cov_samples, MARGIN=c(2, 3), sd)
  rownames(cov_se) <- colnames(draws1_mat)
  colnames(cov_se) <- colnames(draws2_mat)
  
  return(list(cov_samples=cov_samples, cov_se=cov_se))
}
