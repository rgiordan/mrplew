###################################################
# Various ways to use the MrP weights

library(einsum)
library(tidyverse)
library(brms)
library(tidybayes)


#' @importFrom brms posterior_epred
#' @importFrom brms posterior_linpred
SafeGetYhatDraws <- function(mrplew_list, post, survey_df) {
    if ("yhat_draws" %in% names(mrplew_list)) {
        yhat_draws <- mrplew_list$yhat_draws
    } else {
        yhat_draws <- posterior_epred(post, newdata=survey_df)
    }
    return(yhat_draws)
}


#' @importFrom brms posterior_epred
#' @importFrom brms posterior_linpred
SafeGetEtaDraws <- function(mrplew_list, post, survey_df) {
    if ("eta_draws" %in% names(mrplew_list)) {
        eta_draws <- mrplew_list$eta_draws
    } else {
        eta_draws <- posterior_linpred(post, newdata=survey_df)
    }
    return(eta_draws)
}



#' Use MrPaw to check covariate balance.
#' 
#' @param mrpaw_list The output of one of the Get*MCMCWeights functions
#' @param x_sur A matrix of some regressors in the survey data
#' @param x_pop A matrix of the same regressors in the population data
#' @param pop_w (Optional) A vector of population weights. Taken to be all one if NULL.
#'
#' @return The weighted average regressors in each regressor matrix, where the
#' weights in the survey are given by the MrP affine weights.
#' The regressor matrices might be made with a call to model.matrix using the
#' same formula for each dataset.
#' 
#' @export
CheckCovariateBalance <- function(mrpaw_list, x_sur, x_pop, pop_w=NULL) {
  pop_w <- GetPopulationWeights(x_pop, pop_w)
  stopifnot(ncol(x_pop) == ncol(x_sur))
  stopifnot(nrow(x_sur) == length(mrpaw_list$w))
  if (any(colnames(x_pop) != colnames(x_sur))) {
    warning_text <- paste0(
      "The column names of the covariates do not match: ",
      paste0("(", colnames(x_pop), ")", collapse=", "),
      " versus ",
      paste0("(", colnames(x_sur), ")", collapse=", ")
    )
    warn(warning_text)
  }
  return(list(
    x_mean_pop=colSums(pop_w * x_pop),
    x_mean_sur=colSums(mrpaw_list$w * x_sur)
  ))
}

