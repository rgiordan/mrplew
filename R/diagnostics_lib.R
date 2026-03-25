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


#' Use mrplew to check invariance to weighting datapoints by covariate
#' (i.e. distribution shift aligned with a covariate).
#' 
#' @param infl_list The output of the EvalInfluenceFunction function
#' @param x_sur A matrix of some regressors in the survey data
#'
#' @return The derivative of MrP in the direction of a distribution
#' shift aligned with each regressor in x_sur.
#' 
#' @export
CheckDatapointWeighting <- function(infl_list, x_sur) {
    stopifnot(nrow(x_sur) == length(infl_list$infl_vec))
    return(colSums(x_sur * infl_list$infl_vec))
}
