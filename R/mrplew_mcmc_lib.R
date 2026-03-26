library(tidyverse)
library(brms)


##############################################
# MCMC

#' Get mrplew weights for MCMC estimators.
#'
#' @param yhat_pop_draws MCMC draws x Population size matrix of yhat posterior draws
#' @param dloglikdy_survey_draws MCMC draws x Survey size matrix of draws of d log p(y_i | theta) / dy_i
#' @param pop_w Optional.  The weight given to each row of pop_df.  Defaults to ones.
#'
#' @return Draws from the MrP estimate, and the weight vector
#' whose i-th entry is d E[MrP | X, Y] / d y_i where i indexes survey observations.
#'
#'@export
get_mrplew_mcmc <- function(yhat_pop_draws, linpred_survey_draws, pop_w=NULL) {
    stopifnot(is.matrix(yhat_pop_draws))
    stopifnot(is.matrix(dloglikdy_survey_draws))
    stopifnot(nrow(yhat_pop_draws) == nrow(dloglikdy_survey_draws)) # Number of MCMC draws should be the same
    
    if (is.null(pop_w)) {
        pop_w <- rep(1, ncol(yhat_pop_draws))
    } else {
        stopifnot(length(pop_w) == ncol(yhat_pop_draws))
    }
    mrp_draws <- yhat_pop %*% pop_w
    mrplew_w <- cov(mrp_draws, dloglikdy_survey_draws)[1,]

    result_list <- list(
        mrp_draws=mrp_draws,
        mrplew_w=mrplew_w
    )
}


