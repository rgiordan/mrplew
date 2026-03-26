###################################################
# Various ways to use the MrP weights

library(einsum)
library(tidyverse)
library(brms)
library(tidybayes)


#' @importFrom brms posterior_epred
#' @importFrom brms posterior_linpred
safe_get_yhat_draws <- function(mrplew_list, post, survey_df) {
    if ("yhat_draws" %in% names(mrplew_list)) {
        yhat_draws <- mrplew_list$yhat_draws
    } else {
        yhat_draws <- posterior_epred(post, newdata=survey_df)
    }
    return(yhat_draws)
}


#' @importFrom brms posterior_epred
#' @importFrom brms posterior_linpred
safe_get_eta_draws <- function(mrplew_list, post, survey_df) {
    if ("eta_draws" %in% names(mrplew_list)) {
        eta_draws <- mrplew_list$eta_draws
    } else {
        eta_draws <- posterior_linpred(post, newdata=survey_df)
    }
    return(eta_draws)
}
