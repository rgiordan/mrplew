library(tidyverse)


#'@export
logit <- function(x) {
  return(log(x / (1 - x)))
}


#'@export
expit <- function(x) {
  return(exp(x) / (1 + exp(x)))
}


# Use pop_w for weights if specfied, otherwise use 
# a vector of ones as long as pop_df.
get_population_weights <- function(pop_df, pop_w=NULL) {
  if (is.null(pop_w)) {
    pop_w <- rep(1, nrow(pop_df)) / nrow(pop_df)
  }
  
  weight_sum <- sum(pop_w)
  if (abs(weight_sum - 1) > 1e-6) {
    warning(sprintf("The population weights do not sum to one: %f", weight_sum))
  }
  return(pop_w)
}



check_logit_family <- function(logit_fit) {
    logit_family <- family(logit_fit)
    if (!(logit_family$family %in% c("binomial", "bernoulli"))) {
      warning(sprintf("Family is not binomial or bernoulli (%s)", logit_family$family))
    }
    if (logit_family$link != "logit") {
      warning(sprintf("Link is not logit (%s)", logit_family$link))
    }
}



check_ols_family <- function(lin_post) {
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
get_response <- function(post) {
    return(as.numeric(standata(post)$Y))
}



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
