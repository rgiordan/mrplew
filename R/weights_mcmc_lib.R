library(tidyverse)
library(brms)


##############################################
# MCMC


#' Get mrplew weights for the logistic MCMC estimator.
#' @param logit_post The output of `brm(..., survey_df, family=binomial(link="logit"))`
#' @param survey_df The survey dataframe
#' @param pop_df The population dataframe
#' @param pop_w Optional.  The weight given to each row of pop_df.  Defaults to ones.
#' @param save_preds Optional.  If true, save the posterior predictions for re-use.
#' @param re_formula Optional.  Formula containing group-level effects to be considered in the prediction. If `NULL` (default), include all group-level effects; if `NA`, include no group-level effects.
#' @param allow_new_levels Optional.  If true, allow new levels of group-level effects in prediction stage.
#'
#' @return Draws from the MrP estimate, and the weight vector
#' whose n-th entry is d E[MrP | X, Y] / d y_n.
#'
#' @importFrom brms posterior_epred
#' @importFrom brms posterior_linpred
#'@export
GetLogitMCMCWeights <- function(logit_post, survey_df, pop_df, pop_w=NULL, 
                                save_preds=FALSE, re_formula=NULL,
                                allow_new_levels=FALSE) {
    stopifnot(class(logit_post) == "brmsfit")

    CheckLogitFamily(logit_post)

    pop_w <- GetPopulationWeights(pop_df, pop_w)

    # posterior_epred should be yhat.
    # posterior_linpred should be theta^T x_n.  
    # Draws are in rows and observations in columns.
    yhat_pop <- posterior_epred(logit_post, newdata=pop_df,
                                re_formula=re_formula,
                                allow_new_levels=allow_new_levels)
    mrp_draws_logit <- yhat_pop %*% pop_w

    # Get the influence scores for the logit model
    # The log likelihood derivative for the n^th datapoint is just the theta^T x_n
    # TODO: optionally return the linpred for further diagnostics
    ll_grad_draws_logit <- posterior_linpred(logit_post, newdata=survey_df)
    w_logit <- cov(mrp_draws_logit, ll_grad_draws_logit)[1,]

    result_list <- list(
        mrp_draws=mrp_draws_logit,
        w=w_logit
    )

    if (save_preds) {
        yhat_draws <- posterior_epred(logit_post, newdata=survey_df)
        eta_draws <- posterior_linpred(logit_post, newdata=survey_df)

        result_list$yhat_pop <- yhat_pop
        result_list$yhat_draws <- yhat_draws
        result_list$eta_draws <- eta_draws
    }
    return(result_list)
}



GetOLSLikelihoodComponentDraws <- function(lin_post, survey_df) {
    # posterior_epred should be yhat.
    # posterior_linpred should be theta^T x_n.  
    # Draws are in rows and observations in columns.

    # get_variables(lin_post)
    sigma_draws <- lin_post %>% spread_draws(sigma) %>% pull(sigma)
    yhat_draws <- posterior_linpred(lin_post, newdata=survey_df)
    stopifnot(ncol(yhat_draws) == nrow(survey_df))
    y <- GetResponse(lin_post)
    resid_draws <- (y - t(yhat_draws)) %>% t()
    return(list(resid_draws=resid_draws, sigma_draws=sigma_draws, yhat_draws=yhat_draws))
}




#' Get mrplew weights for the logistic MCMC estimator.
#' @param lin_post The output of `brm(..., survey_df, family=gaussian())`
#' @param survey_df The survey dataframe
#' @param pop_df The population dataframe
#' @param pop_w Optional.  The weight given to each row of pop_df.  Defaults to ones.
#' @param re_formula Optional.  Formula containing group-level effects to be considered in the prediction. If `NULL` (default), include all group-level effects; if `NA`, include no group-level effects.
#' @param allow_new_levels Optional.  If true, allow new levels of group-level effects in prediction stage.
#'
#' @return Draws from the MrP estimate, and the weight vector
#' whose n-th entry is d E[MrP | X, Y] / d y_n.
#'
#' @importFrom brms posterior_epred
#' @importFrom brms posterior_linpred
#'@export
GetOLSMCMCWeights <- function(lin_post, survey_df, pop_df, pop_w=NULL, 
                              re_formula=NULL, allow_new_levels=FALSE) {
    stopifnot(class(lin_post) == "brmsfit")
    CheckOLSFamily(lin_post)

    # TODO: here and elsewhere check that "y" is actually the response in the
    # formula, or else extract it.

    # This used to be necessary but is now covered by GetResponse()
    # stopifnot("y" %in% names(survey_df))
    # print("Warning: the response variable must be y for this function to work!")

    pop_w <- GetPopulationWeights(pop_df, pop_w)
    yhat_pop <- posterior_epred(lin_post, newdata=pop_df, re_formula=re_formula,
                                allow_new_levels=allow_new_levels)

    mrp_draws_lin <- yhat_pop %*% pop_w
    ols_ll_draws <- GetOLSLikelihoodComponentDraws(lin_post, survey_df)

    # The log likelihood derivative for the n^th datapoint is
    # sigma^{-2} (y_n - \hat{y}_n)
    ll_grad_draws_lin <- -1 * ols_ll_draws$resid_draws / (ols_ll_draws$sigma_draws^2)
    w_lin <- cov(mrp_draws_lin, ll_grad_draws_lin)[1,]

    return(list(
        mrp_draws=mrp_draws_lin,
        w=w_lin
    ))
}