

#' Get mrplew weights for the brms logistic MCMC estimator.
#'
#' @param brms_post The output of `brm(..., survey_df, family=binomial(link="logit"))`
#' @param pop_df The population dataframe
#' @param pop_frac  The weight given to each row of pop_df.  Defaults to 1/N
#' @param ...  Additional arguments passed to posterior_epred
#'
#' @return Draws from the MrP estimate
#'
#' @importFrom brms posterior_epred
#'@export
get_mrp_draws_brms <- function(brms_post, pop_df=NULL, pop_frac=NULL, ...) {
    stopifnot(class(brms_post) == "brmsfit")
    yhat_pop_draws <- posterior_epred(brms_post, newdata=pop_df, ...)
    mrp_draws <- get_mrp_draws(yhat_pop_draws=yhat_pop_draws, pop_frac=pop_frac)
    return(mrp_draws)
}

#' Get mrplew weights for the brms logistic MCMC estimator.
#'
#' @param brms_post The output of `brm(..., survey_df, family=binomial(link="logit"))`
#' @param survey_df The survey dataframe
#' @param mrp_draws Optional. Draws of the mrp estimate from the same posterior, brms_post.
#' @param pop_df Optional.  The population dataframe
#' @param pop_frac Optional.  The weight given to each row of pop_df.  Defaults to 1/N.
#' @param save_draws Optional.  If true, save the posterior predictions for re-use.
#'
#' @return Draws from the MrP estimate, and the weight vector
#' whose n-th entry is d E[MrP | X, Y] / d y_n.
#'
#' @importFrom brms posterior_epred
#' @importFrom brms posterior_linpred
#'@export
get_mrplew_logistic_brms <- function(brms_post, survey_df, 
                                     mrp_draws=NULL, 
                                     pop_df=NULL, pop_frac=NULL, 
                                     save_draws=FALSE) {
    stopifnot(class(brms_post) == "brmsfit")
    check_logit_family(brms_post)

    if (is.null(mrp_draws) & is.null(pop_df)) {
        stop("You must specify either mrp_draws or pop_df.")
    }
    if (is.null(mrp_draws)) {
        mrp_draws <- get_mrp_draws_brms(brms_post=brms_post, 
                                        pop_df=pop_df,
                                        pop_frac=pop_frac)
    }

    # Draws are in rows and observations in columns.
    # d log p(y | theta) / d y_i = theta^T x_i
    dloglikdy_survey_draws <- posterior_linpred(brms_post, newdata=survey_df)

    result_list <- get_mrplew_mcmc(
        mrp_draws=mrp_draws,
        dloglikdy_survey_draws=dloglikdy_survey_draws)

    if (save_draws) {
        yhat_survey_draws <- expit(dloglikdy_survey_draws)
        result_list$yhat_survey_draws <- yhat_survey_draws
        result_list$dloglikdy_survey_draws <- dloglikdy_survey_draws
    }
    return(result_list)
}





######################
# OLS

# Get the likelihood components of  posterior draws from a linear brms model
get_ols_likelihood_component_draws <- function(brms_post, survey_df) {
    # posterior_epred should be yhat.
    # posterior_linpred should be theta^T x_n.  
    # Draws are in rows and observations in columns.

    stopifnot(class(brms_post) == "brmsfit")
    check_ols_family(brms_post)

    # get_variables(lin_post)
    sigma_draws <- brms_post %>% spread_draws(sigma) %>% pull(sigma)
    yhat_draws <- posterior_linpred(brms_post, newdata=survey_df)
    stopifnot(ncol(yhat_draws) == nrow(survey_df))
    y <- get_response(brms_post)
    resid_draws <- (y - t(yhat_draws)) %>% t()
    return(list(resid_draws=resid_draws, sigma_draws=sigma_draws, yhat_draws=yhat_draws))
}

#'@export
get_gaussian_dloglikdy_draws <- function(resid_draws, sigma_draws) {
    # The log likelihood derivative for the n^th datapoint is
    # -sigma^{-2} (y_n - \hat{y}_n)
    dloglikdy_survey_draws <- -1 * resid_draws / (sigma_draws^2)
    return(dloglikdy_survey_draws)
}


#' Get mrplew weights for the brms logistic MCMC estimator.
#'
#' @param brms_post The output of `brm(..., survey_df, family=gaussian())`
#' @param survey_df The survey dataframe
#' @param pop_df The population dataframe
#' @param pop_frac Optional.  The weight given to each row of pop_df.  Defaults to 1/N.
#'
#' @return Draws from the MrP estimate, and the weight vector
#' whose n-th entry is d E[MrP | X, Y] / d y_n.
#'
#' @importFrom brms posterior_epred
#' @importFrom brms posterior_linpred
#'@export
get_ols_mcmc_weights <- function(brms_post, survey_df, pop_df, pop_frac=NULL, 
                              re_formula=NULL, allow_new_levels=FALSE) {
    stopifnot(class(brms_post) == "brmsfit")
    check_ols_family(brms_post)

    yhat_pop_draws <- posterior_epred(brms_post, newdata=pop_df)
    ols_ll_draws <- get_ols_likelihood_component_draws(brms_post, survey_df)
    dloglikdy_survey_draws <- get_gaussian_dloglikdy_draws(
        resid_draws=ols_ll_draws$resid_draws, sigma_draws=ols_ll_draws$sigma_draws)
    mrp_draws <- get_mrp_draws(yhat_pop_draws, pop_frac)

    result_list <- get_mrplew_mcmc(
        mrp_draws=mrp_draws,
        dloglikdy_survey_draws=dloglikdy_survey_draws)

    return(result_list)
}