library(tidyverse)




CheckColumnNames <- function(x_ols, x_pop) {
    # Right now we assume that the model matrices are the same.
    # If they aren't we throw an error.
    if (!all(colnames(x_ols) == colnames(x_pop))) {
        print("Warning: the column names do not match.")
        print("Columns in the population but not sample: ")
        print(setdiff(colnames(x_pop), colnames(x_ols)))
        print("Columns in the sample but not population: ")
        print(setdiff(colnames(x_ols), colnames(x_pop)))
        print("It is also possible that the columns are ordered differently.")
        stop()
    }
}


#' Get mrplew weights for the OLS estimator.
#' @param lm_fit The output of `lm`
#' @param survey_df The survey dataframe
#' @param pop_df The population dataframe
#'
#' @return The MrP OLS estimate for sanity checking, and the weight vector
#' whose n-th entry is d MrP / d y_n.
#'
#'@export
GetOLSWeights <- function(lm_fit, survey_df, pop_df, pop_w=NULL) {
    stopifnot(class(lm_fit) == "lm")
    pop_w <- GetPopulationWeights(pop_df, pop_w)

    # Strip the response from the formula since it might be missing
    # in the population dataframe.
    reg_form <- update(formula(lm_fit), NULL ~ .)
    x_ols <- model.matrix(reg_form, survey_df)
    x_pop <- model.matrix(reg_form, pop_df)

    CheckColumnNames(x_ols, x_pop)

    # The OLS model weights have a closed form.  (The logistic ones do
    # too, but I'll just put this in for now)
    # d \hat\beta / dy_n = (X^T X)^{-1} x_n
    # MrP = \sum_s w_s x_s^T \hat\beta

    betahat <- coefficients(lm_fit)
    yhat_pop <- x_pop %*% betahat
    xtx <- t(x_ols) %*% x_ols
    mrp_ols <- t(pop_w) %*% yhat_pop %>% as.numeric()
    w_ols <- t(pop_w) %*% x_pop %*% solve(xtx, t(x_ols)) %>% as.numeric()
    return(list(mrp=mrp_ols, w=w_ols))
}


# TODO: combine this and the above function with a class checker?

#' Get mrplew weights for the logistic estimator.
#' @param logit_fit The output of `glm(..., survey_df, family=binomial(link="logit"))`
#' @param survey_df The survey dataframe
#' @param pop_df The population dataframe
#'
#' @return The MrP Logistic estimate for sanity checking, and the weight vector
#' whose n-th entry is d MrP / d y_n.
#'
#'@export
GetLogitWeights <- function(logit_fit, survey_df, pop_df, pop_w=NULL) {
    stopifnot(class(logit_fit) == c("glm", "lm"))
    CheckLogitFamily(logit_fit)
    pop_w <- GetPopulationWeights(pop_df, pop_w)

    # Strip the response from the formula since it might be missing
    # in the population dataframe.
    reg_form <- update(formula(logit_fit), NULL ~ .)

    phat <- predict(logit_fit, survey_df, type="response")
    phat_pop <- predict(logit_fit, pop_df, type="response")
    vhat <- phat * (1 - phat)
    vhat_pop <- phat_pop * (1 - phat_pop)

    x_ols <- model.matrix(reg_form, survey_df)
    x_pop <- model.matrix(reg_form, pop_df)

    CheckColumnNames(x_ols, x_pop)

    # Mrp = pi^T expit(x_pop beta)
    # d Mrp / d betahat = pi^T (v_pop . x_pop)
    # hess = v_survey . x_ols
    # dbetahat_dy <- solve(hess, t(x_ols)) =>
    # d Mrp / d y = w^T (v_pop . x_pop) hess^{-1} x_survey
    hess <- t(x_ols * vhat) %*% x_ols
    mrp_chain_rule_term <- colSums(pop_w * vhat_pop * x_pop)
    w_logit <- t(mrp_chain_rule_term) %*% solve(hess, t(x_ols)) %>% as.numeric()

    mrp_logit <- sum(pop_w * phat_pop)

    return(list(mrp=mrp_logit, w=w_logit))
}










