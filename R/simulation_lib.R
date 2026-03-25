library(tidyverse)
library(mvtnorm)

#'@export
Logit <- function(x) {
  return(log(x / (1 - x)))
}


#'@export
Expit <- function(x) {
  return(exp(x) / (1 + exp(x)))
}


#' Draw correlated group memberships with a censored normal
#' @param n_obs The number of rows to draw
#' @param latent_mean A length `p` mean of the latent variable
#' @param latent_cov A `p` by `p` covariance of the latent variable
#'
#' @return An `nobs` by `p` matrix of group membership.  A row
#' is a member of the group if the corresponding component of its
#' latent normal is greater than zero.
#'
#'@export
DrawCorrelatedGroups <- function(n_obs, latent_mean, latent_cov) {
  p <- length(latent_mean)
  stopifnot(dim(latent_cov) == c(p, p))
  x <- rmvnorm(n_obs, mean=latent_mean, sigma=latent_cov)
  g <- 1 * (x >= 0)
  colnames(g) <- paste0("g", 1:p)
  return(g)
}



#' Draw group membership effect sizes
#' @param n_groups The number of distinct groups
#' @param degree The degree of interaction
#'
#' @return A list containing a number of dataframes.
#' - beta_df: The coefficients used to generate the effects
#' - regressor_df: The mapping between regressors (group interactions) and categories
#' - group_df: The mapping between individual group membership and categories  
#' - reg_form: The formula string used to generate the model matrix.
#' - g_cols: The names of the groups.
#' The dataframe shared columns are as follows:
#' - s: An integer identifying each unique combination of groups
#' - name: A readable name of a regressor, i.e. a group interaction term
#'@export
DrawGroupEffects <- function(n_groups, degree) {
  # Expand the columns into all possible combinations (i.e., categories s)
  g_cols <- paste0("g", 1:n_groups)
  group_df <-
    lapply(1:n_groups, \(g) c(0, 1)) %>%
    do.call(expand.grid, .) %>%
    mutate(s=1:n())
  names(group_df)  <- c(g_cols, "s")


  # Get a long dataframe with each regressor's value within
  # each category.  A single category contains multiple regressors.#
  g_sum <- paste0(g_cols, collapse=" + ")
  if (degree > 1) {
    reg_form <- sprintf("~ (%s)^%d - 1", g_sum, degree)
  } else {
    # Raising to degree one is not allowed in model.matrix
    reg_form <- sprintf("~ %s - 1", g_sum)
  }
  regressor_df <- 
    model.matrix(formula(reg_form), group_df) %>%
    as.data.frame() %>%
    inner_join(group_df, by=g_cols) %>%
    pivot_longer(cols=-s) %>%
    mutate(order=str_count(name, pattern=":") + 1)

  # Construct a coefficient for each regressor.
  beta_df <-
    regressor_df %>%
    group_by(name, order) %>%
    summarize(count=n(), .groups="drop") %>%
    mutate(beta=rnorm(n(), 0, sd=1.5)) %>%
    mutate(beta=beta / order) %>%
    mutate(beta=beta - mean(beta))
  
  return(list(
    beta_df=beta_df, 
    regressor_df=regressor_df, 
    group_df=group_df, 
    reg_form=reg_form,
    g_cols=g_cols))
}


#' Accumulate responses of group membership interactions
#' @param group_effects A list of group effects, e.g. as returned by `DrawGroupEffects`
#'
#' @return A dataframe with regressors, groups, and ey values suitable for
#' joining with raw data.
AccumulateInteractionEffects <- function(group_effects) {
  # Compute the effect size for each category.
  # Join the coefficients with the category memberships,
  # and add up all the effects that contribute to a category.
  
  effect_df <- 
    select(group_effects$beta_df, name, beta) %>%
    inner_join(group_effects$regressor_df, by="name") %>%
    mutate(gbeta=value * beta) %>%
    group_by(s) %>%
    summarize(logit_ey=sum(gbeta), .groups="drop") %>%
    mutate(ey=Expit(logit_ey)) %>%
    inner_join(group_effects$group_df, by="s")
  
  attr(effect_df, "reg_form") <- group_effects$reg_form
  return(effect_df)
}



#' Draw binary responses for a particular set of group memberships
#' @param g_matrix A matrix of group memberships, e.g. as returned by `DrawCorrelatedGroups`
#' @param effect_df A dataframe containing the columns in `g_matrix` as well as
#' expected response `ey`, e.g. as returned by `AccumulateInteractionEffects`.
#'
#' @return A dataframe containing group memberships `s` and responses `y` corresponding
#' to the rows of `g_matrix`.
#'@export
DrawResponse <- function(g_matrix, effect_df) {
  g_cols <- names(g_matrix)
  stopifnot(all(g_cols %in% names(effect_df)))
  reg_form <- attr(effect_df, "reg_form")
  response_df <- 
    model.matrix(formula(reg_form), data.frame(g_matrix)) %>%
    as.data.frame() %>%
    mutate(n=1:n()) %>%
    left_join(effect_df, by=g_cols) %>%
    mutate(y=1 * (runif(n()) <= ey)) %>%
    arrange(n)
  return(response_df)
}




#' Draw survey data with some reasonable defaults
#' @param n_groups The number of distinct groups
#' @param n_obs The number of survey observations
#' @param n_obs The number of population observations
#' @param degree The maximum degree of true interactions
#'
#' @return A list of simualted data
#'@export
SimulateSurveyData <- function(n_groups, n_obs, n_obs_pop, degree) {
  survey_mean <- 1.5 * 1:n_groups / (n_groups)
  survey_mean <- survey_mean - mean(survey_mean)
  survey_cov <- diag(n_groups)

  pop_mean <- rev(survey_mean)
  pop_cov <- diag(n_groups)

  g_survey <- DrawCorrelatedGroups(n_obs, survey_mean, survey_cov)
  g_pop <- DrawCorrelatedGroups(n_obs_pop, pop_mean, pop_cov)

  group_effects <- DrawGroupEffects(n_groups, degree)
  effect_df <- AccumulateInteractionEffects(group_effects)

  survey_df <- DrawResponse(g_survey, effect_df)
  pop_df <- DrawResponse(g_pop, effect_df)

  return(list(
    group_effects=group_effects,
    effect_df=effect_df,
    survey_df=survey_df,
    pop_df=pop_df
  ))
}





#' Aggregate simulation data by group.
#' @param sim_data Output of `SimulateSurveyData`
#'
#' @return A list containing aggregated and joined simulation data
#'@export
AggregateSimulationData <- function(sim_data) {

  ##########################################################
  # Aggregate across groups to check 

  group_effects <- sim_data$group_effects
  g_cols <- group_effects$g_cols

  survey_agg_df <-
    sim_data$survey_df %>%
    group_by(pick(all_of(c(g_cols, "s")))) %>%
    summarize(count=n(), ybar=mean(y), ey=mean(ey), .groups="drop") %>%
    mutate(w=count / sum(count))

  pop_agg_df <-
    sim_data$pop_df %>%
    group_by(pick(all_of(c(g_cols, "s")))) %>%
    summarize(count=n(), ybar=mean(y), ey=mean(ey), .groups="drop") %>%
    mutate(w=count / sum(count))

  # Compare the averages and weights directly from the population
  # to the sample.
  joint_df <- left_join(
    pop_agg_df,
    survey_agg_df,
    by=c(g_cols, "s"),
    suffix=c("_pop", "_sur")) %>%
    mutate(w_opt=w_pop / count_sur)

  return(list(joint_df=joint_df, survey_agg_df=survey_agg_df, pop_agg_df=pop_agg_df))
}


