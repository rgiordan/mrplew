library(tidyverse)
library(lme4)
library(brms)
library(tidybayes)
library(mrplew)
library(gridExtra)

AEqBLine <- function() {
  geom_abline(aes(slope=1, intercept=0))
}

##########################################
# Simulate some categories
# 

set.seed(25338)

# induce group probabilities with a truncated normal
n_groups <- 2
degree <- 2
n_obs <- 1000
n_obs_pop <- 100000

# Simulate some data.
sim_data <- SimulateSurveyData(n_groups, n_obs, n_obs_pop, degree=degree)

# Survey data:
survey_df <- sim_data$survey_df

# Population data:
pop_df <- sim_data$pop_df

mrp_true <- with(pop_df, mean(ey))
print(mrp_true)

# If there is imbalance, this should differ from the true mrp.
with(survey_df, mean(ey))

# Accumulate data within group.  This helps speed up MCMC prediction.
agg_list <- AggregateSimulationData(sim_data)
survey_agg_df <- agg_list$survey_agg_df
pop_agg_df <- agg_list$pop_agg_df

# The optimal weights are the ratio of the population and sample
# weights within the group.  Get the optimal weight for each row
# of the survey for comparison with implicit weights.
joint_df <- agg_list$joint_df
w_opt <-
  survey_df %>%
  inner_join(select(joint_df, s, w_opt), by="s") %>%
  pull(w_opt)


#####################################
# Run logistic and OLS regression


g_sum <- paste(sim_data$group_effects$g_cols, collapse=" + ")
reg_form <-sprintf( "y ~ 1 + (%s)^%d", g_sum, degree)

logit_fit <- glm(formula(reg_form), survey_df, family=binomial(link="logit"))
lm_fit <- lm(formula(reg_form), survey_df)

coefficients(lm_fit)
coefficients(logit_fit)

mrp_ols_weights <- get_ols_weights(lm_fit, survey_df, pop_df)
mrp_logit_weights <- get_logit_weights(logit_fit, survey_df, pop_df)

cat(paste(
  mrp_ols_weights$mrp, 
  mrp_logit_weights$mrp, 
  mrp_true, collapse=", "), "\n")
if (FALSE) {
  grid.arrange(
    qplot(mrp_ols_weights$w, w_opt) + AEqBLine(),
    qplot(mrp_logit_weights$w, w_opt) + AEqBLine()
  )
}




##########################################
# Run posterior samplers

num_draws <- 5000

# Get the MrP posterior draws and weights for the logit model
stan_time <- Sys.time()
logit_post <- brm(formula(reg_form), survey_df, family=bernoulli(link="logit"),
                  chains=4, cores=4, seed=1543, warmup=500, iter=num_draws,
                  file="logit_posterior")
stan_time <- Sys.time() - stan_time
print(stan_time)

if (FALSE) {
  # Sanity check that the logistic regression and posterior match
  plot(fixef(logit_post)[, "Estimate"], coefficients(logit_fit)); abline(0,1)
}

# get_logit_mcmc_weights also computes draws of MrP so there is no
# ambiguity about how we are estimating it.
logit_mcmc_mrp <- get_logit_mcmc_weights(
  logit_post, survey_df, pop_agg_df, pop_w=pop_agg_df$w)

cat(mean(logit_mcmc_mrp$mrp_draws), ", ", mrp_true, "\n")
cat(mean(logit_mcmc_mrp$mrp_draws), ", ", mrp_ols_weights$mrp, "\n")
if (FALSE) {
  # Compare the logistic regression weights to the MrPlew MCMC weights
  qplot(logit_mcmc_mrp$w, w_opt) + AEqBLine()
}



# Get the MrP posterior draws and weights for the normal model
stan_time <- Sys.time()
lin_post <- brm(formula(reg_form), survey_df, family=gaussian(),
                chains=4, cores=4, seed=1543, warmup=500, iter=num_draws,
                file="ols_posterior")
stan_time <- Sys.time() - stan_time
print(stan_time)

# get_logit_mcmc_weights also computes draws of MrP so there is no
# ambiguity about how we are estimating it.
lin_mcmc_mrp <- get_ols_mcmc_weights(
  lin_post, survey_df, pop_agg_df, pop_w=pop_agg_df$w)

cat(mean(lin_mcmc_mrp$mrp_draws), ", ", mrp_true, "\n")
cat(mean(lin_mcmc_mrp$mrp_draws), ", ", mrp_ols_weights$mrp, "\n")
if (FALSE) {
  # Compare the OLS equivalent weights to the MrPlew MCMC weights
  qplot(lin_mcmc_mrp$w, w_opt) + AEqBLine()
}





