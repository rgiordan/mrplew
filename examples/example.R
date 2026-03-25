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

mrp_ols_weights <- GetOLSWeights(lm_fit, survey_df, pop_df)
mrp_logit_weights <- GetLogitWeights(logit_fit, survey_df, pop_df)

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

logit_post_file <- "logit_post.rds"

if (file.exists(logit_post_file)) {
  rds_load <- readRDS(logit_post_file)
  logit_post <- rds_load$logit_post
  sim_data <- rds_load$sim_data

  agg_list <- AggregateSimulationData(sim_data)
  survey_agg_df <- agg_list$survey_agg_df
  pop_agg_df <- agg_list$pop_agg_df
  joint_df <- agg_list$joint_df

} else {
  stan_time <- Sys.time()
  logit_post <- brm(formula(reg_form), survey_df, family=bernoulli(link="logit"),
                    chains=4, cores=4, seed=1543, warmup=500, iter=num_draws)
  stan_time <- Sys.time() - stan_time
  print(stan_time)
  saveRDS(list(
    logit_post=logit_post, 
    sim_data=sim_data), 
    file=logit_post_file)  
}

if (FALSE) {
  # Sanity check
  plot(fixef(logit_post)[, "Estimate"], coefficients(logit_fit)); abline(0,1)
}

# GetLogitMCMCWeights also computes draws of MrP so there is no
# ambiguity about how we are estimating it.
logit_mcmc_mrp <- GetLogitMCMCWeights(
  logit_post, survey_df, pop_agg_df, pop_w=pop_agg_df$w)

cat(mean(logit_mcmc_mrp$mrp_draws), ", ", mrp_true, "\n")
cat(mean(logit_mcmc_mrp$mrp_draws), ", ", mrp_ols_weights$mrp, "\n")
if (FALSE) {
  qplot(logit_mcmc_mrp$w, w_opt) + AEqBLine()
}



################################################################
# Get the MrP posterior draws and weights for the normal model


ols_post_file <- "ols_post.rds"

if (file.exists(ols_post_file)) {
  rds_load <- readRDS(ols_post_file)
  lin_post <- rds_load$lin_post
  sim_data <- rds_load$sim_data

  agg_list <- AggregateSimulationData(sim_data)
  survey_agg_df <- agg_list$survey_agg_df
  pop_agg_df <- agg_list$pop_agg_df
  joint_df <- agg_list$joint_df

} else {
  stan_time <- Sys.time()
  lin_post <- brm(formula(reg_form), survey_df, family=gaussian(),
                  chains=4, cores=4, seed=1543, warmup=500, iter=num_draws)
  stan_time <- Sys.time() - stan_time
  print(stan_time)
  saveRDS(list(
    lin_post=lin_post, 
    sim_data=sim_data), 
    file=ols_post_file)  
}
GetOLSMCMCWeights
# GetLogitMCMCWeights also computes draws of MrP so there is no
# ambiguity about how we are estimating it.
lin_mcmc_mrp <- GetOLSMCMCWeights(
  lin_post, survey_df, pop_agg_df, pop_w=pop_agg_df$w)

cat(mean(lin_mcmc_mrp$mrp_draws), ", ", mrp_true, "\n")
cat(mean(lin_mcmc_mrp$mrp_draws), ", ", mrp_ols_weights$mrp, "\n")
if (FALSE) {
  qplot(lin_mcmc_mrp$w, w_opt) + AEqBLine()
}





