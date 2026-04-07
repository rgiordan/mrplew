library(tidyverse)
library(lme4)
library(brms)
library(tidybayes)
library(mrplew)
library(gridExtra)

##########################################
# Simulate some MrP data

set.seed(25338)

# induce group probabilities with a truncated normal
n_groups <- 4
degree <- 2
n_obs <- 2000
n_obs_pop <- 5000

# Simulate some data.
sim_data <- 
  mrplew::simulate_survey_data(n_groups, n_obs, n_obs_pop, degree=degree)

# Survey data:
survey_df <- sim_data$survey_df

# Population data:
pop_df <- sim_data$pop_df

mrp_true <- with(pop_df, mean(ey))
print(mrp_true)
mrp_df <- data.frame(truth=mrp_true)

# If there is imbalance, this should differ from the true mrp.
with(survey_df, mean(ey))

# Accumulate data within group.  This helps speed up MCMC prediction.
agg_list <- aggregate_simuilation_data(sim_data)
survey_agg_df <- agg_list$survey_agg_df
pop_agg_df <- agg_list$pop_agg_df

# The optimal weights are the ratio of the population and sample
# weights within the group.  Get the optimal weight for each row
# of the survey for comparison with implicit weights.
joint_df <- agg_list$joint_df
w_opt <-
  survey_df %>%
  left_join(select(joint_df, s, w_opt), by="s") %>%
  pull(w_opt)

weights_df <- data.frame(w_opt=w_opt)
if (FALSE) {
  # Look at the dispersion of the optimal weights
  hist(w_opt, 100)
}

#####################################
# Run OLS regression and get equivalent weights

g_sum <- paste(sim_data$group_effects$g_cols, collapse=" + ")
reg_form <-sprintf( "y ~ 1 + (%s)^%d", g_sum, degree)

lm_fit <- lm(formula(reg_form), survey_df)
ols_mrplew <- get_mrplew_lm(lm_fit, survey_df, pop_df)
mrp_df$ols <- ols_mrplew$mrp

weights_df$w_ols <- ols_mrplew$mrplew_w
if (FALSE) {
  ggplot(aes(x=w_opt, y=w_ols), data=weights_df) +
    geom_point() + geom_abline()
}



##########################################
# Run logistic MCMC posterior samplers

num_draws <- 5000

# Get the MrP posterior draws and weights for the logit model
stan_time <- Sys.time()
logit_post <- brm(formula(reg_form), survey_df, family=bernoulli(link="logit"),
                  chains=4, cores=4, seed=1543, warmup=500, iter=num_draws,
                  file="logit_posterior")
stan_time <- Sys.time() - stan_time
print(stan_time)

# get_logit_mcmc_weights also computes draws of MrP so there is no
# ambiguity about how we are estimating it.
logit_mcmc_mrplew <- get_mrplew_logistic_brms(
  logit_post, survey_df, pop_agg_df, pop_w=pop_agg_df$w)
mrp_df$mcmc <- mean(logit_mcmc_mrplew$mrp_draws)
weights_df$w_mcmc <- logit_mcmc_mrplew$mrplew_w

print(mrp_df)

if (FALSE) {
  # Compare the logistic regression weights to the MrPlew MCMC weights
  ggplot(aes(x=w_opt), data=weights_df) +
    geom_point(aes(y=w_ols, color="ols")) + 
    geom_point(aes(y=w_mcmc, color="mcmc logistic")) + 
    geom_abline()
}



##########################################
# Check balance

balance_reg_form <-
  formula(sprintf( "y ~ 1 + (%s)^%d", g_sum, degree + 2))

ols_balance_list <- check_covariate_balance(
  mrplew_list=ols_mrplew, 
  survey_df=survey_df, 
  pop_df=pop_agg_df, 
  reg_form=balance_reg_form, 
  pop_w=pop_agg_df$w)
ols_balance_df <-
  ols_balance_list$balance_df %>%
  mutate(pct=100 * difference / mrp_true)

mcmc_balance_df <- check_covariate_balance(
  mrplew_list=logit_mcmc_mrplew, 
  survey_df=survey_df, 
  pop_df=pop_agg_df, 
  reg_form=balance_reg_form, 
  pop_w=pop_agg_df$w)$balance_df %>%
  mutate(pct=100 * difference / mrp_true)


opt_balance_df <- get_balance_df(
  x1=ols_balance_list$x_pop,
  x2=ols_balance_list$x_survey, 
  w1=pop_agg_df$w, 
  w2=w_opt) %>%
  rename(pop=x1bar, survey=x2bar) %>%
  mutate(pct=100 * difference / mrp_true)


balance_df <- rbind(
  ols_balance_df %>% mutate(name="ols"),
  mcmc_balance_df %>% mutate(name="logistic mcmc"),
  opt_balance_df %>% mutate(name="optimal"))
jitter_amount <- max(abs(balance_df$pct)) * 5e-4
ggplot(balance_df) +
  geom_bar(
    aes(fill=name, y=pct + jitter_amount, x=reg),
    position="dodge", stat="identity") +
  coord_flip() +
  xlab(NULL) +
  labs(fill="Method") +
  ylab("Imbalance (% of MrP truth)")


compute_mrplew_variance_estimate <- function(w, y, yhat) {
  stopifnot(length(y) == length(w))
  stopifnot(length(yhat) == length(w))
  resid <- y - yhat
  weps <- resid * w
  weps_bar <- mean(weps)
  return(mean(weps^2) - weps_bar^2)
}

compute_mrplew_variance_estimate(ols_mrplew$mrplew_w, y=survey_df$y, yhat)


