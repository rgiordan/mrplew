library(testthat)


GetTestData <- function() {
  set.seed(42)

  n_groups <- 2
  degree <- 2
  n_obs <- 1000
  n_obs_pop <- 10 * n_obs

  sim_data <- SimulateSurveyData(n_groups, n_obs, n_obs_pop, degree=degree)
  return(sim_data)
}


AggregateSimulationData <- function(sim_data, resp="y") {

  ##########################################################
  # Aggregate across groups to check 

  group_effects <- sim_data$group_effects
  g_cols <- group_effects$g_cols

  survey_agg_df <-
    sim_data$survey_df %>%
    group_by(pick(all_of(c(g_cols, "s")))) %>%
    summarize(count=n(), ybar=mean(.data[[resp]]), ey=mean(ey), .groups="drop") %>%
    mutate(w=count / sum(count))

  pop_agg_df <-
    sim_data$pop_df %>%
    group_by(pick(all_of(c(g_cols, "s")))) %>%
    summarize(count=n(), ybar=mean(.data[[resp]]), ey=mean(ey), .groups="drop") %>%
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


AssertNearlyEqual <- function(x, y, tol=1e-9, desc=NULL) {
  diff_norm <- max(abs(x - y))
  if (is.null(desc)) {
    info_str <- sprintf("%e > %e", diff_norm, tol)
  } else {
    info_str <- sprintf("%s: %e > %e", desc, diff_norm, tol)
  }
  expect_true(diff_norm < tol, info=info_str)
}


AssertNearlyZero <- function(x, tol=1e-15, desc=NULL) {
  x_norm <- max(abs(x))
  if (is.null(desc)) {
    info_str <- sprintf("%e > %e", x_norm, tol)
  } else {
    info_str <- sprintf("%s: %e > %e", desc, x_norm, tol)
  }
  expect_true(x_norm < tol, info=info_str)
}






################################################################
# Load saved posterior samples (or run if needed)

RunAndCachePosterior <- function(
  sim_data, reg_form, post_filename, family, num_draws) {
    stan_time <- Sys.time()
    post <- brm(formula(reg_form), sim_data$survey_df, family=family,
                chains=4, cores=4, seed=1543, 
                iter=num_draws)
    stan_time <- Sys.time() - stan_time
    print(stan_time)
    saveRDS(list(
      post=post, 
      sim_data=sim_data),
      file=post_filename)
}


# Load a posterior for a particular test method, and run MCMC
# if the file is not available.
SafeLoadPosterior <- function(method) {
  GetDefaultConfig <- function() {
    config <- list()
    sim_data <- GetTestData()
    config$sim_data <- sim_data
    group_effects <- sim_data$group_effects
    g_cols <- group_effects$g_cols
    g_sum <- paste(group_effects$g_cols, collapse=" + ")
    config$reg_form <-sprintf( "y ~ 1 + (%s)^%d", g_sum, degree=2)
    return(config)
  }

  if (method == "ols") {
    config <- GetDefaultConfig()
    config$family <- gaussian()
    config$post_filename <- test_path("mcmc_cache/ols_post_test.rds")
    config$num_draws <- 5000
  } else if (method == "logit") {
    config <- GetDefaultConfig()
    config$family <- bernoulli(link="logit")
    config$post_filename <- test_path("mcmc_cache/logit_post_test.rds")
    config$num_draws <- 5000
  } else if (method == "ols_response_name") {
    config <- GetDefaultConfig()
    config$family <- gaussian()
    config$post_filename <- test_path("mcmc_cache/ols_post_response_name_test.rds")
    config$num_draws <- 100
    config$sim_data$survey_df <- rename(config$sim_data$survey_df, resp=y)
    config$sim_data$pop_df <- rename(config$sim_data$pop_df, resp=y)
    config$reg_form <- sub("^y", "resp", config$reg_form)
  } else if (method == "logit_response_name") {
    config <- GetDefaultConfig()
    config$family <- bernoulli(link="logit")
    config$post_filename <- test_path("mcmc_cache/logit_post_response_name_test.rds")
    config$num_draws <- 100
    config$sim_data$survey_df <- rename(config$sim_data$survey_df, resp=y)
    config$sim_data$pop_df <- rename(config$sim_data$pop_df, resp=y)
    config$reg_form <- sub("^y", "resp", config$reg_form)
  } else {
    expect_true(FALSE, sprintf("SafeLoadPosterior: Unknown method %s", method))
  }

  print(config$post_filename)

  if (!file.exists(config$post_filename)) {
    print(sprintf("Could not locate file %s", config$post_filename))
    print("Re-running MCMC and saving result to rds file")
    RunAndCachePosterior(
      config$sim_data, config$reg_form, config$post_filename, 
      family=config$family, num_draws=config$num_draws)
  } else {
    print(sprintf("Found cached posterior samples in %s", config$post_filename))
  }

  print(sprintf("Reading posterior from file %s", config$post_filename))
  rds_load <- readRDS(config$post_filename)
  return(rds_load)
}

