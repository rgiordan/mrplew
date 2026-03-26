
#' Encode two vectors to factors ensuring that the levels are the same.
#' @param vec1 A vector to be converted to a factor
#' @param vec2 A vector to be converted to a factor
#'
#' @return A list containing consistent factor versions of vec1 and vec2
convert_to_consistent_factors <- function(vec1, vec2) {
  all_unique_levels <- sort(unique(c(as.character(vec1), as.character(vec2))))
  factor1 <- factor(vec1, levels = all_unique_levels)
  factor2 <- factor(vec2, levels = all_unique_levels)
  
  return(list(factor1 = factor1, factor2 = factor2))
}


#' A version of model.matrix that respects na.action.
#'
#' @param form A formula
#' @param df A dataframe
#'
#' @return The output of model.matrix(form, df) but respecting na.action.
safe_model_matrix <- function(form, df, na.action=na.pass) {
  # It appears that model.matrix does not actually respect na.fail or na.pass.
  # foo <- data.frame(x=c(1, 2, NA))
  # model.matrix(~x , foo, na.action="na.fail")
  # model.matrix(~x , foo, na.action=na.fail)
  # model.matrix(~x , foo, na.action="na.pass")
  # model.matrix(~x , foo, na.action=na.pass)
  
  terms_obj <- terms(form, data=df, na.action=na.action)
  mf <- model.frame(terms_obj, data=df, na.action=na.action)
  return(model.matrix(terms_obj, data=mf))
}


#' Create regressor matrices with consistent encodings.
#'
#' @param form A model formula
#' @param df1 A dataframe for model matrix
#' @param df2 A dataframe for model matrix
#' @param frame Use model.frame rather than model.matrix.
#'
#' @return A list containing regressor dataframes x1 and x2 for df1 and df2
#' respectively with consistent encodings.  NAs will be included.
#'
#' @export
get_consistent_regressors <- function(form, df1, df2, frame=FALSE) {
  # Generate design matrices from a formula in a way
  # that ensures the same encoding conventions are used for
  # each of two dataframes
  
  # Make a randomly named column to keep track of which df is which
  id_col <- paste0("df_source_", sample(1e12, 1))
  
  if (id_col %in% c(names(df1), names(df2))) {
    sprintf("ID column %d is already present, which should never happen.",
            id_cols)
  }
  df1[[id_col]] <- 1
  df2[[id_col]] <- 2
  df_comb <- bind_rows(df1, df2)
  
  if (frame) {
    x <- model.frame(form, df_comb, na.action=na.pass) %>% as.data.frame()
  } else {
    x <- safe_model_matrix(form, df_comb, na.action=na.pass) %>%
      as.data.frame()
  }
  x[[id_col]] <- df_comb[[id_col]]
  x1 <- x[ x[[id_col]] == 1, ]
  x2 <- x[ x[[id_col]] == 2, ]
  x1[[id_col]] <- NULL
  x2[[id_col]] <- NULL
  return(list(x1=x1, x2=x2))
}


check_balance_matrices <- function(x1, x2) {
  stopifnot(ncol(x1) == ncol(x2))
  if (any(colnames(x2) != colnames(x2))) {
    warning_text <- paste0(
      "The column names of the covariates do not match: ",
      paste0("(", colnames(x1), ")", collapse=", "),
      " versus ",
      paste0("(", colnames(x2), ")", collapse=", ")
    )
    warn(warning_text)
  }
}

#'@export
get_balance_df <- function(x1, x2, w1, w2) {
  check_balance_matrices(x1, x2)
  stopifnot(nrow(x1) == length(w1))
  stopifnot(nrow(x2) == length(w2))
  x1bar <- t(x1) %*% w1
  x2bar <- t(x2) %*% w2
  balance_df <- data.frame(reg=colnames(x1), x1bar=x1bar, x2bar=x2bar, difference=x1bar - x2bar)
  rownames(balance_df) <- NULL
  return(balance_df)
}


#' @export
check_covariate_balance <- function(mrpaw_list, survey_df, pop_df, reg_form, pop_w=NULL) {
  x_balance <- get_consistent_regressors(form=reg_form, df1=pop_df, df2=survey_df)
  x_pop <- x_balance$x1
  x_survey <- x_balance$x2
  stopifnot(sum(is.na(x_poststrat)) == 0)
  stopifnot(sum(is.na(x_survey)) == 0)
  pop_w <- get_population_weights(pop_df, pop_w)
  balance_df <- 
    get_balance_df(x1=x_pop, x2=x_survey, w1=pop_w, w2=mrpaw_list$w) %>%
    rename(poststrat=x1bar, survey=x2bar)
  return(list(balance_df=balance_df, x_pop=x_pop, x_survey=x_survey))
}
