# Compute the maximum permissible delta for this perturbation and e_y.
# We want to draw binary data with mean e_y + delta * x
#
# We require
# 0 <= e_y + delta * x <= 1
# Solving gives
# - e_y / |x| <= delta * sign(x) <= (1 - e_y) / |x|.
#
# One of the two bounds is always trivial since delta >= 0,
# and which bound matters depends on the sign of x.
#
#' @param e_y A vector of expectations in [0,1]
#' @param x A direction of perturbation, where e_y_new = e_y + delta * x
#'
#'@export
get_maximum_delta_for_e_y_perturbation <- function(e_y, x) {
  delta_max_per_obs <- case_when(
    x == 0 ~ Inf,
    x > 0 ~ (1 - e_y) / x,
    x < 0 ~ e_y / abs(x),
    TRUE ~ NA)
  delta_max <- min(delta_max_per_obs)
  return(delta_max)
}


#' Draw binary data with e_y_new, conditional on original binary draws 
#' y_orig and an estimate of their expectation, e_y.
#'
#' @param y_orig A vector of binary {0, 1} draws
#' @param e_y A vector of estimated expectatoins of y_orig
#' @param e_y_new A vector of expectations of the new binary data
#' @param unif_base (Optional) Uniform [0,1] random numbers to condition on.
#'                  If unspecified, new draws are taken.
#'
#' @return A list containing
#'   - ytil: A new draw of binary random variables.  Marginally, ytil has
#'           mean e_y_new if y_orig had mean e_y
#'   - unif_base: The uniform random numbers used for the analysis
#'   - unif_cond: A rescaled version of unif_base to represent a
#'                draw of u | y.
#'
#'@export
draw_conditional_binary_data <- function(y_orig, e_y, e_y_new, unif_base=NULL) {
  if (is.null(unif_base)) {
    unif_base <- runif(n_obs)
  }

  # Draw from the uniform conditional on y_orig and e_y
  rescale_unif <- function(u, lb, ub) {
    return(lb + u * (ub - lb))
  }
  cond_unif <- case_when(
    y_orig == 1 ~ rescale_unif(unif_base, 0, e_y),
    y_orig == 0 ~ rescale_unif(unif_base, e_y, 1),
    TRUE ~ NA
  )
  stopifnot(all(y_orig == as.integer(cond_unif <= e_y)))
  stopifnot(all(e_y_new <= 1))
  stopifnot(all(e_y_new >= 0))
  ytil <- as.integer(cond_unif <= e_y_new)
  return(list(unif_base=unif_base, unif_cond=cond_unif, ytil=ytil))
}

