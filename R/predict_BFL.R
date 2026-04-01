#' Posterior predictive draws from a fitted BFL model (internal)
#'
#' Called internally by \code{score_BFL()}. Users should use
#' \code{score_BFL()} directly rather than calling this function.
#'
#' @param global_fit Object of class \code{"BFL"} returned by \code{run_BFL()}.
#' @param seed Optional integer seed for reproducible sampling.
#'
#' @return A list with components:
#' \describe{
#'   \item{draws_int}{Matrix (S x N) of integer class draws in 1..C.}
#'   \item{draws}{Matrix (S x N) of cause-label draws (character).}
#'   \item{causes}{Character vector of class labels (length C).}
#'   \item{pi}{Posterior draws of class prevalence (S x C).}
#'   \item{row_hash}{Character vector of row hashes (length N).}
#' }
#' @noRd
predict_BFL <- function(global_fit, seed = NULL) {

  stopifnot(
    !is.null(global_fit$phi),
    !is.null(global_fit$lambda),
    !is.null(global_fit$pi),
    !is.null(global_fit$causes)
  )

  pred_res <- bfl_post_pred_sampler(
    posterior_phi    = global_fit$phi,
    posterior_lambda = global_fit$lambda,
    posterior_pi     = global_fit$pi,
    seed             = seed
  )

  draws_int <- pred_res$posterior_pred_Y  # S x N
  causes    <- as.character(global_fit$causes)

  draws_lab <- matrix(
    causes[draws_int],
    nrow = nrow(draws_int),
    ncol = ncol(draws_int)
  )

  list(
    draws_int = draws_int,
    draws     = draws_lab,
    causes    = causes,
    pi        = global_fit$pi,
    row_hash  = global_fit$row_hash
  )
}
