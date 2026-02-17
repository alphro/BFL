#' Predict causes from a fitted global BFL model
#'
#' Returns posterior predictive draws for each target observation. Draws are
#' returned both as integer class indices (1..C) and as cause labels.
#'
#' @param global_fit Output of \code{run_BFL()}.
#' @param seed Optional integer seed for reproducible sampling.
#'
#' @return A list with components:
#' \describe{
#'   \item{draws_int}{Matrix (S x N) of integer class draws in 1..C.}
#'   \item{draws}{Matrix (S x N) of cause label draws (character).}
#'   \item{causes}{Character vector of class labels (length C).}
#'   \item{pi}{Posterior draws of class prevalence (S x C).}
#'   \item{model}{User-specified BFL variant.}
#'   \item{row_hash}{Character vector of row hashes (length N), if present.}
#' }
#' @export
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

  out <- list(
    draws_int = draws_int,
    draws     = draws_lab,
    causes    = causes,
    pi        = global_fit$pi,
    model     = global_fit$model
  )

  if (!is.null(global_fit$row_hash)) {
    out$row_hash <- global_fit$row_hash
  }

  out
}
