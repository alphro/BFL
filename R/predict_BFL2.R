#' Posterior predictive sampling — apply pipeline formulation
#'
#' Computes posterior predictive draws using the sapply/apply pipeline.
#' The score for person \eqn{i}, cause \eqn{c}, draw \eqn{s} is:
#'
#' \deqn{
#'   \text{score}[i,c,s]
#'   = \pi_{s,c} \sum_{m=1}^{M} \hat\phi_{i,c,m} \lambda_{s,c,m}
#' }
#'
#' Scores are normalised over causes and a single cause is sampled per
#' person per draw.  The normalisation is applied via
#' \code{apply(..., c(1,3), ...)} which reorders the array dimensions
#' from \eqn{I \times C \times S} to \eqn{C \times I \times S}; the
#' downstream \code{apply(..., c(2,3), sample.int)} and transpose recover
#' the final \eqn{S \times I} draw matrix.
#'
#' @param global_fit Object of class \code{"BFL"} returned by \code{run_BFL()}.
#' @param seed Optional integer seed for reproducible sampling.
#' @export
predict_BFL2 <- function(global_fit, seed = NULL) {

  if (!is.null(seed)) {
    set.seed(seed)
  }

  stopifnot(
    !is.null(global_fit$phi),
    !is.null(global_fit$lambda),
    !is.null(global_fit$pi),
    !is.null(global_fit$causes)
  )

  phi     <- global_fit$phi
  lambda  <- global_fit$lambda
  pi_draw <- global_fit$pi
  causes  <- as.character(global_fit$causes)

  S <- dim(lambda)[1]
  C <- dim(phi)[2]

  # ------------------------------------------------------------------
  # Step 1: build raw score array  (I x C x S)
  #   For each draw s, apply() computes score[i,c] = sum_m phi[i,c,m] *
  #   lambda[s,c,m] for every (i,c) pair, then scales by pi[s,c].
  # ------------------------------------------------------------------
  posterior_pred_Y_prob <- sapply(seq_len(S), function(s) {
    apply(phi, c(1, 2), function(phi_icm, c_idx) {
      sum(phi_icm * lambda[s, c_idx, ])
    }, c_idx = seq_len(C)) * pi_draw[s, ]
  }, simplify = "array")

  # ------------------------------------------------------------------
  # Step 2: normalise over causes for each (i, s)
  #   apply() over c(1,3) passes the length-C slice for each (i,s) pair.
  #   Because the function returns a vector, apply() promotes it to the
  #   first dimension: output shape is C x I x S.
  # ------------------------------------------------------------------
  posterior_pred_Y_prob <- apply(
    posterior_pred_Y_prob,
    c(1, 3),
    function(x) x / sum(x)
  )

  # ------------------------------------------------------------------
  # Step 3: sample one cause per (i, s)
  #   posterior_pred_Y_prob is now C x I x S.
  #   apply() over c(2,3) passes the length-C probability vector for each
  #   (i,s) pair and samples one cause index.  Output shape: I x S.
  #   Transpose to S x I.
  # ------------------------------------------------------------------
  posterior_pred_Y <- apply(
    posterior_pred_Y_prob,
    c(2, 3),
    function(x) sample.int(C, 1L, prob = x)
  )

  draws_int <- t(posterior_pred_Y)

  draws_lab <- matrix(
    causes[draws_int],
    nrow = nrow(draws_int),
    ncol = ncol(draws_int)
  )

  out <- list(
    draws_int = draws_int,
    draws     = draws_lab,
    causes    = causes,
    pi        = pi_draw,
    model     = global_fit$model
  )

  if (!is.null(global_fit$row_hash)) {
    out$row_hash <- global_fit$row_hash
  }

  out
}
