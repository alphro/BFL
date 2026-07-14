#' Posterior predictive sampling for global BFL
#'
#' Computes posterior predictive class probabilities and sampled class draws
#' for each observation given posterior draws of \code{lambda} and \code{pi}.
#'
#' @param posterior_phi Array (I x C x M) of per-observation class scores.
#' @param posterior_lambda Array (S x C x M) of posterior model-weight draws.
#' @param posterior_pi Matrix (S x C) of posterior cause-fraction draws.
#' @param seed Optional integer seed for reproducible sampling.
#'
#' @return A list with:
#' \describe{
#'   \item{posterior_pred_Y_prob_mean}{Matrix (I x C) of posterior predictive
#'     class probabilities averaged over all S draws.}
#'   \item{posterior_pred_Y}{Matrix (S x I) of sampled class indices (1..C).}
#' }
#' @keywords internal
bfl_post_pred_sampler <- function(posterior_phi, posterior_lambda, posterior_pi,
                                  seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  I <- dim(posterior_phi)[1]
  C <- dim(posterior_phi)[2]
  M <- dim(posterior_phi)[3]
  S <- dim(posterior_lambda)[1]

  post_pred_sample_cpp(
    phi_vec    = as.vector(posterior_phi),
    lambda_vec = as.vector(posterior_lambda),
    pi_vec     = as.vector(posterior_pi),
    I = I, C = C, M = M, S = S
  )
}
