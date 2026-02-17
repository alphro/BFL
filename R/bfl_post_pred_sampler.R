#' Posterior predictive sampling for global BFL
#'
#' Computes posterior predictive class probabilities and sampled class draws
#' for each observation given posterior draws of \code{lambda} and \code{pi}.
#'
#' @param posterior_phi Array (I x C x M) of per-observation class scores/likelihoods
#'   from local models, aligned to global classes.
#' @param posterior_lambda Array (S x C x M) of posterior draws for per-class model
#'   weights.
#' @param posterior_pi Matrix (S x C) of posterior draws for class prevalence.
#' @param seed Optional integer seed for reproducible sampling.
#'
#' @return A list with components:
#' \describe{
#'   \item{posterior_pred_Y_prob}{Array (I x C x S) of posterior predictive class
#'   probabilities.}
#'   \item{posterior_pred_Y}{Matrix (S x I) of sampled class indices from the
#'   posterior predictive distribution.}
#' }
#'
#' @keywords internal
bfl_post_pred_sampler <- function(
    posterior_phi,
    posterior_lambda,
    posterior_pi,
    seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)

  # Dimensions
  S <- dim(posterior_lambda)[1]
  I <- dim(posterior_phi)[1]
  C <- dim(posterior_phi)[2]

  # Compute unnormalized posterior predictive probs for each draw s:
  # p_s(i,c) ∝ pi_s(c) * sum_m lambda_s(c,m) * phi(i,c,m)
  posterior_pred_Y_prob <- sapply(seq_len(S), function(s) {
    apply(posterior_phi, c(1, 2), function(phi_icm) {
      sum(phi_icm * posterior_lambda[s, , ])
    }) * posterior_pi[s, ]
  }, simplify = "array")  # I x C x S

  # Normalize across classes for each (i, s)
  posterior_pred_Y_prob <- apply(
    posterior_pred_Y_prob,
    c(1, 3),
    function(x) x / sum(x)
  ) # still I x C x S

  # Sample Y for each (i, s)
  posterior_pred_Y <- apply(
    posterior_pred_Y_prob,
    c(2, 3),
    function(x) sample.int(C, 1, prob = x)
  ) # I x S

  list(
    posterior_pred_Y_prob = posterior_pred_Y_prob,
    posterior_pred_Y      = t(posterior_pred_Y) # S x I
  )
}
