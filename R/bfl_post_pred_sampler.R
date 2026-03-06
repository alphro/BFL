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
bfl_post_pred_sampler <- function(posterior_phi, posterior_lambda, posterior_pi, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  S <- dim(posterior_lambda)[1]
  I <- dim(posterior_phi)[1]
  C <- dim(posterior_phi)[2]
  M <- dim(posterior_phi)[3]

  stopifnot(
    length(dim(posterior_phi)) == 3,
    length(dim(posterior_lambda)) == 3,
    nrow(posterior_pi) == S,
    ncol(posterior_pi) == C,
    dim(posterior_lambda)[2] == C,
    dim(posterior_lambda)[3] == M
  )

  posterior_pred_Y_prob <- sapply(seq_len(S), function(s) {
    # score: I x C
    score <- Reduce(`+`, lapply(seq_len(M), function(m) {
      sweep(posterior_phi[, , m], 2, posterior_lambda[s, , m], `*`)
    }))
    score <- sweep(score, 2, posterior_pi[s, ], `*`)
    score <- score / rowSums(score)
    score
  }, simplify = "array")  # I x C x S

  # sample class for each (i,s) -> I x S, then transpose to S x I
  posterior_pred_Y <- apply(posterior_pred_Y_prob, c(1,3), function(p) {
    sample.int(C, 1, prob = p)
  })

  list(
    posterior_pred_Y_prob = posterior_pred_Y_prob,   # I x C x S
    posterior_pred_Y      = t(posterior_pred_Y)      # S x I
  )
}
