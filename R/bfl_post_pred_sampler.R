#' Posterior predictive sampling for global BFL
#'
#' Computes posterior predictive class probabilities and sampled class draws
#' for each observation given posterior draws of \code{lambda} and \code{pi}.
#'
#' Memory-efficient implementation: draws are processed one at a time.
#' Only the running mean (I x C) and integer samples (S x I) are retained,
#' never the full I x C x S probability array.  This keeps peak memory at
#' O(I * C + S * I) rather than O(I * C * S).
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
#'   \item{posterior_pred_Y_prob_mean}{Matrix (I x C) of posterior predictive
#'     class probabilities averaged over all S draws.}
#'   \item{posterior_pred_Y}{Matrix (S x I) of sampled class indices from the
#'     posterior predictive distribution.}
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

  # Accumulators — never allocate I x C x S
  prob_sum      <- matrix(0, nrow = I, ncol = C)   # running sum → mean at end
  draws_int     <- matrix(0L, nrow = S, ncol = I)  # integer class samples

  for (s in seq_len(S)) {
    # score: I x C
    score <- Reduce(`+`, lapply(seq_len(M), function(m) {
      sweep(posterior_phi[, , m], 2, posterior_lambda[s, , m], `*`)
    }))
    score <- sweep(score, 2, posterior_pi[s, ], `*`)
    rs    <- rowSums(score)
    rs[rs == 0] <- 1  # guard against all-zero rows
    score <- score / rs

    prob_sum      <- prob_sum + score
    draws_int[s, ] <- apply(score, 1L, function(p) sample.int(C, 1L, prob = p))
  }

  list(
    posterior_pred_Y_prob_mean = prob_sum / S,  # I x C
    posterior_pred_Y           = draws_int       # S x I
  )
}
