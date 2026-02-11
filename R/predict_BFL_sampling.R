#' Posterior predictive sampling for BFL (paper-faithful)
#'
#' @param global_fit Output of run_BFL()
#' @param ndraws Number of posterior draws to use (subsampled from Stan draws)
#' @param seed Optional seed for reproducibility
#' @param return_probs Also return posterior-mean probs (useful for CSMF)
#'
#' @keywords internal
predict_BFL_sampling <- function(global_fit, ndraws = NULL, seed = NULL, return_probs = TRUE) {
  validate_global_fit(global_fit)

  phi <- global_fit$phi
  d <- dim(phi)
  if (length(d) != 3) stop("global_fit$phi must be an array N x C x M")
  N <- d[1]; C <- d[2]; M <- d[3]

  pi_draws  <- global_fit$pi
  lam_draws <- global_fit$lambda

  S <- nrow(pi_draws)
  if (!is.null(ndraws)) {
    ndraws <- as.integer(ndraws)
    if (ndraws < 1 || ndraws > S) stop("ndraws must be between 1 and nrow(pi)")
    if (!is.null(seed)) set.seed(seed)
    draw_idx <- sample.int(S, ndraws, replace = FALSE)
  } else {
    draw_idx <- seq_len(S)
  }
  S_use <- length(draw_idx)

  pi_use <- pi_draws[draw_idx, , drop = FALSE]
  pi_use <- pi_use / rowSums(pi_use)
  lam_use <- lam_draws[draw_idx, , , drop = FALSE]

  row_normalize <- function(mat) {
    rs <- rowSums(mat)
    rs[rs == 0] <- 1
    mat / rs
  }

  if (!is.null(seed)) set.seed(seed)
  y_draws <- matrix(NA_integer_, nrow = S_use, ncol = N)
  probs_mean <- matrix(0, nrow = N, ncol = C)

  for (k in seq_len(S_use)) {
    score <- matrix(0, nrow = N, ncol = C)
    for (m in seq_len(M)) {
      score <- score + phi[, , m] * matrix(rep(lam_use[k, , m], each = N), nrow = N, ncol = C)
    }
    score <- score * matrix(rep(pi_use[k, ], each = N), nrow = N, ncol = C)

    prob_k <- row_normalize(score)
    probs_mean <- probs_mean + prob_k

    for (i in seq_len(N)) {
      y_draws[k, i] <- sample.int(C, size = 1, prob = prob_k[i, ])
    }
  }

  probs_mean <- probs_mean / S_use
  colnames(probs_mean) <- global_fit$causes

  # majority vote across draws (matches their get_Y_pred spirit)
  pred_idx_vote <- apply(y_draws, 2, function(col) {
    tab <- tabulate(col, nbins = C)
    which.max(tab)
  })
  pred_vote <- global_fit$causes[pred_idx_vote]

  out <- list(pred = pred_vote, draws = y_draws)
  if (return_probs) out$prob <- probs_mean
  out
}
