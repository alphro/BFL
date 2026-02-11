#' Deterministic BFL prediction (posterior-mean probabilities)
#'
#' @keywords internal
predict_BFL_deterministic <- function(global_fit, return_probs = TRUE) {
  validate_global_fit(global_fit)

  phi <- global_fit$phi
  d <- dim(phi)
  if (length(d) != 3) stop("global_fit$phi must be an array N x C x M")

  N <- d[1]; C <- d[2]; M <- d[3]

  pi_draws  <- global_fit$pi      # S x C
  lam_draws <- global_fit$lambda  # S x C x M

  if (is.null(dim(pi_draws)) || length(dim(pi_draws)) != 2) {
    stop("global_fit$pi must be a matrix S x C")
  }
  if (is.null(dim(lam_draws)) || length(dim(lam_draws)) != 3) {
    stop("global_fit$lambda must be an array S x C x M")
  }

  S <- nrow(pi_draws)
  pi_draws <- pi_draws / rowSums(pi_draws)

  # posterior-mean predictive probs: E_s[ normalize( pi_s * sum_m lam_s * phi_m ) ]
  probs_mean <- matrix(0, nrow = N, ncol = C)

  row_normalize <- function(mat) {
    rs <- rowSums(mat)
    rs[rs == 0] <- 1
    mat / rs
  }

  for (s in seq_len(S)) {
    score <- matrix(0, nrow = N, ncol = C)
    for (m in seq_len(M)) {
      score <- score + phi[, , m] * matrix(rep(lam_draws[s, , m], each = N), nrow = N, ncol = C)
    }
    score <- score * matrix(rep(pi_draws[s, ], each = N), nrow = N, ncol = C)
    probs_mean <- probs_mean + row_normalize(score)
  }

  probs_mean <- probs_mean / S
  colnames(probs_mean) <- global_fit$causes

  pred_idx <- max.col(probs_mean, ties.method = "first")
  pred <- global_fit$causes[pred_idx]

  if (return_probs) return(list(pred = pred, prob = probs_mean))
  list(pred = pred)
}
