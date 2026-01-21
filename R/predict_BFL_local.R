#' Predict causes locally using a global BFL fit
#'
#' Applies global BFL parameters to local symptom data.
#' This function runs locally and never shares individual-level data.
#'
#' @param global_fit Output of run_BFL_global()
#' @param X Numeric matrix (N x P) of local symptoms
#' @param return_probs Logical; return per-individual probabilities
#'
#' @export
predict_BFL_local <- function(global_fit, X = NULL, return_probs = TRUE) {
  validate_global_fit(global_fit)

  phi <- global_fit$phi
  N <- nrow(phi); C <- ncol(phi)

  if (!is.null(X) && nrow(X) != N) {
    stop("X and phi dimension mismatch.")
  }

  pi <- global_fit$pi
  pi <- pi / sum(pi)

  # unnormalized post probs
  probs <- phi * matrix(rep(pi, each = N), nrow = N, ncol = C)
  probs <- probs / rowSums(probs)

  pred_idx <- max.col(probs, ties.method = "first")
  pred <- global_fit$causes[pred_idx]

  if (return_probs) {
    colnames(probs) <- global_fit$causes
    return(list(pred = pred, prob = probs))
  } else {
    return(list(pred = pred))
  }
}
