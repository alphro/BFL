#' Predict causes locally using a global BFL fit
#'
#' Applies global BFL parameters to local symptom data.
#' This function runs locally and never shares individual-level data.
#'
#' @param global_fit Output of run_BFL_global()
#' @param return_probs Logical; return per-individual probabilities
#'
#' @export
predict_BFL_local <- function(global_fit, return_probs = TRUE) {
  validate_global_fit(global_fit)

  phi <- global_fit$phi
  N <- nrow(phi); C <- ncol(phi)

  pi <- global_fit$pi
  pi <- pi / sum(pi)

  # log posterior up to constant
  logp <- log(pmax(phi, 1e-300)) + matrix(log(pi), nrow = N, ncol = C, byrow = TRUE)

  # stable row softmax
  logp <- logp - apply(logp, 1, max)
  probs <- exp(logp)
  probs <- probs / rowSums(probs)

  pred_idx <- max.col(probs, ties.method = "first")
  pred <- global_fit$causes[pred_idx]

  if (return_probs) {
    colnames(probs) <- global_fit$causes
    list(pred = pred, prob = probs)
  } else {
    list(pred = pred)
  }
}
