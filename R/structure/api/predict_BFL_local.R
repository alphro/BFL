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
predict_BFL_local <- function(
    global_fit,
    X,
    return_probs = FALSE
) {

  validate_global_fit(global_fit)

  X <- as.matrix(X)
  storage.mode(X) <- "numeric"

  log_pi <- safe_log(global_fit$pi)

  # phi must be provided via global_fit
  # expected shape: P x C
  phi <- global_fit$phi

  if (ncol(X) != nrow(phi)) {
    stop("X and phi dimension mismatch.")
  }

  # Compute log posterior up to normalization
  log_post <- X %*% safe_log(phi)
  log_post <- sweep(log_post, 2, log_pi, "+")

  # Normalize
  probs <- apply(log_post, 1, function(x) {
    softmax(x)
  })
  probs <- t(probs)

  y_pred <- apply(probs, 1, argmax)
  csmf_hat <- colMeans(probs)

  out <- list(
    y_pred = y_pred,
    csmf_hat = csmf_hat
  )

  if (return_probs) {
    out$probs <- probs
  }

  out
}
