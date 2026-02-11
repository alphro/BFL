#' Predict causes using a global BFL fit
#'
#' @param global_fit Output of run_BFL()
#' @param sampling Logical; if TRUE, use posterior predictive sampling.
#' @param ndraws Number of posterior draws to use when sampling (NULL = all).
#' @param seed Optional seed for reproducibility (only used when sampling = TRUE or ndraws not NULL).
#' @param return_probs Logical; return posterior-mean probabilities.
#'
#' @export
predict_BFL <- function(global_fit,
                        sampling = FALSE,
                        ndraws = NULL,
                        seed = NULL,
                        return_probs = TRUE) {
  if (isTRUE(sampling)) {
    predict_BFL_sampling(global_fit, ndraws = ndraws, seed = seed, return_probs = return_probs)
  } else {
    predict_BFL_deterministic(global_fit, return_probs = return_probs)
  }
}
