#' Cause-Specific Mortality Fraction (CSMF) Accuracy
#'
#' Computes CSMF accuracy between estimated and true cause-specific
#' mortality fractions. The metric is defined as:
#'
#' \deqn{
#'   1 - \frac{\sum_c |\hat{\pi}_c - \pi_c|}{2 (1 - \min(\pi))}
#' }
#'
#' where \eqn{\hat{\pi}} is the estimated prevalence vector and
#' \eqn{\pi} is the true prevalence vector.
#'
#' @param pi_hat Numeric vector of estimated class prevalences (length C).
#'   Must sum to 1.
#' @param pi Numeric vector of true class prevalences (length C).
#'   Must sum to 1.
#'
#' @return Numeric scalar giving CSMF accuracy.
#'   Values closer to 1 indicate better agreement.
#'
#' @details
#' Both \code{pi_hat} and \code{pi} must be valid probability
#' distributions of equal length. The function checks that both
#' vectors have the same length and sum approximately to 1.
#'
#' @export
CSMF_acc <- function(pi_hat, pi) {
  # Ensure pi_hat and pi are of the same length and that they are proper probability distributions
  if(length(pi_hat) != length(pi) || sum(pi_hat) > 1.0001 || sum(pi) > 1.0001) {
    stop("Invalid input: pi_hat and pi must be of the same length and sum to 1.")
  }
  
  # Calculate the CSMF accuracy
  C <- length(pi)
  csmf_accuracy <- 1 - (sum(abs(pi_hat - pi)) / (2 * (1 - min(pi))))
  
  return(csmf_accuracy)
}
