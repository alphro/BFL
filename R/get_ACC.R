#' Top-1 Accuracy from BFL posterior draws
#'
#' Computes Top-1 accuracy using the posterior mode across draws.
#'
#' @param pred_Yt Output of \code{predict_BFL()}.
#'   Must contain \code{draws_int} (S x N) and \code{causes}.
#' @param true_Yt True labels (length N). May be character, factor, or integer.
#'
#' @return Numeric scalar giving Top-1 accuracy.
#' @export
get_ACC <- function(pred_Yt, true_Yt) {
  
  stopifnot(!is.null(pred_Yt$draws_int))
  stopifnot(length(true_Yt) == ncol(pred_Yt$draws_int))
  
  draws <- pred_Yt$draws_int
  causes <- pred_Yt$causes
  
  # Map true labels to integer indices 1..C
  cause_to_idx <- setNames(seq_along(causes), causes)
  true_Yt_int <- unname(cause_to_idx[as.character(true_Yt)])
  
  acc <- rep(NA_real_, length(true_Yt_int))
  
  for (t in seq_along(true_Yt_int)) {
    
    if (is.na(true_Yt_int[t])) next
    
    cause_frequencies <- table(draws[, t])
    most_probable_cause_number <- as.integer(
      names(cause_frequencies)[which.max(cause_frequencies)]
    )
    
    acc[t] <- ifelse(most_probable_cause_number == true_Yt_int[t], 1, 0)
  }
  
  mean(acc, na.rm = TRUE)
}