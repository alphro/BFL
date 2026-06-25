#' Score a BFL fit using the apply-pipeline prediction formula
#'
#' Calls \code{predict_BFL2()} and computes top-1 accuracy and balanced
#' accuracy.  The modal cause across posterior draws is taken as the point
#' prediction for each person.  Balanced accuracy is macro-average recall
#' computed over all causes in \code{fit$causes}, with causes absent from
#' the truth labels contributing \code{NA} (excluded from the mean).
#'
#' @param fit      Object of class \code{"BFL"} returned by \code{run_BFL()}.
#' @param Y_eval   Character or factor vector of true labels.
#' @param eval_idx Optional integer indices selecting the evaluation subset
#'   from the columns of the draw matrix.
#' @param seed     Optional integer seed for reproducible sampling.
#'
#' @return List with elements \code{top1_acc} and \code{balanced_acc}.
#' @export
score_BFL2 <- function(fit,
                       Y_eval,
                       eval_idx = NULL,
                       seed = NULL) {

  pred   <- predict_BFL2(fit, seed = seed)
  causes <- as.character(fit$causes)

  if (is.null(eval_idx)) {
    col_idx    <- seq_len(ncol(pred$draws_int))
    ytrue_eval <- as.character(Y_eval)
  } else {
    col_idx    <- match(eval_idx, fit$stan_idx)
    ytrue_eval <- as.character(Y_eval[eval_idx])
  }

  draws_eval <- pred$draws_int[, col_idx, drop = FALSE]

  # Modal cause across draws for each person
  yhat_idx <- apply(draws_eval, 2, function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  })
  yhat <- causes[yhat_idx]

  # Confusion matrix: rows = truth, cols = predicted
  conf_mat <- table(
    truth = factor(ytrue_eval, levels = causes),
    pred  = factor(yhat,       levels = causes)
  )

  # Balanced accuracy: macro-average recall (per-row recall, row = true class)
  row_totals <- rowSums(conf_mat)
  recall     <- diag(conf_mat) / row_totals
  recall[row_totals == 0] <- NA

  list(
    top1_acc     = mean(yhat == ytrue_eval),
    balanced_acc = mean(recall, na.rm = TRUE),
    conf_mat     = conf_mat
  )
}
