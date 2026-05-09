score_BFL2 <- function(fit,
                       Y_eval,
                       eval_idx = NULL,
                       seed = NULL) {

  pred <- predict_BFL2(fit, seed = seed)

  causes <- as.character(fit$causes)

  if (is.null(eval_idx)) {
    col_idx  <- seq_len(ncol(pred$draws_int))
    ytrue_eval <- as.character(Y_eval)
  } else {
    col_idx <- match(eval_idx, fit$stan_idx)
    ytrue_eval <- as.character(Y_eval[eval_idx])
  }

  draws_eval <- pred$draws_int[, col_idx, drop = FALSE]

  yhat_idx <- apply(draws_eval, 2, function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  })

  yhat_top1 <- causes[yhat_idx]

  conf_mat <- table(
    truth = factor(ytrue_eval, levels = causes),
    pred  = factor(yhat_top1,  levels = causes)
  )

  row_totals       <- rowSums(conf_mat)
  per_class_recall <- diag(conf_mat) / row_totals
  per_class_recall[row_totals == 0] <- NA

  list(
    top1_acc     = mean(yhat_top1 == ytrue_eval),
    balanced_acc = mean(per_class_recall, na.rm = TRUE),
    conf_mat     = conf_mat
  )
}
