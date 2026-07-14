#' Score a fitted BFL model against ground-truth labels
#'
#' Generates posterior predictive draws, then computes top-1 accuracy,
#' balanced accuracy, CSMF accuracy, a confusion matrix, and per-cause
#' prevalence errors.
#'
#' \strong{Balanced accuracy} is defined as the macro-average recall over
#' causes that actually appear in \code{ytrue_eval} (the evaluation truth
#' labels).  Causes absent from the evaluation set are excluded from the
#' average, matching the convention used by
#' \code{caret::confusionMatrix(...)$overall["Balanced Accuracy"]} and
#' \code{caret::confusionMatrix(...)$byClass[,"Recall"]}.
#'
#' Concretely: a 2-D confusion matrix is built with \strong{Prediction} as
#' rows and \strong{Reference} (truth) as columns, restricted to
#' \code{sort(unique(ytrue_eval))}.  Per-class recall is
#' \code{diag / colSums}; classes with zero reference count receive
#' \code{NA} and are excluded via \code{na.rm = TRUE}.  The full
#' all-cause confusion matrix is retained separately for diagnostics.
#'
#' CSMF computation depends on \code{fit$label_shift}:
#' \itemize{
#'   \item \strong{No shift} (\code{label_shift = FALSE}): \code{pi_true} is the
#'     full-dataset prevalence from \code{Y_eval}. A correction is applied when
#'     held labels were supplied to \code{run_BFL(Y_add = ...)} (Domain/Mix),
#'     using \code{fit$nLc} (cause counts of \code{Y_add}) and \code{fit$n_add}.
#'     Formula: \code{(n_total * pi_hat_raw + nLc) / (n_total + n_add)}.
#'   \item \strong{Label shift} (\code{label_shift = TRUE}): \code{pi_true} is
#'     the prevalence of the \emph{unlabeled evaluation subset} (\code{eval_idx}).
#'     No correction is applied; \code{eval_idx} must be provided.
#' }
#'
#' @param fit Object of class \code{"BFL"} returned by \code{run_BFL()}.
#' @param Y_eval Character or factor vector of true labels (length
#'   \code{fit$n_total}).
#' @param pi_true Optional named numeric vector of true class prevalences. You
#'   normally leave this \code{NULL}: with no shift it defaults to the
#'   full-target prevalence — the model-row truth from \code{Y_eval} plus the
#'   held-row truth from \code{fit$nLc} (the \code{Y_add} counts) — and under
#'   \code{label_shift = TRUE} it is taken from \code{Y_eval[eval_idx]}.
#' @param eval_idx Optional integer indices selecting which rows to score for
#'   top-1 and balanced accuracy. If \code{NULL}, \strong{all} rows of
#'   \code{Y_eval} (every row the model scored) are evaluated. \emph{Required}
#'   when \code{fit$label_shift = TRUE}.
#' @param seed Optional integer seed for reproducible posterior sampling.
#'
#' @return A list with components:
#' \describe{
#'   \item{top1_acc}{Top-1 accuracy on \code{eval_idx} rows.}
#'   \item{balanced_acc}{Balanced accuracy — macro-average recall over causes
#'     observed in the evaluation truth labels.  Matches
#'     \code{caret::confusionMatrix(...)$byClass[,"Recall"]} averaged over
#'     observed classes.}
#'   \item{csmf_acc}{Overall CSMF accuracy (with automatic correction if
#'     Stan saw a subset of rows).}
#'   \item{conf_mat}{Full diagnostic confusion matrix (truth \eqn{\times}
#'     predicted) over all \code{fit$causes}; rows = truth, cols = predicted.}
#'   \item{conf_mat_bal}{Observed-cause confusion matrix used to compute
#'     \code{balanced_acc}; rows = Prediction, cols = Reference (caret
#'     orientation), restricted to \code{sort(unique(ytrue_eval))}.}
#'   \item{pi_hat}{Posterior mean class prevalences after CSMF correction.}
#'   \item{pi_true}{True class prevalences used for evaluation.}
#'   \item{csmf_error_by_cause}{Named numeric: \code{pi_hat - pi_true}
#'     per cause.}
#' }
#'
#' @seealso \code{\link{run_BFL}}, \code{\link{plot_score_BFL}}
#' @export
score_BFL <- function(fit,
                      Y_eval,
                      pi_true  = NULL,
                      eval_idx = NULL,
                      seed     = NULL) {

  stopifnot(inherits(fit, "BFL"))
  stopifnot(!is.null(Y_eval))

  causes <- as.character(fit$causes)
  N      <- fit$n_total
  Y_eval <- as.character(Y_eval)

  if (length(Y_eval) != N) {
    stop("Y_eval must have length equal to fit$n_total (", N, "); got ", length(Y_eval), ".")
  }

  # ------------------------------------------------------------------
  # 1. Posterior predictive draws (internal)
  # ------------------------------------------------------------------
  pred <- predict_BFL(fit, seed = seed)

  # ------------------------------------------------------------------
  # 2. True prevalence + CSMF target
  #
  # No shift / balanced:
  #   target = full dataset prevalence
  #   correction applied when Stan only saw a subset of rows
  #
  # Label shift / unbalanced:
  #   target = unlabeled evaluation subset prevalence
  #   no correction (target of inference is pi on eval_idx)
  # ------------------------------------------------------------------
  pi_hat_raw        <- colMeans(pred$pi)
  names(pi_hat_raw) <- causes
  n_stan <- length(fit$stan_idx)

  if (isTRUE(fit$label_shift)) {
    if (is.null(eval_idx)) {
      stop("For label_shift=TRUE, eval_idx must be provided so CSMF is ",
           "computed on the unlabeled target subset.")
    }
    # true prevalence on evaluation subset only
    if (is.null(pi_true)) {
      tab     <- table(factor(Y_eval[eval_idx], levels = causes))
      pi_true <- as.numeric(tab / sum(tab))
    }
    pi_true_named <- setNames(as.numeric(pi_true), causes)
    # no correction under label shift: target is unlabeled subset
    pi_hat <- pi_hat_raw
  } else {
    # ── No shift: truth + estimate are over the FULL target (model + held rows) ──
    nLc_vec <- setNames(numeric(length(causes)), causes)
    if (!is.null(fit$nLc)) {
      shared          <- intersect(names(fit$nLc), causes)
      nLc_vec[shared] <- as.numeric(fit$nLc[shared])
    }
    n_add <- if (!is.null(fit$n_add)) fit$n_add else sum(nLc_vec)

    # default true prevalence = full target = model-row truth (Y_eval) + held truth (nLc).
    # So Y_add held labels are folded in automatically; you rarely pass pi_true.
    if (is.null(pi_true)) {
      tab     <- as.numeric(table(factor(Y_eval, levels = causes)))
      pi_true <- (tab + nLc_vec) / (length(Y_eval) + n_add)
    }
    pi_true_named <- setNames(as.numeric(pi_true), causes)

    # CSMF correction: fold the held labels into pi_hat (same full-target denominator)
    if (!is.null(fit$nLc)) {
      pi_hat <- (N * pi_hat_raw + nLc_vec) / (N + n_add)
    } else {
      pi_hat <- pi_hat_raw
    }
  }

  csmf_acc <- CSMF_acc(pi_hat, pi_true_named)

  # ------------------------------------------------------------------
  # 4. Top-1 predictions (modal draw per observation)
  # ------------------------------------------------------------------
  # pred$draws_int is S x n_total. fit$stan_idx == 1..N (all rows enter the
  # model), so eval_idx (n_total space) maps to columns directly and the
  # match() below is effectively a no-op; it is kept for robustness.
  if (is.null(eval_idx)) {
    col_idx  <- seq_len(n_stan)   # all Stan rows
    eval_idx <- fit$stan_idx      # corresponding n_total positions for Y_eval
  } else {
    col_idx <- match(eval_idx, fit$stan_idx)
    if (anyNA(col_idx))
      stop("Some eval_idx rows were not seen by Stan (not in fit$stan_idx). ",
           "Top-1 accuracy can only be computed for rows Stan modelled.")
  }

  draws_eval <- pred$draws_int[, col_idx, drop = FALSE]  # S x |eval_idx|
  ytrue_eval <- Y_eval[eval_idx]

  # Modal class index per observation (1-indexed, across S draws)
  yhat_idx  <- modal_class_cpp(draws_eval, length(causes))
  yhat_top1 <- causes[yhat_idx]

  # ------------------------------------------------------------------
  # 5. Top-1 accuracy
  # ------------------------------------------------------------------
  top1_acc <- mean(yhat_top1 == ytrue_eval, na.rm = TRUE)

  # ------------------------------------------------------------------
  # 6. Full diagnostic confusion matrix (all fit$causes)
  #
  # Rows = truth, columns = predicted, levels = all fit$causes.
  # Retained for diagnostics / plotting; NOT used for balanced accuracy.
  # ------------------------------------------------------------------
  conf_mat <- table(
    truth = factor(ytrue_eval, levels = causes),
    pred  = factor(yhat_top1,  levels = causes)
  )

  # ------------------------------------------------------------------
  # 7. Balanced accuracy — macro-average recall over OBSERVED causes
  #
  # Convention: matches caret::confusionMatrix(...)$byClass[,"Recall"].
  #
  # Implementation details:
  #   - Restrict both factor levels to lvls_eval = sort(unique(ytrue_eval))
  #     so only causes that actually appear in the evaluation truth set
  #     contribute to the average.
  #   - Orientation follows caret: Prediction = rows, Reference = cols.
  #     colSums(conf_mat_bal) gives the per-class reference (truth) totals
  #     used as the recall denominator.
  #   - Predictions for causes NOT in lvls_eval become NA in the factor and
  #     are excluded from the confusion matrix, matching caret behaviour.
  #   - Classes with zero reference count receive NA and are excluded via
  #     na.rm = TRUE (should not arise here since lvls_eval = unique(ytrue_eval),
  #     but guarded defensively).
  # ------------------------------------------------------------------
  lvls_eval <- sort(unique(ytrue_eval))

  conf_mat_bal <- table(
    Prediction = factor(yhat_top1,  levels = lvls_eval),
    Reference  = factor(ytrue_eval, levels = lvls_eval)
  )

  ref_totals       <- colSums(conf_mat_bal)           # truth counts per observed cause
  per_class_recall <- diag(conf_mat_bal) / ref_totals # TP / (TP + FN) per cause
  per_class_recall[ref_totals == 0] <- NA             # defensive guard (should not fire)
  balanced_acc     <- mean(per_class_recall, na.rm = TRUE)

  # ------------------------------------------------------------------
  # 8. Per-cause CSMF error
  # ------------------------------------------------------------------
  csmf_error_by_cause <- pi_hat - pi_true_named

  structure(
    list(
      top1_acc            = top1_acc,
      balanced_acc        = balanced_acc,
      csmf_acc            = csmf_acc,
      conf_mat            = conf_mat,      # full C×C diagnostic matrix (truth × pred)
      conf_mat_bal        = conf_mat_bal,  # observed-cause matrix used for balanced_acc (Pred × Ref)
      pi_hat              = pi_hat,
      pi_true             = pi_true_named,
      csmf_error_by_cause = csmf_error_by_cause
    ),
    class = "BFL_score"
  )
}
