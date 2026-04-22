#' Score a fitted BFL model against ground-truth labels
#'
#' Generates posterior predictive draws, then computes top-1 accuracy,
#' balanced accuracy, CSMF accuracy, a confusion matrix, and per-cause
#' prevalence errors.
#'
#' CSMF correction is applied automatically when \code{fit$stan_idx}
#' covers fewer rows than \code{fit$n_total} (Domain/Mix variants).
#' The correction uses \code{fit$nLc} — labeled-row counts for rows
#' outside Stan — which are stored by \code{run_BFL()} at fit time.
#'
#' @param fit Object of class \code{"BFL"} returned by \code{run_BFL()}.
#' @param Y_eval Character or factor vector of true labels (length
#'   \code{fit$n_total}).
#' @param pi_true Optional named numeric vector of true class prevalences.
#'   If \code{NULL}, computed empirically from \code{Y_eval}.
#' @param eval_idx Optional integer indices defining the evaluation subset for
#'   top-1 and balanced accuracy (e.g. held-out unlabeled rows).
#'   If \code{NULL}, all \code{N} rows are used.
#' @param seed Optional integer seed for reproducible posterior sampling.
#'
#' @return A list with components:
#' \describe{
#'   \item{top1_acc}{Top-1 accuracy on \code{eval_idx} rows.}
#'   \item{balanced_acc}{Balanced accuracy (mean per-class recall) on
#'     \code{eval_idx} rows.}
#'   \item{csmf_acc}{Overall CSMF accuracy (with automatic correction if
#'     Stan saw a subset of rows).}
#'   \item{conf_mat}{Confusion matrix (truth \eqn{\times} predicted).}
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
  # 2. True prevalence
  # ------------------------------------------------------------------
  if (is.null(pi_true)) {
    tab     <- table(factor(Y_eval, levels = causes))
    pi_true <- as.numeric(tab / sum(tab))
  }
  pi_true_named <- setNames(as.numeric(pi_true), causes)

  # ------------------------------------------------------------------
  # 3. CSMF correction
  #    Correction is needed when Stan only saw a subset of N rows
  #    (Domain/Mix). fit$nLc holds labeled-row counts outside Stan.
  #    Formula: (n_stan * colMeans(pi) + nLc) / n_total
  # ------------------------------------------------------------------
  pi_hat_raw        <- colMeans(pred$pi)   # posterior mean from Stan
  names(pi_hat_raw) <- causes

  n_stan           <- length(fit$stan_idx)
  needs_correction <- n_stan < N

  if (needs_correction) {
    nLc_vec           <- numeric(length(causes))
    names(nLc_vec)    <- causes
    if (!is.null(fit$nLc)) {
      shared              <- intersect(names(fit$nLc), causes)
      nLc_vec[shared]     <- as.numeric(fit$nLc[shared])
    }
    pi_hat <- (n_stan * pi_hat_raw + nLc_vec) / N
  } else {
    pi_hat <- pi_hat_raw
  }

  csmf_acc <- CSMF_acc(pi_hat, pi_true_named)

  # ------------------------------------------------------------------
  # 4. Top-1 predictions (modal draw per observation)
  # ------------------------------------------------------------------
  # pred$draws_int is S x n_stan (Stan-row space, NOT n_total space).
  # eval_idx is expressed in n_total (1..N) space, so we must map it
  # to column positions in draws_int via fit$stan_idx.
  # For Partial/Base (n_stan == N), fit$stan_idx == 1..N so match() is a no-op.
  # For Domain/Mix (n_stan < N), this converts e.g. missing_indices → 1..n_stan cols.
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
  # 6. Confusion matrix
  # ------------------------------------------------------------------
  conf_mat <- table(
    truth = factor(ytrue_eval, levels = causes),
    pred  = factor(yhat_top1,  levels = causes)
  )

  # ------------------------------------------------------------------
  # 7. Balanced accuracy (mean per-class recall)
  # ------------------------------------------------------------------
  per_class_recall <- diag(conf_mat) / pmax(rowSums(conf_mat), 1L)
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
      conf_mat            = conf_mat,
      pi_hat              = pi_hat,
      pi_true             = pi_true_named,
      csmf_error_by_cause = csmf_error_by_cause
    ),
    class = "BFL_score"
  )
}
