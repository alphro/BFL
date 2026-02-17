#' Score predictions from a fitted BFL model
#'
#' Computes Top-1 accuracy, CSMF accuracy, a confusion matrix, and
#' per-cause prevalence differences for a set of BFL predictions.
#'
#' Evaluation can be performed either:
#' \itemize{
#'   \item Using a \code{mask} object (recommended; produced by a masking helper), or
#'   \item By directly supplying \code{Y_eval} (and optionally \code{missing_idx}, \code{stan_idx}).
#' }
#'
#' For \code{"Domain"} and \code{"Mix"} models, \code{stan_idx} is required
#' to reconstruct full-sample CSMF estimates.
#'
#' @param pred Output from \code{predict_BFL()}.
#' @param mask Optional mask object containing:
#'   \itemize{
#'     \item \code{Y_true}
#'     \item \code{missing_idx}
#'     \item \code{stan_idx}
#'   }
#' @param Y_eval Optional vector of true labels (length N).
#' @param pi_true Optional numeric vector of true class prevalences.
#'   If NULL, computed empirically from \code{Y_eval}.
#' @param missing_idx Optional integer indices defining evaluation subset
#'   (typically unlabeled rows).
#' @param stan_idx Optional integer indices used by Stan (required for
#'   \code{"Domain"} and \code{"Mix"} models if not using \code{mask}).
#'
#' @return A list with components:
#' \describe{
#'   \item{top1_acc}{Top-1 classification accuracy.}
#'   \item{csmf_acc}{Overall CSMF accuracy (global prevalence metric).}
#'   \item{conf_mat}{Confusion matrix (truth x predicted).}
#'   \item{pi_hat}{Posterior mean class prevalence.}
#'   \item{pi_true}{True class prevalence used for evaluation.}
#'   \item{csmf_error_by_cause}{Named numeric vector of per-cause
#'     prevalence differences (\code{pi_hat - pi_true}).}
#' }
#'
#' @export
score_BFL <- function(pred,
                      mask = NULL,
                      Y_eval = NULL,
                      pi_true = NULL,
                      missing_idx = NULL,
                      stan_idx = NULL) {
  
  stopifnot(!is.null(pred$model))
  stopifnot(!is.null(pred$draws_int))
  stopifnot(!is.null(pred$pi))
  stopifnot(!is.null(pred$causes))
  
  # ---------------------------------------
  # Resolve evaluation labels + indices
  # ---------------------------------------
  if (!is.null(mask)) {
    Y_eval      <- mask$Y_true
    missing_idx <- mask$missing_idx
    stan_idx    <- mask$stan_idx
  }
  
  if (is.null(Y_eval)) stop("Must supply either mask or Y_eval.")
  
  causes <- as.character(pred$causes)
  N <- length(Y_eval)
  
  # ---------------------------------------
  # Compute pi_true if not supplied
  # ---------------------------------------
  if (is.null(pi_true)) {
    tab <- table(factor(as.character(Y_eval), levels = causes))
    pi_true <- as.numeric(tab / sum(tab))
  }
  pi_true_named <- pi_true
  names(pi_true_named) <- causes
  
  # posterior mean prevalence
  pi_hat <- colMeans(pred$pi)
  names(pi_hat) <- causes
  
  # ---------------------------------------
  # Top-1 accuracy (+ confusion matrix)
  # Evaluate on missing_idx if provided; otherwise on all rows
  # ---------------------------------------
  if (is.null(missing_idx)) {
    eval_idx <- seq_len(N)
    top1_acc <- get_ACC(pred, Y_eval)
    # top-1 predictions as labels
    yhat_top1 <- causes[apply(pred$draws_int, 2, function(x) {
      tab <- table(x)
      as.integer(names(tab)[which.max(tab)])
    })]
  } else {
    eval_idx <- missing_idx
    
    pred_missing <- list(
      draws_int = pred$draws_int[, eval_idx, drop = FALSE],
      causes    = causes
    )
    top1_acc <- get_ACC(pred_missing, Y_eval[eval_idx])
    
    yhat_top1 <- causes[apply(pred_missing$draws_int, 2, function(x) {
      tab <- table(x)
      as.integer(names(tab)[which.max(tab)])
    })]
  }
  
  ytrue_eval <- as.character(Y_eval[eval_idx])
  conf_mat <- table(
    truth = factor(ytrue_eval, levels = causes),
    pred  = factor(yhat_top1,  levels = causes)
  )
  
  # ---------------------------------------
  # CSMF logic (global metric)
  # ---------------------------------------
  if (pred$model %in% c("Partial", "Base")) {
    csmf_acc <- CSMF_acc(pi_hat, pi_true_named)
  } else {
    if (is.null(stan_idx)) stop("stan_idx must be supplied for Domain or Mix.")
    
    n0 <- N
    nU <- length(stan_idx)
    
    not_stan_idx <- setdiff(seq_len(n0), stan_idx)
    tab_labeled <- table(factor(as.character(Y_eval[not_stan_idx]), levels = causes))
    nLc <- as.numeric(tab_labeled)
    
    pi_hat_full <- (nU * pi_hat + nLc) / n0
    csmf_acc <- CSMF_acc(pi_hat_full, pi_true_named)
  }
  
  # ---------------------------------------
  # Per-cause prevalence error (NOT "CSMF per cause")
  # ---------------------------------------
  csmf_error_by_cause <- pi_hat - pi_true_named
  
  list(
    top1_acc = top1_acc,
    csmf_acc = csmf_acc,
    conf_mat = conf_mat,
    pi_hat = pi_hat,
    pi_true = pi_true_named,
    csmf_error_by_cause = csmf_error_by_cause
  )
}