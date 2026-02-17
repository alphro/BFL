#' Summary for BFL predictions
#'
#' Computes top-1 accuracy and CSMF accuracy. If `pred_out$draws` is present
#' (sampling mode), top-1 uses majority vote across draws (paper style).
#'
#' @param pred_out Output of predict_BFL(return_probs = TRUE). Must contain
#'   `pred` and `prob`. May also contain `draws` (S x N integer indices).
#' @param y_true Character vector of true causes (length N).
#' @param cause_level Logical; if TRUE, compute per-cause table (slow).
#'
#' @export
summary_BFL <- function(pred_out, y_true, cause_level = FALSE) {
  if (is.null(pred_out$pred)) stop("pred_out must contain $pred.")
  if (is.null(pred_out$prob)) stop("pred_out must contain $prob (use return_probs=TRUE).")

  pred   <- as.character(pred_out$pred)
  prob   <- pred_out$prob
  y_true <- as.character(y_true)

  if (length(pred) != length(y_true)) stop("pred and y_true must have same length.")
  N <- length(y_true)

  causes_prob <- colnames(prob)
  if (is.null(causes_prob)) stop("pred_out$prob must have colnames = cause labels.")

  # ----------------------------
  # Top-1 accuracy (paper style if draws exist)
  # ----------------------------
  if (!is.null(pred_out$draws)) {
    draws <- pred_out$draws  # S x N integer indices into causes_prob
    if (!is.matrix(draws)) stop("pred_out$draws must be a matrix (S x N).")
    if (ncol(draws) != N) stop("pred_out$draws must have ncol == length(y_true).")

    # majority vote per individual (matches get_Y_pred/get_acc spirit)
    pred_idx_vote <- apply(draws, 2, function(col) which.max(tabulate(col, nbins = length(causes_prob))))
    pred_vote <- causes_prob[pred_idx_vote]
    top1 <- mean(pred_vote == y_true, na.rm = TRUE)
  } else {
    # deterministic path: pred is already final
    top1 <- mean(pred == y_true, na.rm = TRUE)
  }

  # ----------------------------
  # CSMF accuracy (match paper's CSMF_acc, but aligned)
  # ----------------------------
  # pi_hat over *all* causes in prob columns
  pi_hat <- colMeans(prob, na.rm = TRUE)
  pi_hat <- pi_hat / sum(pi_hat)

  # pi_true aligned to same cause set (missing causes => 0)
  tab_true <- table(factor(y_true, levels = causes_prob))
  pi_true <- as.numeric(tab_true) / sum(tab_true)

  # paper formula + checks
  if (length(pi_hat) != length(pi_true) || sum(pi_hat) > 1.0001 || sum(pi_true) > 1.0001) {
    stop("Invalid pi_hat/pi_true (length mismatch or not summing to 1).")
  }
  csmf <- 1 - (sum(abs(pi_hat - pi_true)) / (2 * (1 - min(pi_true))))

  population <- tibble::tibble(
    N = N,
    unique_causes_in_truth = length(unique(y_true)),
    csmf_acc = csmf,
    top1_acc = top1
  )

  if (!cause_level) return(list(population = population, per_cause = NULL))

  # ----------------------------
  # per-cause table (unchanged)
  # ----------------------------
  causes <- sort(unique(c(y_true, pred)))
  per_cause <- lapply(causes, function(ca) {
    tp <- sum(pred == ca & y_true == ca, na.rm = TRUE)
    fp <- sum(pred == ca & y_true != ca, na.rm = TRUE)
    fn <- sum(pred != ca & y_true == ca, na.rm = TRUE)

    support <- sum(y_true == ca, na.rm = TRUE)
    pred_count <- sum(pred == ca, na.rm = TRUE)

    recall <- if (support == 0) NA_real_ else tp / support
    precision <- if (pred_count == 0) NA_real_ else tp / pred_count
    f1 <- if (is.na(recall) || is.na(precision) || (recall + precision) == 0) NA_real_
    else 2 * recall * precision / (recall + precision)

    tibble::tibble(
      cause = ca,
      support = support,
      prevalence = support / N,
      pred_count = pred_count,
      recall = recall,
      precision = precision,
      f1 = f1
    )
  }) |>
    dplyr::bind_rows() |>
    dplyr::arrange(dplyr::desc(prevalence))
a
  list(population = population, per_cause = per_cause)
}
