#' Scientist-facing summary for BFL predictions
#'
#' @param pred_out Output of predict_BFL(return_probs = TRUE)
#' @param y_true Character vector of true causes (length N)
#' @param csmf_fun Function(pi_hat, pi_true) -> numeric
#' @param cause_level Logical; if TRUE, compute per-cause table (slow)
#'
#' @export
summary_BFL <- function(pred_out,
                        y_true,
                        csmf_fun = csmf_accuracy,
                        cause_level = FALSE) {

  if (is.null(pred_out$pred)) stop("pred_out must contain $pred.")
  if (is.null(pred_out$prob)) stop("pred_out must contain $prob (use return_probs=TRUE).")

  pred <- as.character(pred_out$pred)
  prob <- pred_out$prob
  y_true <- as.character(y_true)

  if (length(pred) != length(y_true))
    stop("pred and y_true must have same length.")

  N <- length(y_true)

  # ---- population-level metrics (FAST) ----
  unique_causes_in_truth <- length(unique(y_true))
  top1 <- mean(pred == y_true, na.rm = TRUE)

  pi_hat <- colMeans(prob, na.rm = TRUE)
  pi_hat <- pi_hat / sum(pi_hat)

  pi_true <- prop.table(table(y_true))
  csmf <- csmf_fun(pi_hat, pi_true)

  population <- tibble::tibble(
    N = N,
    unique_causes_in_truth = unique_causes_in_truth,
    csmf_acc = csmf,
    top1_acc = top1
  )

  # ---- EARLY EXIT (default) ----
  if (!cause_level) {
    return(list(
      population = population,
      per_cause = NULL
    ))
  }

  # ---- per-cause table (SLOW) ----
  causes <- sort(unique(c(y_true, pred)))

  per_cause <- lapply(causes, function(ca) {
    tp <- sum(pred == ca & y_true == ca, na.rm = TRUE)
    fp <- sum(pred == ca & y_true != ca, na.rm = TRUE)
    fn <- sum(pred != ca & y_true == ca, na.rm = TRUE)

    support <- sum(y_true == ca, na.rm = TRUE)
    pred_count <- sum(pred == ca, na.rm = TRUE)

    recall <- if (support == 0) NA_real_ else tp / support
    precision <- if (pred_count == 0) NA_real_ else tp / pred_count
    f1 <- if (is.na(recall) || is.na(precision) || (recall + precision) == 0)
      NA_real_ else 2 * recall * precision / (recall + precision)

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

  list(
    population = population,
    per_cause = per_cause
  )
}
