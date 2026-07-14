#' S3 methods and diagnostics for BFL objects
#'
#' Methods for objects of class \code{"BFL"} (returned by \code{\link{run_BFL}})
#' and class \code{"BFL_score"} (returned by \code{\link{score_BFL}}).
#'
#' \strong{BFL fit methods:}
#' \describe{
#'   \item{\code{print.BFL}}{Compact one-glance overview of the fit.}
#'   \item{\code{summary.BFL}}{Full CSMF table with credible intervals and
#'     posterior mean lambda weights.}
#'   \item{\code{plot.BFL}}{CSMF dot-plot with credible intervals.}
#' }
#'
#' \strong{BFL score methods:}
#' \describe{
#'   \item{\code{print.BFL_score}}{Quick three-number summary: top-1,
#'     balanced accuracy, CSMF accuracy.}
#'   \item{\code{summary.BFL_score}}{Full per-cause table: true prevalence,
#'     estimated prevalence, CSMF error, and recall per cause.}
#'   \item{\code{plot.BFL_score}}{Two panels: per-cause recall bars and
#'     per-cause CSMF error bars.}
#' }
#'
#' @name BFL_methods
#' @aliases print.BFL summary.BFL plot.BFL
#'   print.BFL_score summary.BFL_score plot.BFL_score
NULL


# ============================================================
# print.BFL — compact one-glance overview
# ============================================================

#' @rdname BFL_methods
#' @param x A \code{"BFL"} object.
#' @param ... Ignored.
#' @export
print.BFL <- function(x, ...) {

  C <- length(x$causes)
  M <- length(x$model_names)
  N <- x$n_total
  S <- nrow(x$pi)

  needs_correction <- !is.null(x$nLc)
  variant <- if (!x$has_labels) {
    if (needs_correction) "Domain (no partial labels, subset to Stan)" else "Base (no partial labels)"
  } else {
    if (x$label_shift) {
      if (needs_correction) "Mix / label-shift" else "Partial / label-shift"
    } else {
      if (needs_correction) "Mix" else "Partial"
    }
  }

  cat("BFL Fit\n")
  cat(strrep("\u2500", 40), "\n", sep = "")
  cat(sprintf("  Records   : %d\n", N))
  cat(sprintf("  Causes    : %d\n", C))
  cat(sprintf("  Models    : %d (%s)\n", M, paste(x$model_names, collapse = ", ")))
  cat(sprintf("  Post draws: %d\n", S))
  cat(sprintf("  Variant   : %s\n", variant))
  cat(strrep("\u2500", 40), "\n", sep = "")

  # Top-3 causes by posterior mean prevalence
  pi_mean  <- colMeans(x$pi)
  names(pi_mean) <- x$causes
  top3     <- sort(pi_mean, decreasing = TRUE)[seq_len(min(3, C))]
  top3_str <- paste(sprintf("%s (%.1f%%)", names(top3), top3 * 100), collapse = ", ")
  cat(sprintf("  Top causes: %s\n", top3_str))
  cat("  Run summary() for full CSMF table.\n")

  invisible(x)
}


# ============================================================
# summary.BFL — full CSMF table + lambda weights
# ============================================================

#' @rdname BFL_methods
#' @param object A \code{"BFL"} object.
#' @param ci_prob Credible interval width. Default \code{0.95}.
#' @param ... Ignored.
#' @export
summary.BFL <- function(object, ci_prob = 0.95, ...) {

  x   <- object
  C   <- length(x$causes)
  M   <- length(x$model_names)
  N   <- x$n_total
  S   <- nrow(x$pi)

  alpha     <- (1 - ci_prob) / 2
  pi_mean   <- colMeans(x$pi)
  pi_lower  <- apply(x$pi, 2, quantile, probs = alpha)
  pi_upper  <- apply(x$pi, 2, quantile, probs = 1 - alpha)

  names(pi_mean) <- names(pi_lower) <- names(pi_upper) <- x$causes

  # Sort by mean prevalence descending
  ord      <- order(pi_mean, decreasing = TRUE)
  pi_mean  <- pi_mean[ord]
  pi_lower <- pi_lower[ord]
  pi_upper <- pi_upper[ord]

  # Lambda mean per cause x model
  lambda_mean <- apply(x$lambda, c(2, 3), mean)  # C x M
  rownames(lambda_mean) <- x$causes
  colnames(lambda_mean) <- x$model_names
  lambda_mean <- lambda_mean[ord, , drop = FALSE]

  # Inline bar (max width 12 chars)
  .bar <- function(p, width = 12) {
    n <- max(0L, round(p * width))
    paste0(strrep("\u2588", n), strrep("\u2591", width - n))
  }

  cat("BFL Results\n")
  cat(strrep("\u2500", 60), "\n", sep = "")
  cat(sprintf("  Records: %d  |  Models: %d  |  Causes: %d  |  Draws: %d\n",
              N, M, C, S))
  cat(strrep("\u2500", 60), "\n", sep = "")

  cat(sprintf("\nCause-Specific Mortality Fractions (%d%% CI):\n",
              round(ci_prob * 100)))
  cat(sprintf("  %-22s %6s  %6s  %6s\n", "Cause", "Mean", "Lower", "Upper"))
  cat(sprintf("  %s\n", strrep("-", 55)))

  for (i in seq_along(pi_mean)) {
    cat(sprintf("  %-22s %5.1f%%  %5.1f%%  %5.1f%%  %s\n",
                names(pi_mean)[i],
                pi_mean[i]  * 100,
                pi_lower[i] * 100,
                pi_upper[i] * 100,
                .bar(pi_mean[i])))
  }

  cat(sprintf("\nModel Weights (posterior mean lambda):\n"))
  header <- sprintf("  %-22s", "Cause")
  for (m in x$model_names) header <- paste0(header, sprintf("  %8s", m))
  cat(header, "\n")
  cat(sprintf("  %s\n", strrep("-", 22 + 10 * M)))

  for (i in seq_len(nrow(lambda_mean))) {
    row <- sprintf("  %-22s", rownames(lambda_mean)[i])
    for (m in seq_len(M)) row <- paste0(row, sprintf("  %8.3f", lambda_mean[i, m]))
    cat(row, "\n")
  }

  cat("\n")
  if (!is.null(x$nLc)) {
    cat("  Note: CSMF correction applied (Stan saw a subset of records).\n")
  }

  invisible(list(
    pi_summary   = data.frame(cause = names(pi_mean),
                              mean  = pi_mean,
                              lower = pi_lower,
                              upper = pi_upper,
                              row.names = NULL),
    lambda_table = as.data.frame(lambda_mean)
  ))
}


# ============================================================
# plot.BFL — CSMF dotplot with credible intervals
# ============================================================

#' @rdname BFL_methods
#' @param x A \code{"BFL"} object.
#' @param ci_prob Credible interval width. Default \code{0.95}.
#' @param title Optional plot title.
#' @param ... Ignored.
#' @export
plot.BFL <- function(x, ci_prob = 0.95, title = NULL, ...) {

  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required for plot.BFL().")

  alpha    <- (1 - ci_prob) / 2
  pi_mean  <- colMeans(x$pi)
  pi_lower <- apply(x$pi, 2, quantile, probs = alpha)
  pi_upper <- apply(x$pi, 2, quantile, probs = 1 - alpha)

  df <- data.frame(
    cause = x$causes,
    mean  = pi_mean,
    lower = pi_lower,
    upper = pi_upper
  )
  df <- df[order(df$mean, decreasing = TRUE), ]
  df$cause <- factor(df$cause, levels = rev(df$cause))

  plot_title <- if (!is.null(title)) title else "Cause-Specific Mortality Fractions"

  ggplot2::ggplot(df, ggplot2::aes(x = mean, y = cause)) +
    ggplot2::geom_segment(
      ggplot2::aes(x = lower, xend = upper, yend = cause),
      linewidth = 0.8, colour = "steelblue"
    ) +
    ggplot2::geom_point(size = 3, colour = "steelblue") +
    ggplot2::scale_x_continuous(
      labels = scales::percent_format(accuracy = 1),
      limits = c(0, NA)
    ) +
    ggplot2::labs(
      title = plot_title,
      x     = sprintf("Prevalence (%d%% credible interval)", round(ci_prob * 100)),
      y     = NULL
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())
}


# ============================================================
# print.BFL_score — quick three-number summary
# ============================================================

#' @rdname BFL_methods
#' @param x A \code{"BFL_score"} object.
#' @param ... Ignored.
#' @export
print.BFL_score <- function(x, ...) {
  cat("BFL Score\n")
  cat(strrep("\u2500", 40), "\n", sep = "")
  cat(sprintf("  Top-1 accuracy  : %.1f%%\n", x$top1_acc     * 100))
  cat(sprintf("  Balanced acc.   : %.1f%%\n", x$balanced_acc * 100))
  cat(sprintf("  CSMF accuracy   : %.1f%%\n", x$csmf_acc     * 100))
  cat(strrep("\u2500", 40), "\n", sep = "")
  cat("  Run summary() for per-cause detail, plot() for visuals.\n")
  invisible(x)
}


# ============================================================
# summary.BFL_score — full per-cause breakdown
# ============================================================

#' @rdname BFL_methods
#' @param object A \code{"BFL_score"} object.
#' @param ... Ignored.
#' @export
summary.BFL_score <- function(object, ...) {

  x      <- object
  causes <- names(x$pi_true)

  # Per-cause recall from confusion matrix
  conf   <- as.matrix(x$conf_mat)
  recall <- diag(conf) / pmax(rowSums(conf), 1L)
  recall <- recall[causes]

  # Inline bar helper (max width 10)
  .bar <- function(p, width = 10) {
    n <- max(0L, round(p * width))
    paste0(strrep("\u2588", n), strrep("\u2591", width - n))
  }

  cat("BFL Score\n")
  cat(strrep("\u2500", 65), "\n", sep = "")
  cat(sprintf("  Top-1 accuracy  : %.1f%%   Balanced: %.1f%%   CSMF: %.1f%%\n",
              x$top1_acc * 100, x$balanced_acc * 100, x$csmf_acc * 100))
  cat(strrep("\u2500", 65), "\n", sep = "")

  cat(sprintf("\nPer-Cause Detail:\n"))
  cat(sprintf("  %-22s %7s %7s %8s %7s\n",
              "Cause", "TruePrev", "EstPrev", "CSMFerr", "Recall"))
  cat(sprintf("  %s\n", strrep("-", 62)))

  # Sort by true prevalence descending
  ord <- order(x$pi_true, decreasing = TRUE)
  for (i in ord) {
    cat(sprintf("  %-22s %6.1f%% %6.1f%% %+7.1f%% %6.1f%%  %s\n",
                causes[i],
                x$pi_true[i]              * 100,
                x$pi_hat[i]               * 100,
                x$csmf_error_by_cause[i]  * 100,
                recall[i]                 * 100,
                .bar(recall[i])))
  }
  cat("\n")

  invisible(list(
    per_cause = data.frame(
      cause      = causes,
      pi_true    = as.numeric(x$pi_true),
      pi_hat     = as.numeric(x$pi_hat),
      csmf_error = as.numeric(x$csmf_error_by_cause),
      recall     = as.numeric(recall),
      row.names  = NULL
    ),
    overall = list(
      top1_acc     = x$top1_acc,
      balanced_acc = x$balanced_acc,
      csmf_acc     = x$csmf_acc
    )
  ))
}


# ============================================================
# plot.BFL_score — per-cause recall + CSMF error
# ============================================================

#' @rdname BFL_methods
#' @param x A \code{"BFL_score"} object.
#' @param title Optional title prefix for both panels.
#' @param ... Ignored.
#' @export
plot.BFL_score <- function(x, title = NULL, ...) {

  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required for plot.BFL_score().")

  gg           <- ggplot2::ggplot
  title_prefix <- if (is.null(title) || !nzchar(title)) "" else paste0(title, " \u2014 ")
  causes       <- names(x$pi_true)

  # Per-cause recall from confusion matrix diagonal
  conf   <- as.matrix(x$conf_mat)
  recall <- diag(conf) / pmax(rowSums(conf), 1L)

  # ------------------------------------------------------------------
  # 1) Per-cause recall bar chart (sorted descending)
  # ------------------------------------------------------------------
  df_rec <- data.frame(cause  = causes,
                       recall = as.numeric(recall[causes]))
  df_rec <- df_rec[order(df_rec$recall, decreasing = TRUE), , drop = FALSE]
  df_rec$cause <- factor(df_rec$cause, levels = rev(df_rec$cause))

  p_recall <- gg(df_rec, ggplot2::aes(x = recall, y = cause)) +
    ggplot2::geom_col(fill = "steelblue") +
    ggplot2::scale_x_continuous(
      labels = scales::percent_format(accuracy = 1),
      limits = c(0, 1)
    ) +
    ggplot2::labs(
      title = paste0(title_prefix, "Top-1 accuracy by cause (recall)"),
      x     = "Recall",
      y     = NULL
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())

  # ------------------------------------------------------------------
  # 2) Per-cause CSMF error bar chart (sorted by magnitude, two-toned)
  # ------------------------------------------------------------------
  err    <- x$csmf_error_by_cause
  df_err <- data.frame(
    cause = if (!is.null(names(err))) names(err) else paste0("Cause", seq_along(err)),
    error = as.numeric(err)
  )
  df_err <- df_err[order(abs(df_err$error), decreasing = TRUE), , drop = FALSE]
  df_err$cause <- factor(df_err$cause, levels = rev(df_err$cause))

  p_err <- gg(df_err, ggplot2::aes(x = error, y = cause, fill = error > 0)) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.4) +
    ggplot2::geom_col(show.legend = FALSE) +
    ggplot2::scale_fill_manual(
      values = c("TRUE" = "firebrick", "FALSE" = "steelblue")
    ) +
    ggplot2::scale_x_continuous(
      labels = scales::percent_format(accuracy = 0.1)
    ) +
    ggplot2::labs(
      title = paste0(title_prefix, "CSMF error by cause (\u03c0\u0302 \u2212 \u03c0)"),
      x     = NULL,
      y     = NULL
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())

  list(recall = p_recall, csmf_error = p_err)
}
