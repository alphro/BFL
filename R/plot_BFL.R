# ============================================================
# plot_BFL.R
# ============================================================
# PURPOSE:
#   Plot diagnostics from score_BFL(): overall metrics,
#   confusion matrix, and per-cause CSMF errors.
# ============================================================

#' Plot BFL scoring outputs
#'
#' Creates standard plots from a \code{score_BFL()} result:
#' (1) overall Top-1 accuracy + CSMF accuracy,
#' (2) confusion matrix heatmap (Top-1 predictions),
#' (3) per-cause CSMF error (\code{pi_hat - pi_true}).
#'
#' @param score Output list from \code{score_BFL()}.
#' @param normalize_conf Logical; if TRUE, confusion matrix is row-normalized
#'   (each true cause sums to 1). Default FALSE (raw counts).
#' @param conf_max_labels Integer; if number of causes exceeds this, the confusion
#'   matrix plot is omitted (still returned as NULL) to avoid unreadable plots.
#' @param title Optional title prefix used in plot titles.
#'
#' @return A named list of ggplot objects:
#' \describe{
#'   \item{metrics}{Bar chart for overall metrics.}
#'   \item{confusion}{Confusion matrix heatmap (or NULL if too many causes).}
#'   \item{csmf_error}{Per-cause CSMF error bar chart.}
#' }
#'
#' @examples
#' \dontrun{
#' sc <- score_BFL(pred, mask = make_label_mask(Y, model = "Partial"))
#' plots <- plot_BFL(sc, normalize_conf = TRUE, title = "Mexico / seed 1001")
#' print(plots$metrics)
#' print(plots$confusion)
#' print(plots$csmf_error)
#' }
#'
#' @export
plot_BFL <- function(score,
                     normalize_conf = FALSE,
                     conf_max_labels = 40,
                     title = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plot_BFL().")
  }
  if (!requireNamespace("reshape2", quietly = TRUE)) {
    stop("Package 'reshape2' is required for plot_BFL() (for melting matrices).")
  }
  
  stopifnot(is.list(score))
  stopifnot(!is.null(score$top1_acc), !is.null(score$csmf_acc))
  stopifnot(!is.null(score$conf_mat))
  stopifnot(!is.null(score$csmf_error_by_cause))
  
  gg <- ggplot2::ggplot
  
  title_prefix <- if (is.null(title) || !nzchar(title)) "" else paste0(title, " — ")
  
  # ------------------------------------------------------------
  # 1) Overall metrics
  # ------------------------------------------------------------
  df_metrics <- data.frame(
    metric = c("Top-1 accuracy", "CSMF accuracy"),
    value  = c(as.numeric(score$top1_acc), as.numeric(score$csmf_acc))
  )
  
  p_metrics <- gg(df_metrics, ggplot2::aes(x = metric, y = value)) +
    ggplot2::geom_col() +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(
      title = paste0(title_prefix, "Overall metrics"),
      x = NULL, y = NULL
    ) +
    ggplot2::theme_minimal()
  
  # ------------------------------------------------------------
  # 2) Confusion matrix heatmap
  # ------------------------------------------------------------
  conf <- score$conf_mat
  stopifnot(is.matrix(conf))
  
  causes <- colnames(conf)
  if (is.null(causes)) causes <- rownames(conf)
  
  p_conf <- NULL
  if (!is.null(causes) && length(causes) <= conf_max_labels) {
    conf_plot <- conf
    
    if (isTRUE(normalize_conf)) {
      rs <- rowSums(conf_plot)
      rs[rs == 0] <- 1
      conf_plot <- conf_plot / rs
    }
    
    conf_long <- reshape2::melt(conf_plot)
    colnames(conf_long) <- c("true", "pred", "value")
    
    p_conf <- gg(conf_long, ggplot2::aes(x = pred, y = true, fill = value)) +
      ggplot2::geom_tile() +
      ggplot2::labs(
        title = paste0(
          title_prefix,
          "Confusion matrix",
          if (normalize_conf) " (row-normalized)" else " (counts)"
        ),
        x = "Predicted cause", y = "True cause", fill = NULL
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        panel.grid = ggplot2::element_blank()
      )
  }
  
  # ------------------------------------------------------------
  # 3) Per-cause CSMF error (pi_hat - pi_true)
  # ------------------------------------------------------------
  err <- score$csmf_error_by_cause
  stopifnot(is.numeric(err))
  
  df_err <- data.frame(
    cause = names(err) %||% paste0("Cause", seq_along(err)),
    error = as.numeric(err)
  )
  
  # show biggest deviations first
  df_err <- df_err[order(abs(df_err$error), decreasing = TRUE), , drop = FALSE]
  df_err$cause <- factor(df_err$cause, levels = df_err$cause)
  
  p_err <- gg(df_err, ggplot2::aes(x = cause, y = error)) +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::geom_col() +
    ggplot2::labs(
      title = paste0(title_prefix, "CSMF error by cause (pi_hat - pi_true)"),
      x = NULL, y = NULL
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  
  list(
    metrics   = p_metrics,
    confusion = p_conf,
    csmf_error = p_err
  )
}

# small infix helper (kept internal to this file)
`%||%` <- function(x, y) if (!is.null(x)) x else y