#' Run global BFL aggregation
#'
#' Aggregates per-model per-record class scores into a global Bayesian ensemble
#' using Stan. Each local model contributes an \eqn{N \times C_m} matrix over the
#' same target records; BFL learns posterior class prevalence (\code{pi}) and
#' per-class model weights (\code{lambda}).
#'
#' @param local_summaries Named list of local summary objects. Each element must contain:
#'   \itemize{
#'     \item \code{posterior_phi}: numeric matrix (\eqn{N \times C_m}) of per-record class
#'       scores/probabilities for the target set (rows align to \code{X_target}).
#'     \item \code{cause_ids}: character vector (length \eqn{C_m}) giving column labels.
#'     \item \code{target_info$row_hash}: character vector (length \eqn{N}) for row alignment.
#'   }
#' @param X_target Target feature matrix/data.frame. Used only to define canonical
#'   row identity via hashing.
#' @param Y_target Optional target labels (use NA for unknown). If NULL, runs the
#'   no-partial-label model.
#' @param label_shift Logical; if TRUE and \code{Y_target} provided, uses the
#'   unbalanced/label-shift variant.
#' @param model One of \code{"Base"}, \code{"Domain"}, \code{"Partial"}, \code{"Mix"}.
#'   Currently used for bookkeeping (stored in the output).
#' @param stan_args List of arguments passed to Stan sampling.
#'
#' @return A list with components:
#' \describe{
#'   \item{pi}{Posterior draws of class prevalence (S x C).}
#'   \item{lambda}{Posterior draws of per-class model weights (S x C x M).}
#'   \item{phi}{Per-model scores aligned to target rows and global classes (N x C x M).}
#'   \item{causes}{Character vector of global class labels (length C).}
#'   \item{model_names}{Character vector of model names (length M).}
#'   \item{row_hash}{Row hashes for \code{X_target} (length N).}
#'   \item{model}{User-specified BFL variant (bookkeeping).}
#' }
#'
#' @export
run_BFL <- function(
    local_summaries,
    X_target,
    Y_target = NULL,
    label_shift = FALSE,
    model = c("Base","Domain","Partial","Mix"),
    stan_args = list(iter = 2000, chains = 4, seed = 12345)
) {
  model <- match.arg(model)

  validate_local_summaries(local_summaries)

  X_target <- as.matrix(X_target)
  ref_row_hash <- hash_rows(X_target)

  aligned <- align_local_summaries(local_summaries, ref_row_hash = ref_row_hash)

  variant <- if (is.null(Y_target)) {
    "no_partial"
  } else if (isTRUE(label_shift)) {
    "unbalanced"
  } else {
    "balanced"
  }

  stan_data <- if (is.null(Y_target)) {

    build_bfl_stan_data(aligned)

  } else {

    y <- as.character(Y_target)
    stopifnot(length(y) == nrow(X_target))

    Y_known <- as.integer(!is.na(y))

    cause_to_idx <- setNames(seq_along(aligned$global_causes), aligned$global_causes)
    Y_idx <- unname(cause_to_idx[y])
    Y_idx[Y_known == 0] <- 1L

    if (any(Y_known == 1 & is.na(Y_idx))) {
      bad <- unique(y[Y_known == 1 & is.na(Y_idx)])
      stop("Y_target contains labels not in aligned$global_causes: ",
           paste(bad, collapse = ", "))
    }

    if (isTRUE(label_shift)) {
      build_bfl_stan_data_unbalanced(
        aligned,
        Y_known = Y_known,
        Y_idx   = as.integer(Y_idx)
      )
    } else {
      build_bfl_stan_data_balanced(
        aligned,
        Y_known = Y_known,
        Y_idx   = as.integer(Y_idx)
      )
    }
  }

  fit <- run_bfl_stan(stan_data, stan_args, variant = variant)

  phi_list <- lapply(aligned$aligned_phi, as.matrix)
  model_names <- names(phi_list)

  N <- nrow(phi_list[[1]])
  C <- ncol(phi_list[[1]])
  M <- length(phi_list)

  phi_arr <- array(0, dim = c(N, C, M))
  for (m in seq_len(M)) phi_arr[, , m] <- phi_list[[m]]

  list(
    pi = fit$pi,
    lambda = fit$lambda,
    phi = phi_arr,
    causes = aligned$global_causes,
    model_names = model_names,
    row_hash = ref_row_hash,
    model = model
  )
}
