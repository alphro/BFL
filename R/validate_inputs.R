# ============================================================
# Input validation helpers (internal)
# ============================================================

# ------------------------------------------------------------------
# INTERNAL NOTE: validate_local_summaries()
#
# Ensures all local model summaries:
# 1) Were evaluated on the SAME target dataset (same N, same P).
# 2) Contain per-record class scores with consistent column labeling.
# 3) Share the same multiset of row hashes (order may differ; duplicates supported).
#
# Prevents silent failures like:
# - different target datasets mixed together
# - row misalignment across models
# - duplicate-row mismatches during permutation
# - posterior_phi columns not matching cause_ids
# ------------------------------------------------------------------

#' Validate structure of local summaries (internal)
#'
#' Internal helper used by \code{run_BFL()} to sanity-check local model outputs
#' before alignment and Stan aggregation.
#'
#' @param local_summaries List of local summary objects.
#' @return TRUE if validation passes; otherwise throws an error.
#' @noRd
validate_local_summaries <- function(local_summaries) {

  stopifnot(is.list(local_summaries))
  stopifnot(length(local_summaries) >= 1)

  for (x in local_summaries) {
    stopifnot(is.list(x))
    stopifnot(is.matrix(x$posterior_phi))
    stopifnot(!is.null(x$cause_ids))

    cid <- x$cause_ids
    if (is.factor(cid)) cid <- as.character(cid)
    stopifnot(is.character(cid))
    stopifnot(ncol(x$posterior_phi) == length(cid))

    stopifnot(!is.null(x$target_info))
    stopifnot(!is.null(x$target_info$row_hash))
    stopifnot(is.character(x$target_info$row_hash))

    # Strong invariants (helpful when debugging)
    if (!is.null(x$target_info$N)) stopifnot(x$target_info$N == nrow(x$posterior_phi))
    if (!is.null(x$target_info$P)) stopifnot(is.numeric(x$target_info$P) || is.integer(x$target_info$P))
    stopifnot(length(x$target_info$row_hash) == nrow(x$posterior_phi))
  }

  infos <- lapply(local_summaries, `[[`, "target_info")

  # Require N/P in target_info (your pipeline relies on these)
  Ns <- vapply(infos, function(z) z$N, integer(1))
  Ps <- vapply(infos, function(z) z$P, integer(1))
  stopifnot(length(unique(Ns)) == 1)
  stopifnot(length(unique(Ps)) == 1)

  row_hash_lens <- vapply(infos, function(z) length(z$row_hash), integer(1))
  stopifnot(length(unique(row_hash_lens)) == 1)
  stopifnot(unique(row_hash_lens) == unique(Ns))

  # Order-insensitive, duplicate-safe equality check
  ref <- sort(table(infos[[1]]$row_hash))
  ok_same_multiset <- vapply(
    infos,
    function(z) identical(sort(table(z$row_hash)), ref),
    logical(1)
  )
  stopifnot(all(ok_same_multiset))

  TRUE
}

# ------------------------------------------------------------------
# INTERNAL NOTE: validate_global_fit()
#
# Ensures run_BFL() output has the minimum structure needed for prediction/scoring.
# ------------------------------------------------------------------

#' Validate structure of a fitted global BFL object (internal)
#'
#' Internal helper used by downstream prediction/scoring functions.
#'
#' @param global_fit Object returned by \code{run_BFL()}.
#' @return TRUE if validation passes; otherwise throws an error.
#' @noRd
validate_global_fit <- function(global_fit) {
  stopifnot(is.list(global_fit))

  stopifnot(!is.null(global_fit$pi))
  stopifnot(!is.null(global_fit$lambda))
  stopifnot(!is.null(global_fit$phi))
  stopifnot(!is.null(global_fit$causes))
  stopifnot(!is.null(global_fit$model_names))
  stopifnot(!is.null(global_fit$row_hash))

  # Basic shape sanity (won’t over-assume exact dimensions)
  stopifnot(is.matrix(global_fit$pi) || is.array(global_fit$pi))
  stopifnot(is.array(global_fit$lambda))
  stopifnot(is.array(global_fit$phi))

  # phi should be N x C x M
  stopifnot(length(dim(global_fit$phi)) == 3)
  stopifnot(dim(global_fit$phi)[2] == length(global_fit$causes))
  stopifnot(dim(global_fit$phi)[3] == length(global_fit$model_names))
  stopifnot(length(global_fit$row_hash) == dim(global_fit$phi)[1])

  TRUE
}
