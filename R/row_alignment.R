#' Row-hash alignment helpers (internal)
#'
#' Utilities for aligning per-record outputs from local models to a
#' canonical target row order using row-wise hashes.
#'
#' @details
#' These functions assume full row coverage: \code{h_src} must be a
#' permutation of \code{h_ref}, including duplicate counts.
#'
#' They prevent silent row misalignment before BFL aggregation.
#'
#' @name row_alignment_helpers
#' @keywords internal
NULL


# ------------------------------------------------------------
# Compute permutation from source hash to reference hash
# ------------------------------------------------------------
#' Compute row permutation from source to reference hash order
#'
#' Returns an integer permutation vector \code{perm} such that
#' \code{h_src[perm] == h_ref}. Consumes each source position at most once,
#' so duplicate hashes are handled correctly.
#'
#' @param h_src Character vector of row hashes for the source matrix.
#' @param h_ref Character vector of row hashes giving the desired order.
#'   Must be a permutation of \code{h_src} (same multiset).
#'
#' @return Integer vector of length \code{length(h_ref)} such that
#'   \code{h_src[perm] == h_ref}.
#'
#' @keywords internal
row_perm_from_hash <- function(h_src, h_ref) {

  stopifnot(length(h_src) == length(h_ref))
  stopifnot(identical(sort(h_src), sort(h_ref)))

  perm <- integer(length(h_ref))
  used <- logical(length(h_src))

  for (i in seq_along(h_ref)) {

    candidates <- which(h_src == h_ref[i] & !used)

    if (length(candidates) == 0)
      stop("No unused match found for hash.")

    perm[i] <- candidates[1]
    used[candidates[1]] <- TRUE
  }

  stopifnot(identical(h_src[perm], h_ref))

  perm
}


# ------------------------------------------------------------
# Align all local summaries to reference row order
# ------------------------------------------------------------
#' Align all local summaries to a reference row order
#'
#' Re-orders the \code{posterior_phi} rows of every local summary so they
#' match the canonical target row order given by \code{ref_row_hash}. Also
#' updates \code{target_info$row_hash} and \code{target_info$N} in each
#' element.
#'
#' @param local_summaries Named list of local model summary objects. Each
#'   element must contain \code{posterior_phi} (N x C matrix) and
#'   \code{target_info$row_hash} (character vector of length N).
#' @param ref_row_hash Character vector giving the canonical row hash order
#'   (produced by \code{compute_row_hashes()}). Must be a permutation of
#'   each site's \code{target_info$row_hash}.
#'
#' @return The same list structure as \code{local_summaries} with rows of
#'   \code{posterior_phi} permuted to match \code{ref_row_hash}.
#'
#' @keywords internal
align_rows_local_summaries <- function(local_summaries, ref_row_hash) {

  stopifnot(is.character(ref_row_hash))
  stopifnot(length(ref_row_hash) >= 1)

  lapply(local_summaries, function(x) {

    stopifnot(!is.null(x$target_info))
    stopifnot(!is.null(x$target_info$row_hash))

    perm <- row_perm_from_hash(
      h_src = x$target_info$row_hash,
      h_ref = ref_row_hash
    )

    x$posterior_phi <- x$posterior_phi[perm, , drop = FALSE]
    x$target_info$row_hash <- ref_row_hash
    x$target_info$N <- length(ref_row_hash)

    x
  })
}
