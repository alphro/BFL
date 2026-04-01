# PURPOSE:
# Aligns cause indices across local summaries so all sites
# refer to the same global cause set.

#' Align local summaries to a global cause set (internal)
#'
#' Merges the per-site cause vocabularies into a single sorted global cause set
#' and pads each site's \code{posterior_phi} matrix to that full width, filling
#' unobserved causes with zero. Optionally re-orders rows to match a reference
#' row hash before cause alignment.
#'
#' @param local_summaries Named list of local model summary objects. Each element
#'   must contain \code{posterior_phi} (N x C matrix), \code{cause_ids} (character
#'   vector of length C), and \code{target_info$row_hash} (character vector of
#'   length N).
#' @param ref_row_hash Optional character vector of length N giving the canonical
#'   row hash order (produced by \code{compute_row_hashes()}). When supplied,
#'   each site's \code{posterior_phi} rows are permuted to match this order
#'   before cause alignment.
#'
#' @return A list with three elements:
#'   \describe{
#'     \item{global_causes}{Sorted character vector of all cause labels across
#'       all sites.}
#'     \item{aligned_phi}{Named list of N x C_global matrices with columns
#'       ordered by \code{global_causes}; columns absent from a site are zero.}
#'     \item{ref_row_hash}{The \code{ref_row_hash} argument (passed through).}
#'   }
#' @keywords internal
align_local_summaries <- function(local_summaries, ref_row_hash = NULL) {

  if (!is.null(ref_row_hash)) {
    local_summaries <- align_rows_local_summaries(local_summaries, ref_row_hash)
  }

  # Coerce cause_ids safely
  local_summaries <- lapply(local_summaries, function(x) {
    if (is.factor(x$cause_ids)) x$cause_ids <- as.character(x$cause_ids)
    x
  })

  global_causes <- sort(unique(unlist(
    lapply(local_summaries, function(x) x$cause_ids)
  )))

  aligned_phi <- lapply(local_summaries, function(x) {
    N  <- nrow(x$posterior_phi)
    Cg <- length(global_causes)

    phi_aligned <- matrix(
      0,
      nrow = N,
      ncol = Cg,
      dimnames = list(NULL, global_causes)
    )

    idx <- match(x$cause_ids, global_causes)
    if (anyNA(idx)) stop("Some local cause_ids are not in global_causes.")

    phi_aligned[, idx] <- x$posterior_phi
    phi_aligned
  })

  list(
    global_causes = global_causes,
    aligned_phi = aligned_phi,
    ref_row_hash = ref_row_hash
  )
}
