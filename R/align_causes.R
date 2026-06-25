# PURPOSE:
# Aligns cause indices across local summaries so all sites
# refer to the same global cause set.

#' Align local summaries to a global cause set (internal)
#'
#' Merges the per-site cause vocabularies into a single global cause set and
#' pads each site's \code{posterior_phi} matrix to that full width, filling
#' unobserved causes with zero. Optionally re-orders rows to match a reference
#' row hash before cause alignment.
#'
#' Cause order is determined by sorting the union of all sites' cause IDs.
#' When all IDs coerce to numeric without \code{NA}, numeric order is used
#' (e.g., PHMRC integer causes 1..34); otherwise lexicographic order is used.
#' Sorting is required for correct prediction decoding: the Gibbs sampler and
#' posterior-predictive sampler communicate via 1-indexed \emph{positions} into
#' \code{global_causes}, so position \eqn{k} must map to the same cause label
#' regardless of site-combination order.
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
#' @return A list with four elements:
#'   \describe{
#'     \item{global_causes}{Character vector of all cause labels across all
#'       sites, in sorted order (numeric if all IDs are numeric-like,
#'       lexicographic otherwise).}
#'     \item{aligned_phi}{Named list of N x C_global matrices with columns
#'       ordered by \code{global_causes}; columns absent from a site are zero.}
#'     \item{model_presence}{Integer matrix (C_global x M) where entry [c, m]
#'       is 1 if site m structurally includes cause c in its cause list, 0
#'       otherwise. This reflects cause-list membership, NOT phi values — a
#'       cause with all-zero phi due to underflow still gets model_presence=1
#'       so the issue remains visible to downstream diagnostics.}
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

  # Collect all unique cause IDs across sites, then sort into a stable order.
  #
  # Sorting is REQUIRED for correct prediction decoding. The Gibbs sampler
  # and posterior-predictive sampler communicate via 1-indexed POSITIONS into
  # this vector (draws_int), so position k must consistently map to the k-th
  # cause across all runs regardless of site-combination order.
  #
  # Sort rule:
  #   - If every cause ID coerces to a finite numeric without NA, sort
  #     numerically (handles non-contiguous integer IDs, e.g. PHMRC 1..34).
  #   - Otherwise sort lexicographically (handles string cause IDs, e.g.
  #     CHAMPS text labels). Document clearly if mixed/non-standard IDs arise.
  global_causes_raw <- unique(unlist(
    lapply(local_summaries, function(x) x$cause_ids)
  ))
  gc_num <- suppressWarnings(as.numeric(global_causes_raw))
  if (all(!is.na(gc_num))) {
    global_causes <- as.character(sort(gc_num))
  } else {
    global_causes <- sort(global_causes_raw)
  }

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

  # ------------------------------------------------------------------
  # Structural model_presence: C_global x M
  #
  # Entry [c, m] = 1 iff site m's cause_ids includes global_causes[c].
  # This reflects cause-list MEMBERSHIP, not phi values.  A source that
  # structurally contains a cause but whose phi has underflowed to all
  # zeros still gets model_presence=1 — the underflow is then visible as
  # a true-underflow warning in run_bfl_gibbs(), rather than being
  # silently reclassified as structural absence.
  # ------------------------------------------------------------------
  Cg <- length(global_causes)
  M  <- length(local_summaries)
  model_presence <- matrix(
    0L,
    nrow = Cg,
    ncol = M,
    dimnames = list(global_causes, names(local_summaries))
  )
  for (m in seq_len(M)) {
    idx <- match(as.character(local_summaries[[m]]$cause_ids), global_causes)
    if (anyNA(idx)) stop("Some local cause_ids are not in global_causes.")
    model_presence[idx, m] <- 1L
  }

  list(
    global_causes  = global_causes,
    aligned_phi    = aligned_phi,
    model_presence = model_presence,
    ref_row_hash   = ref_row_hash
  )
}
