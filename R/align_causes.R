# PURPOSE:
# Aligns cause indices across local summaries so all sites
# refer to the same global cause set.

align_local_summaries <- function(local_summaries, ref_row_hash = NULL) {

  # Optional row alignment first (future-proof hook)
  if (!is.null(ref_row_hash)) {
    local_summaries <- align_rows_local_summaries(local_summaries, ref_row_hash)
  }

  global_causes <- sort(unique(unlist(
    lapply(local_summaries, function(x) x$cause_ids)
  )))

  aligned_phi <- lapply(local_summaries, function(x) {
    P <- nrow(x$posterior_phi)
    Cg <- length(global_causes)

    phi_aligned <- matrix(
      0,
      nrow = P,
      ncol = Cg,
      dimnames = list(NULL, global_causes)
    )

    idx <- match(x$cause_ids, global_causes)
    if (anyNA(idx)) {
      stop("Some local cause_ids are not in global_causes. Check cause_ids consistency.")
    }

    phi_aligned[, idx] <- x$posterior_phi
    phi_aligned
  })

  list(
    global_causes = global_causes,
    aligned_phi = aligned_phi,
    ref_row_hash = ref_row_hash
  )
}
