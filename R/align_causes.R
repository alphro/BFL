# PURPOSE:
# Aligns cause indices across local summaries so all sites
# refer to the same global cause set.

#' Align local summaries to a global cause set (internal)
#'
#' @param local_summaries list of local summary objects
#' @param ref_row_hash optional reference hash for row alignment
#' @return list(global_causes, aligned_phi, ref_row_hash)
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
