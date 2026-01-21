# PURPOSE:
# Aligns cause indices across local summaries so all sites
# refer to the same global cause set.

align_local_summaries <- function(local_summaries) {

  global_causes <- sort(unique(unlist(
    lapply(local_summaries, function(x) x$cause_ids)
  )))

  aligned_phi <- lapply(local_summaries, function(x) {
    P <- nrow(x$posterior_phi)
    Cg <- length(global_causes)

    phi_aligned <- matrix(0, nrow = P, ncol = Cg)

    idx <- match(x$cause_ids, global_causes)
    phi_aligned[, idx] <- x$posterior_phi

    phi_aligned
  })

  list(
    global_causes = global_causes,
    aligned_phi = aligned_phi
  )
}
