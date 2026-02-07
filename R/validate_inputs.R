validate_local_summaries <- function(local_summaries) {

  # ------------------------------------------------------------------
  # Basic container checks
  # ------------------------------------------------------------------
  stopifnot(is.list(local_summaries))
  stopifnot(length(local_summaries) >= 1)

  # ------------------------------------------------------------------
  # Validate structure of each local summary
  # ------------------------------------------------------------------
  for (x in local_summaries) {
    stopifnot(is.matrix(x$posterior_phi))
    stopifnot(!is.null(x$cause_ids))
    if (is.factor(x$cause_ids)) x$cause_ids <- as.character(x$cause_ids)
    stopifnot(is.character(x$cause_ids))
    stopifnot(ncol(x$posterior_phi) == length(x$cause_ids))
  }

  # ------------------------------------------------------------------
  # Ensure all local summaries were evaluated on the same target dataset
  # BUT allow different row orders:
  # - require a row-wise hash vector (length N)
  # - check that all summaries share the same multiset of row hashes
  # ------------------------------------------------------------------
  stopifnot(all(vapply(
    local_summaries,
    function(x) !is.null(x$target_info),
    logical(1)
  )))

  infos <- lapply(local_summaries, `[[`, "target_info")

  # dimensions must match
  Ns <- vapply(infos, function(z) z$N, integer(1))
  Ps <- vapply(infos, function(z) z$P, integer(1))
  stopifnot(length(unique(Ns)) == 1)
  stopifnot(length(unique(Ps)) == 1)

  # row_hash must exist and be length N
  has_row_hash <- vapply(infos, function(z) !is.null(z$row_hash), logical(1))
  stopifnot(all(has_row_hash))

  row_hash_lens <- vapply(infos, function(z) length(z$row_hash), integer(1))
  stopifnot(length(unique(row_hash_lens)) == 1)
  stopifnot(unique(row_hash_lens) == unique(Ns))

  # check same multiset of hashes across summaries (order-insensitive)
  ref <- sort(infos[[1]]$row_hash)
  ok_same_set <- vapply(infos, function(z) identical(sort(z$row_hash), ref), logical(1))
  stopifnot(all(ok_same_set))

  TRUE
}

validate_global_fit <- function(global_fit) {

  # ------------------------------------------------------------------
  # Minimal structural validation of a fitted BFL object
  # Ensures required components exist for downstream usage
  # ------------------------------------------------------------------
  stopifnot(is.list(global_fit))
  stopifnot(!is.null(global_fit$pi))
  stopifnot(!is.null(global_fit$phi))

  TRUE
}
