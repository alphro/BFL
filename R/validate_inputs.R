validate_local_summaries <- function(local_summaries) {
  stopifnot(is.list(local_summaries))
  stopifnot(length(local_summaries) >= 1)

  Ns <- vapply(local_summaries, function(x) nrow(x$posterior_phi), integer(1))
  stopifnot(length(unique(Ns)) == 1)  # <-- all must match target N

  for (x in local_summaries) {
    stopifnot(is.matrix(x$posterior_phi))
    stopifnot(ncol(x$posterior_phi) == length(x$cause_ids))
    if (is.factor(x$cause_ids)) x$cause_ids <- as.character(x$cause_ids)
    stopifnot(is.character(x$cause_ids))
  }

  TRUE
}

validate_global_fit <- function(global_fit) {
  stopifnot(is.list(global_fit))
  stopifnot(!is.null(global_fit$pi))
  stopifnot(!is.null(global_fit$phi))
}
