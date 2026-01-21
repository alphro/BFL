safe_log <- function(x, eps = 1e-12) {
  log(pmax(x, eps))
}
