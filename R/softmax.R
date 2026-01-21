softmax <- function(x) {
  exp_x <- exp(x - max(x))
  exp_x / sum(exp_x)
}
