predict_BFL2 <- function(global_fit, seed = NULL) {

  if (!is.null(seed)) {
    set.seed(seed)
  }

  stopifnot(
    !is.null(global_fit$phi),
    !is.null(global_fit$lambda),
    !is.null(global_fit$pi),
    !is.null(global_fit$causes)
  )

  phi     <- global_fit$phi
  lambda  <- global_fit$lambda
  pi_draw <- global_fit$pi
  causes  <- as.character(global_fit$causes)

  S <- dim(lambda)[1]
  I <- dim(phi)[1]
  C <- dim(phi)[2]
  M <- dim(phi)[3]

  draws_int <- matrix(NA_integer_, nrow = S, ncol = I)

  for (s in seq_len(S)) {

    for (i in seq_len(I)) {

      score <- numeric(C)

      for (c in seq_len(C)) {

        score[c] <-
          pi_draw[s, c] *
          sum(phi[i, c, ] * lambda[s, c, ])
      }

      if (!is.finite(sum(score)) || sum(score) <= 0) {
        score[] <- 1 / C
      } else {
        score <- score / sum(score)
      }

      # Zoey-style stochastic categorical sampling
      draws_int[s, i] <- sample.int(
        C,
        size = 1,
        prob = score
      )
    }
  }

  draws_lab <- matrix(
    causes[draws_int],
    nrow = nrow(draws_int),
    ncol = ncol(draws_int)
  )

  out <- list(
    draws_int = draws_int,
    draws     = draws_lab,
    causes    = causes,
    pi        = pi_draw,
    model     = global_fit$model
  )

  if (!is.null(global_fit$row_hash)) {
    out$row_hash <- global_fit$row_hash
  }

  out
}
