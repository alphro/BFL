fit_lcva <- function(X, Y = NULL, K = 5, lcva_args = list()) {
  if (!requireNamespace("LCVA", quietly = TRUE)) {
    stop("Package 'LCVA' is required. Please install it.")
  }

  X <- as.matrix(X)
  storage.mode(X) <- "numeric"

  if (is.null(Y)) {
    stop("For now, fit_lcva requires Y (supervised LCVA).")
  }

  # Encode Y as consecutive integers 1..C
  Y <- as.integer(factor(Y))

  # defaults (override via lcva_args)
  Nitr <- if (!is.null(lcva_args$Nitr)) lcva_args$Nitr else 2000
  thin <- if (!is.null(lcva_args$thin)) lcva_args$thin else 2
  seed <- if (!is.null(lcva_args$seed)) lcva_args$seed else 12345
  model <- if (!is.null(lcva_args$model)) lcva_args$model else "S"
  verbose <- if (!is.null(lcva_args$verbose)) lcva_args$verbose else FALSE

  LCVA::LCVA.train(
    X = X,
    Y = Y,
    Domain = rep(1, length(Y)),
    K = K,
    model = model,
    Nitr = Nitr,
    thin = thin,
    seed = seed,
    verbose = verbose
  )
}
