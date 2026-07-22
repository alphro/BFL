# ==============================================================
# make_bfl_cities.R --- build inst/extdata/bfl_cities.rds
#
# The vignette starts from *pre-built local summaries* (the federated premise:
# each source shares only its trained summary, never its records). This script
# produces that bundle from the raw synthetic data (inst/extdata/va_cities.csv):
#
#   * fit a simple naive-Bayes model per SOURCE city   -> P(symptom | cause) = theta
#   * score the TARGET (Los Angeles) records with it   -> posterior_phi (N x C)
#   * attach cause_ids + compute_row_hashes()          -> a valid local_summary
#
# Must be run in R (compute_row_hashes() uses rlang::hash(), which run_BFL()
# recomputes internally, so the hashes have to be made the same way):
#   Rscript data-raw/make_bfl_cities.R      # from the package root, with BFL installed
# ==============================================================
library(BFL)   # compute_row_hashes()

va <- read.csv("data-raw/va_cities.csv", stringsAsFactors = FALSE)

symptoms <- setdiff(names(va), c("site", "cause"))
causes   <- sort(unique(va$cause))
sources  <- c("NewYork", "Seattle", "Austin")
target   <- "LosAngeles"

la       <- va[va$site == target, ]
X_target <- as.matrix(la[, symptoms]); storage.mode(X_target) <- "numeric"
Y_true   <- la$cause

# ---- fit naive-Bayes theta per source: P(symptom = 1 | cause) ----
alpha <- 1  # Laplace smoothing
fit_theta <- function(dat) {
  th <- matrix(0, length(causes), length(symptoms),
               dimnames = list(causes, symptoms))
  for (cc in causes) {
    sub <- dat[dat$cause == cc, symptoms, drop = FALSE]
    th[cc, ] <- (colSums(sub) + alpha) / (nrow(sub) + 2 * alpha)
  }
  th
}

# ---- score the target with a source's theta -> row-normalized posterior_phi ----
score_target <- function(theta, X) {
  logth  <- log(theta); log1th <- log(1 - theta)
  lp <- X %*% t(logth) + (1 - X) %*% t(log1th)   # N x C log-likelihood (uniform prior)
  lp <- lp - apply(lp, 1, max)                   # stabilize before exp
  p  <- exp(lp)
  p / rowSums(p)                                 # N x C, columns = causes
}

theta    <- lapply(setNames(sources, sources), function(s) fit_theta(va[va$site == s, ]))
row_hash <- compute_row_hashes(X_target)

local_summaries <- lapply(setNames(sources, sources), function(s) {
  phi <- score_target(theta[[s]], X_target)      # N x C, cols in `causes` order
  colnames(phi) <- causes
  list(
    posterior_phi = phi,
    cause_ids     = causes,
    target_info   = list(row_hash = row_hash,
                         N = nrow(phi),          # required by validate_local_summaries()
                         P = ncol(X_target))
  )
})

bfl_cities <- list(
  local_summaries = local_summaries,  # what run_BFL() consumes
  X_target        = X_target,         # the target's own records (needed to predict them)
  Y_true          = Y_true,           # LA truth, for masking + scoring in the vignette
  theta           = theta,            # per-source P(symptom | cause), for the intro
  causes          = causes,
  symptoms        = symptoms
)

saveRDS(bfl_cities, "inst/extdata/bfl_cities.rds")
cat("wrote inst/extdata/bfl_cities.rds\n")
str(bfl_cities, max.level = 2)
