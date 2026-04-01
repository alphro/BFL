# ============================================================
# BFL2 — End-to-End Test Script
#
# Usage (from package root):
#   devtools::document()
#   devtools::load_all()
#   source("test_bfl.R")
#
# Tests all 4 BFL variants: Base, Domain, Partial, Mix
#
# Data has real signal: each cause has 5 marker symptoms
# with high probability (0.75) vs background (0.10).
# Target prevalence is skewed to test CSMF estimation.
# ============================================================

cat("╔══════════════════════════════════════════╗\n")
cat("║         BFL2 End-to-End Test             ║\n")
cat("╚══════════════════════════════════════════╝\n\n")


# ============================================================
# 0.  Synthetic data WITH signal
# ============================================================
cat("── Generating synthetic data (with signal) ─\n")

causes <- c("Malaria", "HIV", "TB", "Pneumonia")
C  <- length(causes)
P  <- 20   # 5 marker symptoms per cause, 0 overlap

# Symptom probability matrix: cause x symptom
# Cause k has high prob on symptoms ((k-1)*5+1) : (k*5)
SIG_HI <- 0.75
SIG_LO <- 0.10
sig_mat <- matrix(SIG_LO, nrow = C, ncol = P)
for (k in seq_len(C)) sig_mat[k, ((k-1)*5 + 1):(k*5)] <- SIG_HI

# Helper: draw n rows with cause-specific symptom profiles
make_data <- function(n, cause_probs = rep(1/C, C), seed = 1) {
  set.seed(seed)
  Y <- sample(causes, n, replace = TRUE, prob = cause_probs)
  X <- matrix(0L, nrow = n, ncol = P)
  for (i in seq_len(n)) {
    k      <- which(causes == Y[i])
    X[i, ] <- rbinom(P, 1, sig_mat[k, ])
  }
  list(X = X, Y = Y)
}

# Training sites: roughly uniform cause distribution
N_tr  <- 150
site1 <- make_data(N_tr, seed = 1)
site2 <- make_data(N_tr, seed = 2)
site3 <- make_data(N_tr, seed = 3)

# Target: skewed prevalence — harder test for CSMF estimation
N_tgt        <- 120
true_prev    <- c(Malaria = 0.40, HIV = 0.30, TB = 0.20, Pneumonia = 0.10)
target_data  <- make_data(N_tgt, cause_probs = true_prev, seed = 99)
X_target     <- target_data$X
Y_true       <- target_data$Y

cat(sprintf("  Training: 3 sites x %d rows\n", N_tr))
cat(sprintf("  Target  : %d rows | true prev: %s\n",
            N_tgt,
            paste(sprintf("%s %.0f%%", names(true_prev), true_prev*100),
                  collapse = " | ")))
cat(sprintf("  Signal  : %d marker symptoms per cause (p=%.2f vs background p=%.2f)\n\n",
            P / C, SIG_HI, SIG_LO))

# MCMC settings — small but functional
lcva_args       <- list(K = 3, Nitr = 150, thin = 1, seed = 1, pred_Nitr = 300)
lcva_args_small <- list(K = 2, Nitr = 100, thin = 1, seed = 1, pred_Nitr = 200)
stan_args       <- list(iter = 500, chains = 2, seed = 1)


# ============================================================
# 1.  compute_row_hashes
# ============================================================
cat("── compute_row_hashes ──────────────────────\n")
hashes <- compute_row_hashes(X_target)
stopifnot(length(hashes) == N_tgt, is.character(hashes), !anyNA(hashes))
cat(sprintf("  OK — %d hashes, all unique: %s\n\n",
            length(hashes), length(unique(hashes)) == N_tgt))


# ============================================================
# 2.  Base local_summaries (multi-site call)
# ============================================================
cat("── run_lcva  (3 training sites → X_target) ─\n")
local_summaries <- run_lcva(
  X_train   = list(site1 = site1$X, site2 = site2$X, site3 = site3$X),
  Y_train   = list(site1 = site1$Y, site2 = site2$Y, site3 = site3$Y),
  X_target  = X_target,
  lcva_args = lcva_args
)
stopifnot(is.list(local_summaries), length(local_summaries) == 3)
cat(sprintf("  OK — %d local summaries: %s\n\n",
            length(local_summaries),
            paste(names(local_summaries), collapse = ", ")))


# ============================================================
# 3.  Label masks
# ============================================================
cat("── make_label_mask (all 4 variants) ────────\n")
mask_base <- make_label_mask(Y_true, model = "Base", seed = 1)
mask_dom  <- make_label_mask(Y_true, model = "Domain",  miss_prop = 0.65, seed = 1)
mask_par  <- make_label_mask(Y_true, model = "Partial", miss_prop = 0.65, seed = 1)
mask_mix  <- make_label_mask(Y_true, model = "Mix",
                              miss_prop = 0.65,
                              mix_domain_prop_within_labeled = 0.50,
                              seed = 1)
cat(sprintf("  Base    stan=%d\n",   length(mask_base$stan_idx)))
cat(sprintf("  Domain  labeled=%d  unlabeled=%d\n",
            length(mask_dom$labeled_idx), length(mask_dom$missing_idx)))
cat(sprintf("  Partial labeled=%d  unlabeled=%d\n",
            length(mask_par$labeled_idx), length(mask_par$missing_idx)))
cat(sprintf("  Mix     domain=%d  partial=%d  unlabeled=%d\n\n",
            length(mask_mix$domain_idx),
            length(mask_mix$partial_idx),
            length(mask_mix$missing_idx)))


# ============================================================
# 4.  Domain / Mix local summaries
#
# Domain (Zoey's design):
#   - All sites predict on X_UNLABELED only
#   - Self-site trained on labeled rows, predicts on X_UNLABELED
#   - run_BFL gets X_full + Y_partial so it can compute nLc
#     Stan uses no-partial model; labeled rows are outside stan_idx
#
# Mix:
#   - All sites predict on X_FULL
#   - Self-site trained on domain rows, also predicts on X_FULL
#   - run_BFL gets X_full + Y_partial (partial-label Stan model)
# ============================================================
cat("── run_lcva  (domain / mix local summaries) ─\n")

X_unlabeled_dom <- X_target[mask_dom$missing_idx, ]
X_labeled_dom   <- X_target[mask_dom$labeled_idx, ]
Y_labeled_dom   <- Y_true[mask_dom$labeled_idx]

# Domain: all sites predict on unlabeled rows only
local_summaries_domain <- c(
  run_lcva(
    X_train   = list(site1 = site1$X, site2 = site2$X, site3 = site3$X),
    Y_train   = list(site1 = site1$Y, site2 = site2$Y, site3 = site3$Y),
    X_target  = X_unlabeled_dom,
    lcva_args = lcva_args
  ),
  list(self = run_lcva(X_labeled_dom, Y_labeled_dom,
                       X_unlabeled_dom, lcva_args_small))
)

# Y_partial for Domain: known labels, NA for unlabeled (stan) rows
Y_partial_dom                        <- rep(NA_character_, N_tgt)
Y_partial_dom[mask_dom$labeled_idx] <- Y_true[mask_dom$labeled_idx]

# Mix: all sites predict on full X_target
local_summaries_mix <- c(
  local_summaries,
  list(self = run_lcva(
    X_target[mask_mix$domain_idx, ],
    Y_true[mask_mix$domain_idx],
    X_target,
    lcva_args_small
  ))
)

cat(sprintf("  Domain summaries: %d sites | Mix summaries: %d sites\n\n",
            length(local_summaries_domain), length(local_summaries_mix)))


# ============================================================
# Sanity checks on a BFL_score
# These thresholds are intentionally conservative given small
# Nitr, but should reliably beat random chance (0.25) when
# real signal is present.
# ============================================================
THRESH_TOP1 <- 0.30   # random chance = 0.25 with 4 causes
THRESH_CSMF <- 0.50   # rough lower bound for reasonable CSMF
MAX_CSMF_ERR <- 0.20  # max tolerated |pi_hat - pi_true| per cause

check_score <- function(sc, label) {
  cat(sprintf("  [sanity checks — %s]\n", label))

  # Beat random chance
  if (sc$top1_acc < THRESH_TOP1)
    warning(sprintf("  !! top1_acc = %.1f%% — below threshold %.0f%%",
                    sc$top1_acc * 100, THRESH_TOP1 * 100))
  else
    cat(sprintf("    top1_acc    %.1f%% > %.0f%%  OK\n",
                sc$top1_acc * 100, THRESH_TOP1 * 100))

  # CSMF accuracy reasonable
  if (sc$csmf_acc < THRESH_CSMF)
    warning(sprintf("  !! csmf_acc = %.1f%% — below threshold %.0f%%",
                    sc$csmf_acc * 100, THRESH_CSMF * 100))
  else
    cat(sprintf("    csmf_acc    %.1f%% > %.0f%%  OK\n",
                sc$csmf_acc * 100, THRESH_CSMF * 100))

  # pi_hat sums to ~1
  pi_sum <- sum(sc$pi_hat)
  stopifnot(abs(pi_sum - 1.0) < 0.02)
  cat(sprintf("    sum(pi_hat) %.4f ≈ 1  OK\n", pi_sum))

  # Per-cause CSMF error within tolerance
  worst_cause <- names(which.max(abs(sc$csmf_error_by_cause)))
  worst_err   <- max(abs(sc$csmf_error_by_cause))
  if (worst_err > MAX_CSMF_ERR)
    warning(sprintf("  !! worst CSMF error = %.1f%% for %s (threshold %.0f%%)",
                    worst_err * 100, worst_cause, MAX_CSMF_ERR * 100))
  else
    cat(sprintf("    max CSMF err %.1f%% (%s) < %.0f%%  OK\n",
                worst_err * 100, worst_cause, MAX_CSMF_ERR * 100))

  cat("\n")
}


# ============================================================
# Helper: fit → print/summary/plot/score → sanity check
# ============================================================
results <- list()

test_variant <- function(label, fit, Y_eval, eval_idx = NULL) {

  cat(sprintf("┌─ %s ", label))
  cat(strrep("─", max(1, 42 - nchar(label))), "┐\n", sep = "")

  stopifnot(inherits(fit, "BFL"))

  # print / summary / plot
  cat("\n"); print(fit)
  sm <- summary(fit)
  stopifnot(!is.null(sm$pi_summary), !is.null(sm$lambda_table))
  p_fit <- plot(fit)
  stopifnot(inherits(p_fit, "gg"))
  cat("  plot.BFL OK\n")

  # site similarity
  cor_mat <- site_similarity_BFL(fit)
  stopifnot(is.matrix(cor_mat), nrow(cor_mat) == length(fit$model_names))
  cat(sprintf("  site_similarity_BFL: %dx%d  OK\n", nrow(cor_mat), ncol(cor_mat)))

  # score
  sc <- score_BFL(fit, Y_eval = Y_eval, eval_idx = eval_idx)
  stopifnot(inherits(sc, "BFL_score"))

  cat("\n"); print(sc)
  sm_sc <- summary(sc)
  stopifnot(!is.null(sm_sc$per_cause), !is.null(sm_sc$overall))
  plots <- plot(sc)
  stopifnot(is.list(plots), inherits(plots$recall, "gg"), inherits(plots$csmf_error, "gg"))
  cat("  plot.BFL_score OK\n\n")

  check_score(sc, label)

  cat(sprintf("└─ %s PASSED ", label))
  cat(strrep("─", max(1, 39 - nchar(label))), "┘\n\n", sep = "")

  results[[label]] <<- list(
    top1     = sc$top1_acc,
    balanced = sc$balanced_acc,
    csmf     = sc$csmf_acc,
    n_eval   = if (is.null(eval_idx)) N_tgt else length(eval_idx)
  )

  invisible(sc)
}


# ============================================================
# 5.  BASE
# ============================================================
cat("═══ Fitting BASE ════════════════════════════\n\n")
fit_base <- run_BFL(local_summaries, X_target, stan_args = stan_args)
test_variant("BASE", fit_base, Y_true)


# ============================================================
# 6.  DOMAIN
# ============================================================
cat("═══ Fitting DOMAIN ══════════════════════════\n\n")
# X_full so run_BFL knows n_total; Y_partial_dom so nLc is computed.
# Stan uses no-partial model because labeled rows are outside stan_idx.
fit_dom <- run_BFL(local_summaries_domain, X_target,
                   Y_target  = Y_partial_dom,
                   stan_args = stan_args)
test_variant("DOMAIN", fit_dom, Y_true, eval_idx = mask_dom$missing_idx)


# ============================================================
# 7.  PARTIAL
# ============================================================
cat("═══ Fitting PARTIAL ═════════════════════════\n\n")
fit_par <- run_BFL(local_summaries, X_target,
                   Y_target  = mask_par$Y_obs,
                   stan_args = stan_args)
test_variant("PARTIAL", fit_par, Y_true, eval_idx = mask_par$missing_idx)


# ============================================================
# 8.  MIX
# ============================================================
cat("═══ Fitting MIX ═════════════════════════════\n\n")
fit_mix <- run_BFL(local_summaries_mix, X_target,
                   Y_target  = mask_mix$Y_obs,
                   stan_args = stan_args)
test_variant("MIX", fit_mix, Y_true, eval_idx = mask_mix$missing_idx)


# ============================================================
# 9.  Comparison table across variants
# ============================================================
cat("╔══════════════════════════════════════════════════╗\n")
cat("║              Variant Comparison                  ║\n")
cat("╠══════════════════════════════════════════════════╣\n")
cat(sprintf("║  %-10s  %8s  %10s  %8s  %6s  ║\n",
            "Variant", "Top-1", "Balanced", "CSMF", "N_eval"))
cat("╠══════════════════════════════════════════════════╣\n")
for (nm in names(results)) {
  r <- results[[nm]]
  cat(sprintf("║  %-10s  %7.1f%%  %9.1f%%  %7.1f%%  %6d  ║\n",
              nm, r$top1*100, r$balanced*100, r$csmf*100, r$n_eval))
}
cat("╠══════════════════════════════════════════════════╣\n")
cat(sprintf("║  True prevalence: %s  ║\n",
            paste(sprintf("%s=%.0f%%", names(true_prev), true_prev*100),
                  collapse = " ")))
cat("╠══════════════════════════════════════════════════╣\n")
cat("║  Domain/Partial/Mix eval on unlabeled rows only  ║\n")
cat("╚══════════════════════════════════════════════════╝\n\n")

cat("All variants passed.\n")
