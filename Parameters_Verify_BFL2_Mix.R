############################################################
# BFL2 Parameters Verify — MIX (BALANCED)
#
# Compares BFL2::run_BFL() vs Zoey's execute_balance_mix().
# Standalone: run rm(list=ls()) + gc() then source this file.
#
# Note: Mix requires the domain/partial split of the target
# site AND a separate LCVA fit on domain rows only (the
# self-site model).  Both are computed here from scratch.
############################################################

rm(list = ls()); gc()

# ============================================================
# 0. CONFIG — edit paths before running
# ============================================================
ZOEY_REPO  <- "/Users/toastymac/Desktop/BFL/Zoey's Github/BFL_VA-main"
REPORT_DIR <- "/Users/toastymac/Desktop/BFL/BFL Reports/Report1/BFL Parameters Verify BFL2"

sites     <- c("Mexico", "AP", "Bohol", "Dar", "Pemba", "UP")
test_site <- "AP"
miss_prop <- 0.80
K         <- 5L
STAN_ARGS <- list(iter = 2000, chains = 4, seed = 12345)
SEED      <- 12345L

# ============================================================
# 1. LIBRARIES
# ============================================================
library(rstan)
library(LCVA)
library(BFL2)
library(knitr)

has_kableExtra <- requireNamespace("kableExtra", quietly = TRUE)
has_gridExtra  <- requireNamespace("gridExtra",  quietly = TRUE)

# ============================================================
# 2. ZOEY SOURCES
# ============================================================
source(file.path(ZOEY_REPO, "model", "assist_functions.R"))
source(file.path(ZOEY_REPO, "model", "execute_balance_mix.R"))

# Verify CSMF_acc implementations are numerically identical before removing
# Zoey's version (which masks BFL2::CSMF_acc after sourcing assist_functions.R)
local({
  pi_true <- c(0.3, 0.5, 0.2)
  pi_hat  <- c(0.25, 0.55, 0.2)
  zo_val  <- CSMF_acc(pi_hat, pi_true)
  bfl_val <- BFL2::CSMF_acc(pi_hat, pi_true)
  if (abs(zo_val - bfl_val) >= 1e-10)
    stop(sprintf(
      "CSMF_acc mismatch: Zoey = %.10f  BFL2 = %.10f — implementations differ!",
      zo_val, bfl_val
    ))
  message(sprintf(
    "CSMF_acc check passed: Zoey = %.8f  BFL2 = %.8f  [OK]", zo_val, bfl_val
  ))
})
rm(list = c("CSMF_acc"))

# ============================================================
# 3. UTILITIES
# ============================================================
`%||%` <- function(a, b) if (!is.null(a)) a else b

train_sites <- setdiff(sites, test_site)

zo_phi_compute <- function(lcva_fit_obj, X_test) {
  out <- LCVA::LCVA.pred(
    fit               = lcva_fit_obj,
    X_test            = as.matrix(X_test),
    model             = "C",
    Burn_in           = 500,
    Nitr              = 4000,
    return_likelihood = TRUE,
    verbose           = FALSE
  )
  apply(out$x_given_y_prob[-c(1:2000), , ], c(2, 3), mean)
}

get_cause_ids <- function(s) {
  as.character(as.numeric(sim_data_filtered_list[[s]]$mapping_new_to_origin))
}

make_ls_entry <- function(phi_rows, cause_ids, row_hash, P_target) {
  list(
    posterior_phi = phi_rows,
    cause_ids     = cause_ids,
    target_info   = list(row_hash = row_hash, N = nrow(phi_rows), P = P_target)
  )
}

align_to_zoey <- function(fit, zo_ca) {
  my_ca <- as.character(fit$causes)
  stopifnot(setequal(zo_ca, my_ca))
  idx <- match(zo_ca, my_ca)
  stopifnot(!anyNA(idx))
  list(
    phi    = fit$phi[,    idx, ,  drop = FALSE],
    pi     = fit$pi[,     idx,    drop = FALSE],
    lambda = fit$lambda[, idx, ,  drop = FALSE],
    idx    = idx
  )
}

bfl2_post_pred <- function(fit, zo_idx, row_subset = NULL) {
  phi_use <- if (is.null(row_subset)) fit$phi else fit$phi[row_subset, , , drop = FALSE]
  pr <- BFL2:::bfl_post_pred_sampler(
    posterior_phi    = phi_use,
    posterior_lambda = fit$lambda,
    posterior_pi     = fit$pi,
    seed             = SEED
  )
  list(
    prob  = pr$posterior_pred_Y_prob_mean[, zo_idx, drop = FALSE],
    draws = pr$posterior_pred_Y
  )
}

zo_to_NC_mean <- function(arr, C) {
  d <- dim(arr)
  stopifnot(length(d) == 3)
  if (d[1] == C) arr <- aperm(arr, c(2, 1, 3))
  apply(arr, c(1, 2), mean)
}

build_param_tbl <- function(label,
                             my_phi, zo_phi,
                             my_pi,  zo_pi_draws,
                             my_lam, zo_lam_draws,
                             my_prob, zo_prob,
                             zo_ca) {
  pi_my <- colMeans(my_pi);       names(pi_my) <- zo_ca
  pi_zo <- colMeans(zo_pi_draws); names(pi_zo) <- zo_ca
  lam_my <- apply(my_lam,       c(2, 3), mean)
  lam_zo <- apply(zo_lam_draws, c(2, 3), mean)
  tbl <- rbind(
    data.frame(experiment = label, component = "phi",
               max_abs_diff  = max(abs(my_phi - zo_phi)),
               mean_abs_diff = mean(abs(my_phi - zo_phi)),
               corr_flat     = cor(as.vector(my_phi), as.vector(zo_phi))),
    data.frame(experiment = label, component = "pi_mean",
               max_abs_diff  = max(abs(pi_my - pi_zo)),
               mean_abs_diff = mean(abs(pi_my - pi_zo)),
               corr_flat     = cor(pi_my, pi_zo)),
    data.frame(experiment = label, component = "lambda_mean",
               max_abs_diff  = max(abs(lam_my - lam_zo)),
               mean_abs_diff = mean(abs(lam_my - lam_zo)),
               corr_flat     = cor(as.vector(lam_my), as.vector(lam_zo))),
    data.frame(experiment = label, component = "postpred_mean",
               max_abs_diff  = max(abs(my_prob - zo_prob)),
               mean_abs_diff = mean(abs(my_prob - zo_prob)),
               corr_flat     = cor(as.vector(my_prob), as.vector(zo_prob)))
  )
  tbl[, 3:5] <- lapply(tbl[, 3:5], function(x) signif(x, 6))
  tbl
}

build_score_tbl <- function(label, zo_obj, my_score) {
  zo_top1 <- zo_obj$top1_acc %||% zo_obj$acc %||% NA_real_
  zo_csmf <- zo_obj$csmf_acc %||% NA_real_
  data.frame(
    experiment = label,
    metric     = c("Top1 Accuracy", "CSMF Accuracy"),
    Zoey       = signif(c(zo_top1, zo_csmf), 6),
    BFL2       = signif(c(my_score$top1_acc, my_score$csmf_acc), 6),
    abs_diff   = signif(c(abs(zo_top1 - my_score$top1_acc),
                           abs(zo_csmf - my_score$csmf_acc)), 6)
  )
}

make_kable <- function(df, caption) {
  tab <- knitr::kable(df, caption = caption, align = "l")
  if (has_kableExtra)
    tab <- tab |>
      kableExtra::kable_styling(full_width = FALSE, position = "left") |>
      kableExtra::column_spec(1:2, bold = TRUE)
  tab
}

save_table_png <- function(df, path, w = 1600, h = 500) {
  if (!has_gridExtra) { message("gridExtra not available; skipping: ", path); return(invisible(NULL)) }
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  png(path, width = w, height = h, res = 200)
  grid::grid.newpage()
  gridExtra::grid.table(df)
  dev.off()
  invisible(path)
}

param_plot <- function(pi_zo, pi_my, lam_zo, lam_my, prob_zo, prob_my,
                        title, lam_subtitle = "λ") {
  op <- par(no.readonly = TRUE); on.exit(par(op))
  par(mfrow = c(1, 3), oma = c(0, 0, 3, 0), mar = c(4, 4, 2, 1))
  plot(pi_zo, pi_my,
       xlab = "Zoey π mean", ylab = "BFL2 π mean", main = "π", pch = 19)
  abline(0, 1, col = "red", lwd = 2)
  plot(as.vector(lam_zo), as.vector(lam_my),
       xlab = "Zoey λ mean", ylab = "BFL2 λ mean", main = lam_subtitle, pch = 19)
  abline(0, 1, col = "red", lwd = 2)
  plot(as.vector(prob_zo), as.vector(prob_my),
       xlab = "Zoey mean p(y|x)", ylab = "BFL2 mean p(y|x)",
       main = "Post. mean p(y|x)", pch = 19)
  abline(0, 1, col = "red", lwd = 2)
  mtext(title, outer = TRUE, cex = 1.1, font = 2)
}

# ============================================================
# 4. LOAD + FILTER DATA
# ============================================================
cat("Loading data...\n")
sim_data_list <- lapply(sites, get_sim_data_from_phmrc_by_location)
names(sim_data_list) <- sites

sim_data_filtered_list <- lapply(sim_data_list, function(sd) {
  filter_sparse_causes(sd, threshold = 0)
})

# ============================================================
# 5. CREATE MISSINGNESS IN TARGET SITE
# ============================================================
set.seed(SEED)
Y_full_truth    <- sim_data_filtered_list[[test_site]]$filtered_data$data.truth$Y.t
missing_indices <- sort(unique(generate_missing_Yt(Y_full_truth, miss_prop = miss_prop)))

sim_data_filtered_list[[test_site]]$filtered_data$data$Y.t <-
  sim_data_filtered_list[[test_site]]$filtered_data$data.truth$Y.t
sim_data_filtered_list[[test_site]]$filtered_data$data$Y.t[missing_indices] <- NA

X_target_full <- as.matrix(sim_data_filtered_list[[test_site]]$filtered_data$data$X)
n0 <- nrow(X_target_full)
P  <- ncol(X_target_full)

cat(sprintf("n0 = %d  |  P = %d  |  missing = %d  (%.0f%%)\n\n",
            n0, P, length(missing_indices),
            100 * length(missing_indices) / n0))

# ============================================================
# 6. FIT LCVA FOR ALL SITES
# ============================================================
cat("Fitting LCVA for all sites...\n")
lcva_fits          <- setNames(vector("list", length(sites)), sites)
posterior_phi_full <- setNames(vector("list", length(sites)), sites)

for (s in sites) {
  cat(sprintf("  [%s] ", s))
  sd      <- sim_data_filtered_list[[s]]$filtered_data
  obs_idx <- which(!is.na(sd$data$Y.t))
  X_tr <- as.matrix(sd$data$X[obs_idx, , drop = FALSE])
  Y_tr <- sd$data$Y.t[obs_idx]

  lcva_fits[[s]] <- LCVA::LCVA.train(
    X      = X_tr, Y = Y_tr,
    Domain = rep(1L, length(Y_tr)),
    K = K, model = "S",
    Nitr = 2000, thin = 2, nchain = 3,
    seed = SEED, verbose = FALSE
  )

  posterior_phi_full[[s]] <- zo_phi_compute(lcva_fits[[s]], X_target_full)
  cat(sprintf("phi %d x %d\n", nrow(posterior_phi_full[[s]]),
              ncol(posterior_phi_full[[s]])))
}

rm(lcva_fits); gc()

# ============================================================
# 7. DOMAIN / PARTIAL SPLIT + SELF-SITE LCVA (Mix-specific)
# ============================================================
split_res       <- split_domain_and_partial(sim_data_filtered_list[[test_site]])
domain_indices  <- split_res$domain_indices
partial_indices <- split_res$partial_indices
test_idx        <- c(partial_indices, missing_indices)   # Mix Stan row-space

cat(sprintf("\ndomain = %d  |  partial = %d  |  missing = %d  |  test_idx = %d\n\n",
            length(domain_indices), length(partial_indices),
            length(missing_indices), length(test_idx)))

# Self-site LCVA trained on domain rows only (Mix uses this for the test site)
cat("Fitting domain-only LCVA for test site (Mix self-site)...\n")
sd_test <- sim_data_filtered_list[[test_site]]$filtered_data

lcva_fit_dom_self <- LCVA::LCVA.train(
  X      = as.matrix(sd_test$data$X[domain_indices, , drop = FALSE]),
  Y      = sd_test$data.truth$Y.t[domain_indices],
  Domain = rep(1L, length(domain_indices)),
  K = K, model = "S",
  Nitr = 2000, thin = 2, nchain = 3,
  seed = SEED, verbose = FALSE
)
phi_dom_self <- zo_phi_compute(lcva_fit_dom_self, X_target_full)
rm(lcva_fit_dom_self); gc()
cat(sprintf("  domain-self phi: %d x %d\n\n", nrow(phi_dom_self), ncol(phi_dom_self)))

# phi collection for Mix: swap test-site entry to domain-only version
posterior_phi_mix              <- posterior_phi_full
posterior_phi_mix[[test_site]] <- phi_dom_self

# ============================================================
# 8. GROUND-TRUTH LABELS
# ============================================================
map_new_to_origin <- sim_data_filtered_list[[test_site]]$mapping_new_to_origin
Y_all_new         <- sim_data_filtered_list[[test_site]]$filtered_data$data.truth$Y.t
Y_all_origin      <- as.character(map_new_to_origin[Y_all_new])

Y_obs_new <- sim_data_filtered_list[[test_site]]$filtered_data$data$Y.t
Y_target_origin_full <- rep(NA_character_, n0)
obs_idx2 <- which(!is.na(Y_obs_new))
Y_target_origin_full[obs_idx2] <- as.character(map_new_to_origin[Y_obs_new[obs_idx2]])

# ============================================================
# ============================================================
# SECTION C: MIX (BALANCED)
# ============================================================
# ============================================================
cat("========================================\n")
cat("SECTION C: MIX (BALANCED)\n")
cat("========================================\n\n")
exp_mix <- "BFL Mix"

# ---- C1. Run Zoey ----
gc(); gc()
zo_mix <- execute_balance_mix(
  test_site              = test_site,
  sites                  = sites,
  posterior_phi_full     = posterior_phi_mix,
  sim_data_filtered_list = sim_data_filtered_list,
  domain_indices         = domain_indices,
  partial_indices        = partial_indices,
  missing_indices        = missing_indices
)

need_mix <- c("causes", "pi_draws", "lambda_draws", "phi_arr",
              "pred_prob", "test_idx", "missing_idx", "stan_data")
missing_mix <- setdiff(need_mix, names(zo_mix))
if (length(missing_mix))
  stop("execute_balance_mix() is missing fields: ",
       paste(missing_mix, collapse = ", "))

zo_ca_mix <- as.character(zo_mix$causes)

stopifnot(identical(sort(zo_mix$missing_idx), sort(missing_indices)))
stopifnot(identical(sort(zo_mix$test_idx),    sort(test_idx)))

N_mix <- length(test_idx)

# ---- C2. Build BFL2 local_summaries (test_idx rows, all sites) ----
X_target_mix <- X_target_full[test_idx, , drop = FALSE]
rh_mix       <- compute_row_hashes(X_target_mix)

ls_mix <- setNames(lapply(sites, function(s) {
  make_ls_entry(
    phi_rows  = posterior_phi_mix[[s]][test_idx, , drop = FALSE],
    cause_ids = get_cause_ids(s),
    row_hash  = rh_mix,
    P_target  = P
  )
}), sites)

# Y_target in test_idx-space: partial rows have labels, missing rows are NA
Y_target_mix <- Y_target_origin_full[test_idx]
stopifnot(sum(!is.na(Y_target_mix)) == length(partial_indices))
stopifnot(sum( is.na(Y_target_mix)) == length(missing_indices))

# ---- C3. Run BFL2 ----
fit_mix <- run_BFL(
  local_summaries = ls_mix,
  X_target        = X_target_mix,
  Y_target        = Y_target_mix,
  stan_args       = STAN_ARGS
)

cat(sprintf("Mix  Zoey N = %d  |  BFL2 n_total = %d  |  BFL2 n_stan = %d\n\n",
            zo_mix$stan_data$N, fit_mix$n_total, length(fit_mix$stan_idx)))
stopifnot(fit_mix$n_total == N_mix)
stopifnot(fit_mix$n_total == zo_mix$stan_data$N)
stopifnot(length(fit_mix$stan_idx) == N_mix)

# ---- C4. Align + compare parameters ----
al_mix <- align_to_zoey(fit_mix, zo_ca_mix)

zo_phi_mix <- zo_mix$phi_arr
if (dim(zo_phi_mix)[2] != length(zo_ca_mix)) {
  zo_phi_mix <- zo_phi_mix[, as.integer(zo_ca_mix), , drop = FALSE]
}

# Trim M dimension to common size (Zoey may store fewer sites)
M_min     <- min(dim(zo_phi_mix)[3], dim(al_mix$phi)[3])
zo_phi_mx <- zo_phi_mix[, , seq_len(M_min), drop = FALSE]
my_phi_mx <- al_mix$phi[, , seq_len(M_min), drop = FALSE]
my_lam_mx <- al_mix$lambda[, , seq_len(M_min), drop = FALSE]
zo_lam_mx <- zo_mix$lambda_draws[, , seq_len(M_min), drop = FALSE]

# Positions of missing_indices inside test_idx (for prob subset)
miss_pos_in_test <- match(missing_indices, test_idx)
stopifnot(!anyNA(miss_pos_in_test))

# Subset phi to miss_pos_in_test rows before sampling to save memory
pp_mix <- bfl2_post_pred(fit_mix, al_mix$idx, row_subset = miss_pos_in_test)
my_prob_mix_miss <- pp_mix$prob

zo_prob_mix <- zo_to_NC_mean(zo_mix$pred_prob, length(zo_ca_mix))
stopifnot(all(dim(my_prob_mix_miss) == dim(zo_prob_mix)))

param_tbl_mix <- build_param_tbl(
  label        = exp_mix,
  my_phi       = my_phi_mx,
  zo_phi       = zo_phi_mx,
  my_pi        = al_mix$pi,
  zo_pi_draws  = zo_mix$pi_draws,
  my_lam       = my_lam_mx,
  zo_lam_draws = zo_lam_mx,
  my_prob      = my_prob_mix_miss,
  zo_prob      = zo_prob_mix,
  zo_ca        = zo_ca_mix
)

# Pre-compute scatter summaries, then free large arrays before score_BFL
pi_zo_mix  <- colMeans(zo_mix$pi_draws); names(pi_zo_mix) <- zo_ca_mix
pi_my_mix  <- colMeans(al_mix$pi);       names(pi_my_mix) <- zo_ca_mix
lam_zo_mix <- apply(zo_lam_mx,       c(2, 3), mean)
lam_my_mix <- apply(my_lam_mx,       c(2, 3), mean)
mn_zo_mix  <- zo_prob_mix
mn_my_mix  <- my_prob_mix_miss
rm(pp_mix, zo_prob_mix, my_prob_mix_miss, zo_phi_mx, my_phi_mx,
   zo_phi_mix, zo_lam_mx, my_lam_mx); gc()

# ---- C5. Score BFL2 ----
# Y_eval must be length fit$n_total = N_mix (test_idx space, all true labels)
Y_eval_mix <- Y_all_origin[test_idx]
score_mix  <- score_BFL(
  fit      = fit_mix,
  Y_eval   = Y_eval_mix,
  eval_idx = miss_pos_in_test,
  seed     = SEED
)

score_tbl_mix <- build_score_tbl(exp_mix, zo_mix, score_mix)

cat("Mix — parameter comparison:\n")
print(make_kable(param_tbl_mix, paste0("Parameter equivalence — ", exp_mix)))
cat("\nMix — score comparison:\n")
print(make_kable(score_tbl_mix, paste0("Accuracy metrics — ", exp_mix)))
cat("\n")

# ---- C6. Save plots + tables ----
out_mix <- file.path(REPORT_DIR, "Mix")
dir.create(out_mix, showWarnings = FALSE, recursive = TRUE)

png(file.path(out_mix, "scatter_params.png"), width = 1800, height = 600, res = 200)
param_plot(pi_zo_mix, pi_my_mix, lam_zo_mix, lam_my_mix,
           mn_zo_mix, mn_my_mix, title = exp_mix, lam_subtitle = "λ (test_idx sites)")
dev.off()
cat("Mix scatter plot saved.\n")

save_table_png(param_tbl_mix, file.path(out_mix, "table_param_summary.png"), w = 1600, h = 500)
save_table_png(score_tbl_mix, file.path(out_mix, "table_scores.png"),        w = 1200, h = 400)

# ============================================================
# SAVE WORKSPACE
# ============================================================
save.image(file.path(out_mix, "bfl2_mix_verify.RData"))
cat("\nSaved workspace:", file.path(out_mix, "bfl2_mix_verify.RData"), "\n")
cat("Done.\n")
