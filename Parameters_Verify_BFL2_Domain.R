############################################################
# BFL2 Parameters Verify — DOMAIN (BALANCED)
#
# Compares BFL2::run_BFL() vs Zoey's execute_balance_domain().
# Standalone: run rm(list=ls()) + gc() then source this file.
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
source(file.path(ZOEY_REPO, "model", "execute_balance_domain.R"))

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
lcva_fits                     <- setNames(vector("list", length(sites)), sites)
posterior_phi_full            <- setNames(vector("list", length(sites)), sites)
LCVA_local_model_test_obs_fit <- list()

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

  if (s == test_site) LCVA_local_model_test_obs_fit[[s]] <- lcva_fits[[s]]

  posterior_phi_full[[s]] <- zo_phi_compute(lcva_fits[[s]], X_target_full)
  cat(sprintf("phi %d x %d\n", nrow(posterior_phi_full[[s]]),
              ncol(posterior_phi_full[[s]])))
}

rm(lcva_fits); gc()

# ============================================================
# 7. GROUND-TRUTH LABELS
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
# SECTION A: DOMAIN (BALANCED)
# ============================================================
# ============================================================
cat("========================================\n")
cat("SECTION A: DOMAIN (BALANCED)\n")
cat("========================================\n\n")
exp_dom <- "BFL Domain"

# ---- A1. Run Zoey ----
zo_dom <- execute_balance_domain(
  test_site                     = test_site,
  sites                         = sites,
  posterior_phi_full            = posterior_phi_full,
  LCVA_local_model_test_obs_fit = LCVA_local_model_test_obs_fit,
  sim_data_filtered_list        = sim_data_filtered_list
)
zo_ca_dom <- as.character(zo_dom$causes)

# ---- A2. Build BFL2 local_summaries (missing rows only, all sites) ----
rh_dom <- compute_row_hashes(X_target_full[missing_indices, , drop = FALSE])

ls_dom <- setNames(lapply(sites, function(s) {
  make_ls_entry(
    phi_rows  = posterior_phi_full[[s]][missing_indices, , drop = FALSE],
    cause_ids = get_cause_ids(s),
    row_hash  = rh_dom,
    P_target  = P
  )
}), sites)

# ---- A3. Run BFL2 ----
fit_dom <- run_BFL(
  local_summaries = ls_dom,
  X_target        = X_target_full,
  Y_target        = Y_target_origin_full,
  stan_args       = STAN_ARGS
)

cat(sprintf("Domain  Zoey N = %d  |  BFL2 n_total = %d  |  BFL2 n_stan = %d\n\n",
            nrow(zo_dom$phi_arr), fit_dom$n_total, length(fit_dom$stan_idx)))
stopifnot(length(fit_dom$stan_idx) == length(missing_indices))

# ---- A4. Align + compare parameters ----
al_dom <- align_to_zoey(fit_dom, zo_ca_dom)

zo_phi_dom <- zo_dom$phi_arr
if (dim(zo_phi_dom)[2] != length(zo_ca_dom)) {
  zo_phi_dom <- zo_phi_dom[, as.integer(zo_ca_dom), , drop = FALSE]
}

pp_dom     <- bfl2_post_pred(fit_dom, al_dom$idx)
zo_prob_dom <- zo_to_NC_mean(zo_dom$pred_prob, length(zo_ca_dom))

stopifnot(all(dim(pp_dom$prob) == dim(zo_prob_dom)))

param_tbl_dom <- build_param_tbl(
  label        = exp_dom,
  my_phi       = al_dom$phi,
  zo_phi       = zo_phi_dom,
  my_pi        = al_dom$pi,
  zo_pi_draws  = zo_dom$pi_draws,
  my_lam       = al_dom$lambda,
  zo_lam_draws = zo_dom$lambda_draws,
  my_prob      = pp_dom$prob,
  zo_prob      = zo_prob_dom,
  zo_ca        = zo_ca_dom
)

# Pre-compute scatter summaries, then free large arrays before score_BFL
pi_zo_dom  <- colMeans(zo_dom$pi_draws); names(pi_zo_dom) <- zo_ca_dom
pi_my_dom  <- colMeans(al_dom$pi);       names(pi_my_dom) <- zo_ca_dom
lam_zo_dom <- apply(zo_dom$lambda_draws, c(2, 3), mean)
lam_my_dom <- apply(al_dom$lambda,       c(2, 3), mean)
mn_zo_dom  <- zo_prob_dom
mn_my_dom  <- pp_dom$prob
rm(pp_dom, zo_prob_dom, zo_phi_dom); gc()

# ---- A5. Score BFL2 ----
score_dom <- score_BFL(
  fit      = fit_dom,
  Y_eval   = Y_all_origin,
  eval_idx = missing_indices,
  seed     = SEED
)

score_tbl_dom <- build_score_tbl(exp_dom, zo_dom, score_dom)

cat("Domain — parameter comparison:\n")
print(make_kable(param_tbl_dom, paste0("Parameter equivalence — ", exp_dom)))
cat("\nDomain — score comparison:\n")
print(make_kable(score_tbl_dom, paste0("Accuracy metrics — ", exp_dom)))
cat("\n")

# ---- A6. Save plots + tables ----
out_dom <- file.path(REPORT_DIR, "Domain")
dir.create(out_dom, showWarnings = FALSE, recursive = TRUE)

png(file.path(out_dom, "scatter_params.png"), width = 1800, height = 600, res = 200)
param_plot(pi_zo_dom, pi_my_dom, lam_zo_dom, lam_my_dom,
           mn_zo_dom, mn_my_dom, title = exp_dom, lam_subtitle = "λ (all sites)")
dev.off()
cat("Domain scatter plot saved.\n")

save_table_png(param_tbl_dom, file.path(out_dom, "table_param_summary.png"), w = 1600, h = 500)
save_table_png(score_tbl_dom, file.path(out_dom, "table_scores.png"),        w = 1200, h = 400)

# ============================================================
# SAVE WORKSPACE
# ============================================================
save.image(file.path(out_dom, "bfl2_domain_verify.RData"))
cat("\nSaved workspace:", file.path(out_dom, "bfl2_domain_verify.RData"), "\n")
cat("Done.\n")
