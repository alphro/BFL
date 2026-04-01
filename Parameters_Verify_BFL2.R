############################################################
# BFL2 Parameters Verify — Domain / Partial / Mix
#
# Compares BFL2::run_BFL() posterior estimates against
# Zoey's execute_unbalanced_domain / execute_balance_partial /
# execute_balance_mix.
#
# Design principle: phi matrices are computed identically for
# both sides (same LCVA calls, same burn-in/drop) so any
# differences reflect Stan model behavior only.
#
# BFL2 API changes from the old verify scripts:
#   BFL:::hash_rows(X)              -> compute_row_hashes(X)
#   run_BFL(..., model = "X")       -> run_BFL(...)  [auto-inferred]
#   BFL:::bfl_post_pred_sampler()   -> BFL2:::bfl_post_pred_sampler()
#   score_BFL(pred, Y, miss, stan)  -> score_BFL(fit, Y_eval, eval_idx)
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
source(file.path(ZOEY_REPO, "model", "execute_balance_partial.R"))
source(file.path(ZOEY_REPO, "model", "execute_balance_mix.R"))

# Verify CSMF_acc implementations are numerically identical before removing
# Zoey's version (which masks BFL2::CSMF_acc after sourcing assist_functions.R)
local({
  pi_true <- c(0.3, 0.5, 0.2)
  pi_hat  <- c(0.25, 0.55, 0.2)
  zo_val  <- CSMF_acc(pi_hat, pi_true)           # Zoey's (currently active in env)
  bfl_val <- BFL2::CSMF_acc(pi_hat, pi_true)     # BFL2's (explicit namespace)
  if (abs(zo_val - bfl_val) >= 1e-10)
    stop(sprintf(
      "CSMF_acc mismatch: Zoey = %.10f  BFL2 = %.10f — implementations differ!",
      zo_val, bfl_val
    ))
  message(sprintf(
    "CSMF_acc check passed: Zoey = %.8f  BFL2 = %.8f  [OK]", zo_val, bfl_val
  ))
})
rm(list = c("CSMF_acc"))   # safe to drop — verified equivalent above

# ============================================================
# 3. UTILITIES
# ============================================================
train_sites <- setdiff(sites, test_site)

# phi computation matching Zoey exactly:
#   LCVA.pred(Burn_in=500, Nitr=4000) then drop first 2000 draws
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
  apply(out$x_given_y_prob[-c(1:2000), , ], c(2, 3), mean)  # N x C_site
}

# cause_ids for site s using Zoey's original-ID convention
get_cause_ids <- function(s) {
  as.character(as.numeric(sim_data_filtered_list[[s]]$mapping_new_to_origin))
}

# build a single local_summaries entry (BFL2 format)
make_ls_entry <- function(phi_rows, cause_ids, row_hash, P_target) {
  list(
    posterior_phi = phi_rows,
    cause_ids     = cause_ids,
    target_info   = list(
      row_hash = row_hash,
      N        = nrow(phi_rows),
      P        = P_target
    )
  )
}

# align BFL2 fit arrays (phi, pi, lambda) to Zoey's cause order
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

# posterior predictive MEAN probs from a BFL2 fit, C-dim reordered to Zoey.
# bfl_post_pred_sampler now returns posterior_pred_Y_prob_mean (I x C) rather
# than the full I x C x S array, so peak memory is O(I*C + S*I) not O(I*C*S).
# row_subset: optional integer vector — subset phi to these rows only.
bfl2_post_pred <- function(fit, zo_idx, row_subset = NULL) {
  phi_use <- if (is.null(row_subset)) fit$phi else fit$phi[row_subset, , , drop = FALSE]
  pr <- BFL2:::bfl_post_pred_sampler(
    posterior_phi    = phi_use,
    posterior_lambda = fit$lambda,
    posterior_pi     = fit$pi,
    seed             = SEED
  )
  # prob is now I x C (posterior mean) — no S dimension
  list(
    prob  = pr$posterior_pred_Y_prob_mean[, zo_idx, drop = FALSE],  # I x C
    draws = pr$posterior_pred_Y                                       # S x I
  )
}

# collapse Zoey pred_prob (N x C x S or C x N x S) to N x C mean
zo_to_NC_mean <- function(arr, C) {
  d <- dim(arr)
  stopifnot(length(d) == 3)
  if (d[1] == C) arr <- aperm(arr, c(2, 1, 3))  # -> N x C x S
  apply(arr, c(1, 2), mean)                       # -> N x C
}

# build parameter comparison table (phi, pi, lambda, postpred)
# my_prob and zo_prob are both N x C posterior-mean matrices (NOT N x C x S)
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

  # my_prob and zo_prob already N x C means — use directly
  mn_my <- my_prob
  mn_zo <- zo_prob

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
               max_abs_diff  = max(abs(mn_my - mn_zo)),
               mean_abs_diff = mean(abs(mn_my - mn_zo)),
               corr_flat     = cor(as.vector(mn_my), as.vector(mn_zo)))
  )
  tbl[, 3:5] <- lapply(tbl[, 3:5], function(x) signif(x, 6))
  tbl
}

# build score comparison table (top1 / csmf, Zoey vs BFL2)
build_score_tbl <- function(label, zo_obj, my_score) {
  zo_top1 <- zo_obj$top1_acc %||% zo_obj$acc     %||% NA_real_
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
`%||%` <- function(a, b) if (!is.null(a)) a else b

# kable wrapper
make_kable <- function(df, caption) {
  tab <- knitr::kable(df, caption = caption, align = "l")
  if (has_kableExtra)
    tab <- tab |>
      kableExtra::kable_styling(full_width = FALSE, position = "left") |>
      kableExtra::column_spec(1:2, bold = TRUE)
  tab
}

# save a data.frame as a PNG table
save_table_png <- function(df, path, w = 1600, h = 500) {
  if (!has_gridExtra) { message("gridExtra not available; skipping: ", path); return(invisible(NULL)) }
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  png(path, width = w, height = h, res = 200)
  grid::grid.newpage()
  gridExtra::grid.table(df)
  dev.off()
  invisible(path)
}

# 1x3 scatter plot: pi, lambda, posterior-pred prob
param_plot <- function(pi_zo, pi_my,
                        lam_zo, lam_my,
                        prob_zo, prob_my,
                        title, lam_subtitle = "λ") {
  op <- par(no.readonly = TRUE); on.exit(par(op))
  par(mfrow = c(1, 3), oma = c(0, 0, 3, 0), mar = c(4, 4, 2, 1))

  plot(pi_zo,  pi_my,
       xlab = "Zoey π mean", ylab = "BFL2 π mean",
       main = "π", pch = 19)
  abline(0, 1, col = "red", lwd = 2)

  plot(as.vector(lam_zo), as.vector(lam_my),
       xlab = "Zoey λ mean", ylab = "BFL2 λ mean",
       main = lam_subtitle, pch = 19)
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

# Stamp observed Y onto target then NA out missing rows
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
# 6. FIT LCVA FOR ALL SITES (observed rows → predict on X_full)
#    phi computed with Zoey-identical burn-in/drop
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
    X      = X_tr,
    Y      = Y_tr,
    Domain = rep(1L, length(Y_tr)),
    K      = K, model = "S",
    Nitr   = 2000, thin = 2, nchain = 3,
    seed   = SEED, verbose = FALSE
  )

  # Both execute_balance_domain and execute_unbalanced_domain require this
  if (s == test_site) LCVA_local_model_test_obs_fit[[s]] <- lcva_fits[[s]]

  posterior_phi_full[[s]] <- zo_phi_compute(lcva_fits[[s]], X_target_full)
  cat(sprintf("phi %d x %d\n", nrow(posterior_phi_full[[s]]),
              ncol(posterior_phi_full[[s]])))
}

# Free LCVA fit objects — posterior_phi_full is already computed,
# raw fits are no longer needed and are large.
rm(lcva_fits); gc()

# ============================================================
# 7. DOMAIN / PARTIAL SPLIT (needed for Mix row-space)
# ============================================================
split_res       <- split_domain_and_partial(sim_data_filtered_list[[test_site]])
domain_indices  <- split_res$domain_indices
partial_indices <- split_res$partial_indices
test_idx        <- c(partial_indices, missing_indices)   # Mix Stan row-space

cat(sprintf("\ndomain = %d  |  partial = %d  |  missing = %d  |  test_idx = %d\n\n",
            length(domain_indices), length(partial_indices),
            length(missing_indices), length(test_idx)))

# Test-site LCVA trained on domain rows ONLY (Mix self-site)
cat("Fitting domain-only LCVA for test site (Mix self-site)...\n")
sd_test <- sim_data_filtered_list[[test_site]]$filtered_data

lcva_fit_dom_self <- LCVA::LCVA.train(
  X      = as.matrix(sd_test$data$X[domain_indices, , drop = FALSE]),
  Y      = sd_test$data.truth$Y.t[domain_indices],
  Domain = rep(1L, length(domain_indices)),
  K      = K, model = "S",
  Nitr   = 2000, thin = 2, nchain = 3,
  seed   = SEED, verbose = FALSE
)
phi_dom_self <- zo_phi_compute(lcva_fit_dom_self, X_target_full)  # n0 x C_test
cat(sprintf("  domain-self phi: %d x %d\n\n",
            nrow(phi_dom_self), ncol(phi_dom_self)))

# phi collection for Mix: swap test-site phi to domain-only version
posterior_phi_mix                  <- posterior_phi_full
posterior_phi_mix[[test_site]]     <- phi_dom_self

# ============================================================
# 8. GROUND-TRUTH LABELS
# ============================================================
map_new_to_origin <- sim_data_filtered_list[[test_site]]$mapping_new_to_origin
Y_all_new         <- sim_data_filtered_list[[test_site]]$filtered_data$data.truth$Y.t
Y_all_origin      <- as.character(map_new_to_origin[Y_all_new])   # length n0, no NAs

# Y with NAs at missing rows (for run_BFL Y_target)
Y_obs_new <- sim_data_filtered_list[[test_site]]$filtered_data$data$Y.t
Y_target_origin_full <- rep(NA_character_, n0)
obs_idx2 <- which(!is.na(Y_obs_new))
Y_target_origin_full[obs_idx2] <- as.character(map_new_to_origin[Y_obs_new[obs_idx2]])

# ============================================================
# ============================================================
# SECTION A: DOMAIN
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
# Domain design (Zoey-aligned):
#   - local_summaries has phi for unlabeled (missing) rows only
#   - X_target = X_full (all N rows), Y_target = observed labels + NAs
#   - run_BFL infers stan_idx = missing_indices via row hashing
#   - CSMF correction (n_stan * pi_hat + nLc) / n_total applied automatically
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

zo_phi_dom <- zo_dom$phi_arr   # N_miss x C_max x M — subset C if needed
if (dim(zo_phi_dom)[2] != length(zo_ca_dom)) {
  zo_phi_dom <- zo_phi_dom[, as.integer(zo_ca_dom), , drop = FALSE]
}

pp_dom     <- bfl2_post_pred(fit_dom, al_dom$idx)          # N_miss x C x S
zo_prob_dom <- zo_to_NC_mean(zo_dom$pred_prob, length(zo_ca_dom))  # N_miss x C (mean)

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

# Pre-compute scatter plot summaries from pp_dom now, then free it before
# score_BFL allocates its own N_miss x C x S array internally.
pi_zo_dom  <- colMeans(zo_dom$pi_draws); names(pi_zo_dom) <- zo_ca_dom
pi_my_dom  <- colMeans(al_dom$pi);       names(pi_my_dom) <- zo_ca_dom
lam_zo_dom <- apply(zo_dom$lambda_draws, c(2, 3), mean)
lam_my_dom <- apply(al_dom$lambda,       c(2, 3), mean)
mn_zo_dom  <- zo_prob_dom   # already N_miss x C mean
mn_my_dom  <- apply(pp_dom$prob,         c(1, 2), mean)  # N_miss x C -> small
rm(pp_dom, zo_prob_dom, zo_phi_dom); gc()

# ---- A5. Score BFL2 ----
# Y_eval = Y_all_origin (length n_total = n0)
# eval_idx = missing_indices (positions to evaluate top-1 on)
# CSMF correction applied automatically inside score_BFL
score_dom <- score_BFL(
  fit      = fit_dom,
  Y_eval   = Y_all_origin,
  eval_idx = missing_indices,
  seed     = SEED
)

score_tbl_dom <- build_score_tbl(exp_dom, zo_dom, score_dom)

cat("Domain — parameter comparison:\n")
print(param_tbl_dom)
cat("\nDomain — score comparison:\n")
print(score_tbl_dom)
cat("\n")

# ---- A6. Plot + free Domain memory ----
dir.create(file.path(REPORT_DIR, "Domain"), showWarnings = FALSE, recursive = TRUE)
png(file.path(REPORT_DIR, "Domain", "scatter_params.png"), width = 1800, height = 600, res = 200)
param_plot(pi_zo_dom, pi_my_dom, lam_zo_dom, lam_my_dom,
           mn_zo_dom, mn_my_dom, title = exp_dom, lam_subtitle = "λ (all sites)")
dev.off()
cat("Domain scatter plot saved.\n")

rm(zo_dom, fit_dom, al_dom, ls_dom, rh_dom,
   pi_zo_dom, pi_my_dom, lam_zo_dom, lam_my_dom, mn_zo_dom, mn_my_dom)
gc(); gc()

# ============================================================
# ============================================================
# SECTION B: PARTIAL
# ============================================================
# ============================================================
cat("========================================\n")
cat("SECTION B: PARTIAL\n")
cat("========================================\n\n")
exp_par <- "BFL Partial"

# ---- B1. Run Zoey ----
gc(); gc()   # clear Domain residuals before Zoey's heavy LCVA.pred calls
zo_par <- execute_balance_partial(
  test_site              = test_site,
  sites                  = sites,
  posterior_phi_full     = posterior_phi_full,
  sim_data_filtered_list = sim_data_filtered_list
)

need_par <- c("causes", "pi_draws", "lambda_draws", "phi_arr", "pred_prob", "missing_idx")
missing_par <- setdiff(need_par, names(zo_par))
if (length(missing_par))
  stop("execute_balance_partial() is missing fields: ",
       paste(missing_par, collapse = ", "),
       "\nModify it to return these (see original verify script comments).")

zo_ca_par   <- as.character(zo_par$causes)
zo_miss_par <- zo_par$missing_idx   # row indices (in 1..n0) Zoey scored

# ---- B2. Build BFL2 local_summaries (train_sites only, all N rows) ----
rh_par <- compute_row_hashes(X_target_full)

ls_par <- setNames(lapply(train_sites, function(s) {
  make_ls_entry(
    phi_rows  = posterior_phi_full[[s]],   # all n0 rows
    cause_ids = get_cause_ids(s),
    row_hash  = rh_par,
    P_target  = P
  )
}), train_sites)

# ---- B3. Run BFL2 ----
# Partial design:
#   - train_sites only (no self-site)
#   - local_summaries has phi for all N rows
#   - Y_target has observed labels + NAs for missing
#   - Stan uses partial-label model; stan_idx = 1..n0
fit_par <- run_BFL(
  local_summaries = ls_par,
  X_target        = X_target_full,
  Y_target        = Y_target_origin_full,
  stan_args       = STAN_ARGS
)

cat(sprintf("Partial  Zoey N = %d  |  BFL2 n_total = %d  |  BFL2 n_stan = %d\n\n",
            nrow(zo_par$phi_arr), fit_par$n_total, length(fit_par$stan_idx)))
stopifnot(fit_par$n_total == n0)
stopifnot(length(fit_par$stan_idx) == n0)

# ---- B4. Align + compare parameters ----
al_par <- align_to_zoey(fit_par, zo_ca_par)

zo_phi_par <- zo_par$phi_arr
if (dim(zo_phi_par)[2] != length(zo_ca_par)) {
  zo_phi_par <- zo_phi_par[, as.integer(zo_ca_par), , drop = FALSE]
}

# BFL2 post-pred: subset phi to zo_miss_par rows BEFORE sampling to avoid OOM.
# phi is 1554 x C x M; computing for all N then subsetting allocates I x C x S
# = 1554 x 34 x 4000 ~ 1.7 GB intermediate. Subsetting first reduces I to
# length(zo_miss_par) ~ 310 rows, cutting memory ~5x.
pp_par <- bfl2_post_pred(fit_par, al_par$idx, row_subset = zo_miss_par)
my_prob_par_miss <- pp_par$prob   # already length(zo_miss_par) x C x S

zo_prob_par <- zo_to_NC_mean(zo_par$pred_prob, length(zo_ca_par))  # N_miss x C (mean)
stopifnot(all(dim(my_prob_par_miss) == dim(zo_prob_par)))

param_tbl_par <- build_param_tbl(
  label        = exp_par,
  my_phi       = al_par$phi,
  zo_phi       = zo_phi_par,
  my_pi        = al_par$pi,
  zo_pi_draws  = zo_par$pi_draws,
  my_lam       = al_par$lambda,
  zo_lam_draws = zo_par$lambda_draws,
  my_prob      = my_prob_par_miss,
  zo_prob      = zo_prob_par,
  zo_ca        = zo_ca_par
)

# Pre-compute scatter summaries then free large arrays before score_BFL
pi_zo_par  <- colMeans(zo_par$pi_draws); names(pi_zo_par) <- zo_ca_par
pi_my_par  <- colMeans(al_par$pi);       names(pi_my_par) <- zo_ca_par
lam_zo_par <- apply(zo_par$lambda_draws, c(2, 3), mean)
lam_my_par <- apply(al_par$lambda,       c(2, 3), mean)
mn_zo_par  <- zo_prob_par   # already N_miss x C mean
mn_my_par  <- apply(my_prob_par_miss,    c(1, 2), mean)
rm(pp_par, zo_prob_par, zo_phi_par, my_prob_par_miss); gc()

# ---- B5. Score BFL2 ----
score_par <- score_BFL(
  fit      = fit_par,
  Y_eval   = Y_all_origin,
  eval_idx = missing_indices,
  seed     = SEED
)

score_tbl_par <- build_score_tbl(exp_par, zo_par, score_par)

cat("Partial — parameter comparison:\n")
print(param_tbl_par)
cat("\nPartial — score comparison:\n")
print(score_tbl_par)
cat("\n")

# ---- B6. Plot + free Partial memory ----
dir.create(file.path(REPORT_DIR, "Partial"), showWarnings = FALSE, recursive = TRUE)
png(file.path(REPORT_DIR, "Partial", "scatter_params.png"), width = 1800, height = 600, res = 200)
param_plot(pi_zo_par, pi_my_par, lam_zo_par, lam_my_par,
           mn_zo_par, mn_my_par, title = exp_par, lam_subtitle = "λ (train sites)")
dev.off()
cat("Partial scatter plot saved.\n")

rm(zo_par, fit_par, al_par, ls_par, rh_par,
   pi_zo_par, pi_my_par, lam_zo_par, lam_my_par, mn_zo_par, mn_my_par)
gc(); gc()

# ============================================================
# ============================================================
# SECTION C: MIX
# ============================================================
# ============================================================
cat("========================================\n")
cat("SECTION C: MIX\n")
cat("========================================\n\n")
exp_mix <- "BFL Mix"

# ---- C1. Run Zoey ----
gc(); gc()   # clear Partial residuals before Zoey's heavy LCVA.pred calls
# posterior_phi_mix has domain-only phi for test_site; all-observed phi for sources
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

# Sanity: Zoey and we agree on the shared splits
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

# Y_target in test_idx-space:
#   - partial rows have observed labels
#   - missing rows are NA
Y_target_mix <- Y_target_origin_full[test_idx]
stopifnot(sum(!is.na(Y_target_mix)) == length(partial_indices))
stopifnot(sum( is.na(Y_target_mix)) == length(missing_indices))

# ---- C3. Run BFL2 ----
# Mix design:
#   - X_target = test_idx rows (partial + missing; domain rows excluded)
#   - Y_target has labels for partial, NA for missing
#   - Stan uses partial-label model; stan_idx = 1..N_mix
#   - No CSMF correction (Stan sees all N_mix rows; domain rows out of scope)
fit_mix <- run_BFL(
  local_summaries = ls_mix,
  X_target        = X_target_mix,
  Y_target        = Y_target_mix,
  stan_args       = STAN_ARGS
)

cat(sprintf("Mix  Zoey N = %d  |  BFL2 n_total = %d  |  BFL2 n_stan = %d\n\n",
            zo_mix$stan_data$N, fit_mix$n_total, length(fit_mix$stan_idx)))
stopifnot(fit_mix$n_total    == N_mix)
stopifnot(fit_mix$n_total    == zo_mix$stan_data$N)
stopifnot(length(fit_mix$stan_idx) == N_mix)

# ---- C4. Align + compare parameters ----
al_mix <- align_to_zoey(fit_mix, zo_ca_mix)

zo_phi_mix <- zo_mix$phi_arr
if (dim(zo_phi_mix)[2] != length(zo_ca_mix)) {
  zo_phi_mix <- zo_phi_mix[, as.integer(zo_ca_mix), , drop = FALSE]
}
# Trim M dimension to common size (Zoey may store fewer sites)
M_min     <- min(dim(zo_phi_mix)[3], dim(al_mix$phi)[3])
zo_phi_mx <- zo_phi_mix[,  , seq_len(M_min), drop = FALSE]
my_phi_mx <- al_mix$phi[,  , seq_len(M_min), drop = FALSE]
my_lam_mx <- al_mix$lambda[, , seq_len(M_min), drop = FALSE]
zo_lam_mx <- zo_mix$lambda_draws[, , seq_len(M_min), drop = FALSE]

# Positions of missing_indices inside test_idx (for prob subset)
miss_pos_in_test <- match(missing_indices, test_idx)
stopifnot(!anyNA(miss_pos_in_test))

# BFL2 post-pred: subset phi to miss_pos_in_test rows BEFORE sampling to save memory
pp_mix <- bfl2_post_pred(fit_mix, al_mix$idx, row_subset = miss_pos_in_test)

my_prob_mix_miss <- pp_mix$prob   # already N_miss x C x S (no further subsetting needed)
zo_prob_mix      <- zo_to_NC_mean(zo_mix$pred_prob, length(zo_ca_mix))  # N_miss x C (mean)
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

# Pre-compute scatter summaries then free large arrays before score_BFL
pi_zo_mix  <- colMeans(zo_mix$pi_draws); names(pi_zo_mix) <- zo_ca_mix
pi_my_mix  <- colMeans(al_mix$pi);       names(pi_my_mix) <- zo_ca_mix
lam_zo_mix <- apply(zo_lam_mx,       c(2, 3), mean)
lam_my_mix <- apply(my_lam_mx,       c(2, 3), mean)
mn_zo_mix  <- zo_prob_mix   # already N_miss x C mean
mn_my_mix  <- apply(my_prob_mix_miss, c(1, 2), mean)
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
print(param_tbl_mix)
cat("\nMix — score comparison:\n")
print(score_tbl_mix)
cat("\n")

# ---- C6. Plot + free Mix memory ----
dir.create(file.path(REPORT_DIR, "Mix"), showWarnings = FALSE, recursive = TRUE)
png(file.path(REPORT_DIR, "Mix", "scatter_params.png"), width = 1800, height = 600, res = 200)
param_plot(pi_zo_mix, pi_my_mix, lam_zo_mix, lam_my_mix,
           mn_zo_mix, mn_my_mix, title = exp_mix, lam_subtitle = "λ (test_idx sites)")
dev.off()
cat("Mix scatter plot saved.\n")

rm(zo_mix, fit_mix, al_mix, pp_mix, ls_mix,
   zo_phi_mix, zo_phi_mx, my_phi_mx, zo_prob_mix, my_prob_mix_miss,
   my_lam_mx, zo_lam_mx,
   pi_zo_mix, pi_my_mix, lam_zo_mix, lam_my_mix, mn_zo_mix, mn_my_mix)
gc(); gc()

# ============================================================
# COMBINED REPORT
# ============================================================
cat("========================================\n")
cat("COMBINED REPORT\n")
cat("========================================\n\n")

all_param_tbl <- rbind(param_tbl_dom, param_tbl_par, param_tbl_mix)
all_score_tbl <- rbind(score_tbl_dom, score_tbl_par, score_tbl_mix)

print(make_kable(all_param_tbl, "Parameter equivalence — Domain / Partial / Mix"))
cat("\n")
print(make_kable(all_score_tbl, "Accuracy metrics — Domain / Partial / Mix"))

# Save combined tables
dir.create(REPORT_DIR, showWarnings = FALSE, recursive = TRUE)
save_table_png(all_param_tbl,
               file.path(REPORT_DIR, "table_all_param_summary.png"),
               w = 2000, h = 700)
save_table_png(all_score_tbl,
               file.path(REPORT_DIR, "table_all_scores.png"),
               w = 1600, h = 600)

# Per-variant sub-directories
for (v in list(
  list(param = param_tbl_dom, score = score_tbl_dom, name = "Domain"),
  list(param = param_tbl_par, score = score_tbl_par, name = "Partial"),
  list(param = param_tbl_mix, score = score_tbl_mix, name = "Mix")
)) {
  sub <- file.path(REPORT_DIR, v$name)
  save_table_png(v$param, file.path(sub, "table_param_summary.png"),  w = 1600, h = 500)
  save_table_png(v$score, file.path(sub, "table_scores.png"),         w = 1200, h = 400)
}

# Scatter plots already saved per-section (A6 / B6 / C6) to conserve memory.

# ============================================================
# SAVE WORKSPACE
# ============================================================
save.image(file.path(REPORT_DIR, "bfl2_params_verify_all.RData"))
cat("\nSaved workspace:", file.path(REPORT_DIR, "bfl2_params_verify_all.RData"), "\n")
cat("Done.\n")
