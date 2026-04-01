# ============================================================
# Quick BFL2 test — Domain / Partial / Mix
# BFL2 (P(X|Y)) vs BFL_yz (P(X|Y,Z)), single seed, single points
#
# To run 3-seed boxplot version instead, see run_one_bfl2() below
# (currently commented out).
# ============================================================
library(devtools); load_all()
library(dplyr); library(tidyr); library(ggplot2); library(LCVA)
options(mc.cores = 1)

# Timing helpers
tic <- function(lbl) { message(sprintf("  [timer] %s ...", lbl)); proc.time() }
toc <- function(t0, lbl) {
  s <- round((proc.time() - t0)[["elapsed"]], 1)
  message(sprintf("  [timer] %s done  (%.1f s)", lbl, s)); s
}

# Settings
test_site    <- "AP"
seed_run     <- 1001L
sites_all    <- c("Mexico", "AP", "Bohol", "Dar", "Pemba", "UP")
source_sites <- setdiff(sites_all, test_site)
lcva_args    <- list(K = 5, Nitr = 1000, thin = 2, pred_Nitr = 1000, seed = 12345)
stan_args    <- list(iter = 1000, chains = 2, seed = 12345)

phmrc <- read.csv(system.file("extdata", "phmrc_clean.csv", package = "BFL"),
                  stringsAsFactors = FALSE)
`%||%` <- function(a, b) if (!is.null(a)) a else b
timing_log <- list()

# ============================================================
# PART 1 — BFL2: Domain / Partial / Mix
# ============================================================
message("=== BFL2 (seed ", seed_run, ") ===")
t_part1 <- tic("BFL2 total")
set.seed(seed_run)

df_site       <- phmrc[phmrc$site == test_site, ]
n             <- nrow(df_site)
idx_labeled   <- sample.int(n, floor(0.20 * n))
idx_unlabeled <- setdiff(seq_len(n), idx_labeled)
X_full        <- as.matrix(df_site[, 1:168])
Y_full        <- df_site$cause
X_labeled     <- X_full[idx_labeled,   , drop = FALSE]
Y_labeled     <- Y_full[idx_labeled]
X_unlabeled   <- X_full[idx_unlabeled, , drop = FALSE]
Y_partial     <- replace(Y_full, idx_unlabeled, NA_character_)
X_train_list  <- setNames(lapply(source_sites, function(s) as.matrix(phmrc[phmrc$site==s, 1:168])), source_sites)
Y_train_list  <- setNames(lapply(source_sites, function(s) phmrc[phmrc$site==s, "cause"]),           source_sites)
sa            <- modifyList(stan_args, list(seed = seed_run))

# Domain
t0 <- tic("BFL2 Domain — LCVA src")
local_dom_src  <- run_lcva(X_train_list, Y_train_list, X_unlabeled, lcva_args)
timing_log[["BFL2_Domain_lcva_src"]]  <- toc(t0, "BFL2 Domain — LCVA src")
t0 <- tic("BFL2 Domain — LCVA self")
local_dom_self <- run_lcva(X_labeled, Y_labeled, X_unlabeled, lcva_args)
timing_log[["BFL2_Domain_lcva_self"]] <- toc(t0, "BFL2 Domain — LCVA self")
t0 <- tic("BFL2 Domain — Stan")
fit_domain <- run_BFL(c(local_dom_src, setNames(list(local_dom_self), test_site)),
                      X_full, Y_target = Y_partial, stan_args = sa)
timing_log[["BFL2_Domain_stan"]] <- toc(t0, "BFL2 Domain — Stan")
sc_domain <- score_BFL(fit_domain, Y_full, eval_idx = idx_unlabeled)

# Partial
t0 <- tic("BFL2 Partial — LCVA src")
local_partial <- run_lcva(X_train_list, Y_train_list, X_full, lcva_args)
timing_log[["BFL2_Partial_lcva_src"]] <- toc(t0, "BFL2 Partial — LCVA src")
t0 <- tic("BFL2 Partial — Stan")
fit_partial <- run_BFL(local_partial, X_full, Y_target = Y_partial,
                       label_shift = FALSE, stan_args = sa)
timing_log[["BFL2_Partial_stan"]] <- toc(t0, "BFL2 Partial — Stan")
sc_partial <- score_BFL(fit_partial, Y_full, eval_idx = idx_unlabeled)

# Mix
t0 <- tic("BFL2 Mix — LCVA self")
local_mix_self <- run_lcva(X_labeled, Y_labeled, X_full, lcva_args)
timing_log[["BFL2_Mix_lcva_self"]] <- toc(t0, "BFL2 Mix — LCVA self")
t0 <- tic("BFL2 Mix — Stan")
fit_mix <- run_BFL(c(local_partial, setNames(list(local_mix_self), test_site)),
                   X_full, Y_target = Y_partial, label_shift = FALSE, stan_args = sa)
timing_log[["BFL2_Mix_stan"]] <- toc(t0, "BFL2 Mix — Stan")
sc_mix <- score_BFL(fit_mix, Y_full, eval_idx = idx_unlabeled)

timing_log[["BFL2_total"]] <- toc(t_part1, "BFL2 total")

res_std <- data.frame(
  site = test_site, seed = seed_run, type = "BFL2",
  method = c("Domain", "Partial", "Mix"),
  csmf = c(sc_domain$csmf_acc, sc_partial$csmf_acc, sc_mix$csmf_acc),
  top1 = c(sc_domain$top1_acc, sc_partial$top1_acc, sc_mix$top1_acc)
)

# ============================================================
# [COMMENTED OUT] 3-seed boxplot version
# Uncomment seeds + run_one_bfl2 block to enable.
# ============================================================
# seeds <- 1001:1003
# run_one_bfl2 <- function(seed_run, ...) { ... }
# res_std <- bind_rows(lapply(seeds, function(s) run_one_bfl2(s, ...)))

# ============================================================
# PART 2a — BFL_yz LCVA fits (cached to disk)
#
# Fitting LCVA is slow (~3 min).  Results are saved so Part 2b
# (Stan runs) can be re-run without re-fitting LCVA.
# To re-run only Stan: load the cache file manually then run Part 2b.
# ============================================================
message("\n=== BFL_yz — LCVA fits (seed ", seed_run, ") ===")

bfl_yz_cache <- sprintf("bfl_yz_lcva_%s_s%d.RData", test_site, seed_run)

if (file.exists(bfl_yz_cache)) {
  message("  Loading cached LCVA fits from: ", bfl_yz_cache)
  load(bfl_yz_cache)   # restores: src_lcva_raw, self_raw
} else {
  t0 <- tic("BFL_yz — LCVA source sites")
  src_lcva_raw <- setNames(lapply(source_sites, function(s) {
    fit_lcva(X_train = as.matrix(phmrc[phmrc$site==s, 1:168]),
             Y_train = phmrc[phmrc$site==s, "cause"], lcva_args = lcva_args)
  }), source_sites)
  timing_log[["BFL_yz_lcva_src"]] <- toc(t0, "BFL_yz — LCVA source sites")

  t0 <- tic("BFL_yz — LCVA self-site")
  self_raw <- fit_lcva(X_train = X_labeled, Y_train = Y_labeled, lcva_args = lcva_args)
  timing_log[["BFL_yz_lcva_self"]] <- toc(t0, "BFL_yz — LCVA self-site")

  save(src_lcva_raw, self_raw, file = bfl_yz_cache)
  message("  Saved LCVA cache to: ", bfl_yz_cache)
}

src_fits      <- lapply(src_lcva_raw, `[[`, "fit")
src_cause_ids <- lapply(src_lcva_raw, `[[`, "cause_ids")
self_fits     <- setNames(list(self_raw$fit),       test_site)
self_cids     <- setNames(list(self_raw$cause_ids), test_site)

# ============================================================
# PART 2b — BFL_yz Stan runs  ← re-run from here after cache exists
# ============================================================
message("\n=== BFL_yz — Stan runs (seed ", seed_run, ") ===")

res_yz <- tryCatch({
  t_part2 <- tic("BFL_yz total")

  dom_fits  <- c(src_fits,      self_fits)
  dom_cids  <- c(src_cause_ids, self_cids)
  burn_in   <- 200L
  sa_yz     <- modifyList(sa, list(init = "0"))   # stable init for large M

  # BFL_yz Domain
  t0 <- tic("BFL_yz Domain")
  fit_yz_dom <- run_BFL_yz(dom_fits, dom_cids, X_full, Y_partial,
                            burn_in = burn_in, label_shift = FALSE, stan_args = sa_yz)
  timing_log[["BFL_yz_Domain"]] <- toc(t0, "BFL_yz Domain")
  sc_yz_dom  <- score_BFL(fit_yz_dom, Y_full, eval_idx = idx_unlabeled)

  # BFL_yz Partial
  t0 <- tic("BFL_yz Partial")
  fit_yz_par <- run_BFL_yz(src_fits, src_cause_ids, X_full, Y_partial,
                            burn_in = burn_in, label_shift = FALSE, stan_args = sa_yz)
  timing_log[["BFL_yz_Partial"]] <- toc(t0, "BFL_yz Partial")
  sc_yz_par  <- score_BFL(fit_yz_par, Y_full, eval_idx = idx_unlabeled)

  # BFL_yz Mix
  t0 <- tic("BFL_yz Mix")
  fit_yz_mix <- run_BFL_yz(dom_fits, dom_cids, X_full, Y_partial,
                            burn_in = burn_in, label_shift = FALSE, stan_args = sa_yz)
  timing_log[["BFL_yz_Mix"]] <- toc(t0, "BFL_yz Mix")
  sc_yz_mix  <- score_BFL(fit_yz_mix, Y_full, eval_idx = idx_unlabeled)

  timing_log[["BFL_yz_total"]] <<- toc(t_part2, "BFL_yz total")

  data.frame(
    site = test_site, seed = seed_run, type = "BFL_yz",
    method = c("Domain", "Partial", "Mix"),
    csmf = c(sc_yz_dom$csmf_acc, sc_yz_par$csmf_acc, sc_yz_mix$csmf_acc),
    top1 = c(sc_yz_dom$top1_acc, sc_yz_par$top1_acc, sc_yz_mix$top1_acc)
  )

}, error = function(e) {
  message("\nBFL_yz skipped — ", conditionMessage(e))
  message("Check: str(src_fits[[1]], max.level=1) for phi/loglambda field names.")
  NULL
})

# Results + timing
res_combined <- bind_rows(res_std, res_yz)
print(res_combined)

timing_tbl <- data.frame(step = names(timing_log),
                          secs = round(unlist(timing_log), 1), row.names = NULL)
cat("\n--- Timing summary ---\n"); print(timing_tbl)

# Plot
res_long <- res_combined %>%
  pivot_longer(c(csmf, top1), names_to = "metric", values_to = "value") %>%
  mutate(metric = recode(metric, csmf = "CSMF accuracy", top1 = "Top-1 accuracy"),
         method = factor(method, c("Domain", "Partial", "Mix")),
         type   = factor(type,   c("BFL2", "BFL_yz")))

ggplot(res_long, aes(method, value, colour = type, shape = type)) +
  geom_point(size = 4) +
  scale_colour_manual(values = c(BFL2 = "grey30", BFL_yz = "steelblue"),
                      labels = c("BFL2 P(X|Y)", "BFL_yz P(X|Y,Z)")) +
  scale_shape_manual(values  = c(BFL2 = 16, BFL_yz = 18),
                     labels  = c("BFL2 P(X|Y)", "BFL_yz P(X|Y,Z)")) +
  facet_wrap(~metric, scales = "free_y", nrow = 2) +
  labs(title    = paste("BFL2 quick test —", test_site, "| seed =", seed_run),
       subtitle = "● BFL2 (P(X|Y))   ◆ BFL_yz (P(X|Y,Z))",
       x = NULL, y = NULL, colour = NULL, shape = NULL) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1),
        legend.position = "bottom")
