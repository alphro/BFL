# BFL — Bayesian Federated Learning for Verbal Autopsy

**Version 0.1.0** · R package

BFL implements Bayesian Federated Learning methods for verbal autopsy (VA) cause-of-death assignment. It aggregates cause-of-death probability scores from multiple training sites into a single posterior estimate of cause-specific mortality fractions (CSMF) and per-record top-1 predictions, without sharing raw data across sites.

---

## Installation

```r
# From local source (devtools)
devtools::install_local("/path/to/BFLpkg")

# Or after cloning
install.packages(".", repos = NULL, type = "source")
```

**Dependencies:** `rstan`, `Rcpp`, `ggplot2`, `dplyr`, `tibble`, `scales`, `rlang`

---

## Package structure

```
BFLpkg/
├── DESCRIPTION
├── NAMESPACE
├── R/
│   ├── run_BFL.R                    # Main entry point
│   ├── score_BFL.R                  # Scoring against ground truth
│   ├── predict_BFL.R                # Posterior predictive sampler
│   ├── BFL_methods.R                # S3: print/summary/plot for BFL objects
│   ├── run_bfl_stan.R               # Stan runner (wraps rstan::sampling)
│   ├── build_bfl_stan_data.R        # Stan data construction (3 variants)
│   ├── align_causes.R               # Cause alignment across sites
│   ├── row_alignment.R              # Row hash matching utilities
│   ├── hash_rows.R                  # compute_row_hashes()
│   ├── CSMF_acc.R                   # CSMF_acc() metric
│   ├── bfl_post_pred_sampler.R      # Internal posterior predictive draws
│   └── validate_inputs.R            # Input validation helpers
├── inst/
│   └── stan/
│       ├── no_partial_labels.stan           # Base / Domain model
│       ├── partial_labels_shift_false.stan  # Partial / Mix (no label shift)
│       └── partial_labels_shift_true.stan   # Partial / Mix (label shift)
├── scripts/
│   └── reproducible_example.R
└── tests/
    └── testthat/
        └── test-package_assets.R
```

---

## BFL variants (input-driven)

`run_BFL()` is decided entirely by **what you pass in** — there is no model-type flag and no row inference. The rows of `local_summaries`, `X_target` and `Y_target` must line up **1:1**; you assemble any labeled rows (and the self-model column) yourself before calling. Base/Domain/Partial/Mix are paper names, not code paths.

| Paper name | `X_target` / phi rows | `Y_target` | Self-model column in `local_summaries`? | `Y_add` | `label_shift` | Stan model |
|---|---|---|---|---|---|---|
| **Base** | N | `NULL` | no | `NULL` | — | `no_partial_labels` |
| **Domain** | N | `NULL` | yes (Tgt New) | held labels | — | `no_partial_labels` |
| **Partial** | N+L (NA + labels) | length N+L | no | `NULL` | T/F | `partial_labels_shift_{false,true}` |
| **Mix** | N+P (NA + labels) | length N+P | yes (Tgt New) | held labels | T/F | `partial_labels_shift_{false,true}` |

The package only checks three things: is `Y_target` `NULL`? is `Y_add` supplied? is `label_shift` set? `Y_add` carries held-out labels (for rows **not** in `X_target`) used solely for the CSMF correction `(n_total·π + nLc)/(n_total + n_add)`. Set `label_shift = TRUE` to use the unbalanced/shift Stan variant for Partial/Mix.

---

## Samplers

`run_BFL()` aggregates with either backend via the `sampler` argument:

- **`sampler = "gibbs"`** (default) — fast conjugate **Rcpp Gibbs** sampler. `gibbs_args = list(logistic_normal = FALSE)` uses a Dirichlet prior (`gibbs_dir`); `list(logistic_normal = TRUE, mh_scale = 0.25)` matches Stan's logistic-normal prior via a Metropolis step (`gibbs_ln`).
- **`sampler = "stan"`** — the original Stan/NUTS path.

Both give the same `pi`/`lambda` up to MCMC noise; Gibbs is roughly **10–15× faster**. Shared MCMC controls go in `mcmc_args = list(iter, chains, seed)`. (Partial-label variants always run through Stan-equivalent logic regardless of `sampler`.)

---

## Core workflow

```r
library(BFL)

# 1. Fit a local model per training site (e.g. LCVA) and get per-record scores
#    on the target (phi: N x C matrix of P(X | cause) probabilities)
local_fit <- LCVA::LCVA.train(X_train, Y_train, ...)
phi       <- LCVA::LCVA.pred(local_fit, X_target, ...)$phi

# 2. Build a local_summaries entry for each site
row_hash <- compute_row_hashes(X_target)

local_summaries <- list(
  site_A = list(
    posterior_phi = phi_A,         # N x C matrix
    cause_ids     = cause_ids_A,   # character vector of length C
    target_info   = list(row_hash = row_hash, N = nrow(phi_A), P = ncol(X_target))
  ),
  site_B = list(...)
)

# 3. Run BFL global aggregation
fit <- run_BFL(
  local_summaries = local_summaries,  # rows match X_target 1:1
  X_target        = X_target,
  Y_target        = Y_target,     # NA for unlabeled records; NULL for Base/Domain
  Y_add           = Y_add,        # held labels for the CSMF correction; NULL otherwise
  sampler         = "gibbs",      # "gibbs" (default, fast) or "stan"
  mcmc_args       = list(iter = 2000, chains = 4, seed = 42),
  gibbs_args      = list(logistic_normal = FALSE)   # ignored when sampler = "stan"
)

# 4. Inspect results
print(fit)          # compact overview
summary(fit)        # full CSMF table + lambda weights per site
plot(fit)           # CSMF dot-plot with credible intervals

# 5. Score against ground truth (unlabeled rows only)
score <- score_BFL(
  fit      = fit,
  Y_eval   = Y_eval,       # ground-truth labels, length == fit$n_total
  eval_idx = unlabeled_idx # positions of unlabeled rows (in Y_eval space)
)
print(score)    # top-1, balanced accuracy, CSMF accuracy
summary(score)  # per-cause recall + CSMF error table
plot(score)     # list of two ggplots: recall bars + CSMF error bars
```

---

## Key functions

| Function | Description |
|----------|-------------|
| `run_BFL()` | Main fitting function. Returns a `"BFL"` object. |
| `score_BFL()` | Scores a `"BFL"` fit against ground truth. Returns a `"BFL_score"` object. |
| `predict_BFL()` | Posterior predictive draws from a fitted `"BFL"` model. |
| `compute_row_hashes()` | Produces a character hash per row for cross-site row alignment. |
| `CSMF_acc()` | CSMF accuracy between predicted and true cause prevalence. |

---

## `run_BFL()` return value

An object of class `"BFL"`:

| Field | Type | Description |
|-------|------|-------------|
| `pi` | S × C matrix | Posterior draws of CSMF |
| `lambda` | S × C × M array | Posterior draws of per-cause site weights |
| `phi` | N × C × M array | Aligned local phi scores |
| `causes` | character[C] | Global cause labels |
| `model_names` | character[M] | Site/model names |
| `row_hash` | character[N] | Reference row hashes for `X_target` |
| `stan_idx` | integer | Rows that entered the model (always `1:nrow(X_target)`) |
| `n_total` | integer | Model rows (= `nrow(X_target)`) |
| `n_add` | integer | Held rows behind `nLc` (0 if none) |
| `nLc` | named int or NULL | Cause counts from `Y_add` (for CSMF correction) |
| `has_labels` | logical | Whether partial labels were passed to Stan |
| `label_shift` | logical | Whether the unbalanced shift variant was used |

---

## `score_BFL()` notes

- `Y_eval` must have `length(Y_eval) == fit$n_total` (the model rows). Held `Y_add` rows are **not** in `Y_eval`; they enter only via `fit$nLc`.
- `eval_idx` is in `1..n_total` space; since `stan_idx == 1:n_total`, it maps to draw columns directly.
- CSMF correction is applied automatically when `Y_add` was supplied: `(n_total·π + nLc)/(n_total + n_add)`.

---

## Stan models

| File | Variant | When used |
|------|---------|-----------|
| `no_partial_labels.stan` | Base, Domain | `Y_target` is `NULL` (or all `NA`) |
| `partial_labels_shift_false.stan` | Partial, Mix | `Y_target` has labels; `label_shift = FALSE` |
| `partial_labels_shift_true.stan` | Partial (shift), Mix (shift) | Same but `label_shift = TRUE` |

---

## Developer Log

> **BFL is born! 🍼 → 📦 — baby steps to a big package.** A version-by-version story, newest first. (Minor versions = steady progress; a new *major* version = a big redesign.)

**v2.0 — in progress · "Spring cleaning" 🧹**
- `run_BFL()` is now input-driven: rows are 1:1 (`local_summaries` = `X_target` = `Y_target`), and the variant is decided by `Y_target`/`Y_add` — no more `n_stan` vs `n_total` inference.
- New `Y_add` argument carries held-out labels for the CSMF correction `(n_total·π + nLc)/(n_total + n_add)` (previously inferred from rows outside Stan).
- Removed the hash-based `stan_idx` matcher and the Domain `is.na(Y_target)` path; `stan_idx` is always `1:N`.

**v1.1 — Jun 2026 · "Cracking the hard case" 💪**
- Fixed severe label-shift (cause-ID alignment, so weights spread across all causes).
- Per-batch CHAMPS + PHMRC result pipelines with cluster orchestration.
- Figures 08/09/10 rebuilt to run from saved results.

**v1.0 — Apr–May 2026 · "It's a real package now!" 🚀**
- C++ Gibbs sampler alongside Stan — ~10× faster, same results.
- Added OpenVA / InSilicoVA and pooled/joint LCVA as comparison methods.
- Large no/mild/severe-shift runs on Hummingbird; balanced-accuracy convention settled.

**v0.3 — Mar 2026 · "Auto-pilot" 🤖**
- `run_BFL()` auto-detects the variant (Base/Domain/Partial/Mix) — no model flag.
- S3 `print`/`summary`/`plot` methods; cleaner API.

**v0.2 — Jan–Feb 2026 · "Finding its shape" 🧩**
- Three-step workflow: local LCVA → global BFL → predict.
- Row hashing + automatic alignment of inputs to `X_target`.
- Base/Domain validated on all 6 sites; Partial/Mix debugged.

**v0.1 — Late 2025 · "BFL is born" 🎉**
- Read the BFL paper + reference prototype.
- Built the first minimal single-domain version.

---

## Notes

- The package was renamed from `BFL2` → `BFL`. Install with `library(BFL)`.
- LCVA fitting requires the `LCVA` package to be installed separately.
- For simulation replication scripts, see the `BFL2_Rerun_Simulation/` project directory.
