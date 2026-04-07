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

**Dependencies:** `rstan`, `ggplot2`, `dplyr`, `tibble`, `scales`, `rlang`

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
│   ├── BFL_methods.R                # S3: print/summary/plot + site_similarity_BFL
│   ├── run_bfl_stan.R               # Stan runner (wraps rstan::sampling)
│   ├── build_bfl_stan_data.R        # Stan data construction (3 variants)
│   ├── align_causes.R               # Cause alignment across sites
│   ├── row_alignment.R              # Row hash matching utilities
│   ├── hash_rows.R                  # compute_row_hashes()
│   ├── label_mask.R                 # make_label_mask()
│   ├── label_mask_split_helpers.R   # split_idx_balanced / shift_beta / shift_dirichlet
│   ├── LCVA_helpers.R               # LCVA interface layer (fit_lcva / predict_lcva / run_lcva)
│   ├── CSMF_acc.R                   # CSMF_acc() metric
│   ├── get_ACC.R                    # get_ACC() accuracy helper
│   ├── bfl_post_pred_sampler.R      # Internal posterior predictive draws
│   └── validate_inputs.R            # Input validation helpers
├── inst/
│   └── stan/
│       ├── no_partial_labels.stan           # Base / Domain model
│       ├── partial_labels_shift_false.stan  # Partial / Mix (no label shift)
│       └── partial_labels_shift_true.stan   # Partial / Mix (label shift)
├── vignettes/
│   └── toy_example.Rmd
├── scripts/
│   └── reproducible_example.R
└── tests/
    └── testthat/
        └── test-package_assets.R
```

---

## BFL variants

`run_BFL()` automatically selects the correct Stan model based on how you build `local_summaries` and what you pass as `Y_target`. No model-type flag is needed.

| Variant | Labeled target data? | `local_summaries` built on | `Y_target` | Stan model |
|---------|----------------------|----------------------------|------------|------------|
| **Base** | No | All N rows | `NULL` | `no_partial_labels` |
| **Domain** | Yes (all as training) | Unlabeled rows only | Full N (NA for unlabeled) | `no_partial_labels` |
| **Partial** | Yes (some as signal to Stan) | All N rows | Full N (NA for unlabeled) | `partial_labels_shift_false` |
| **Mix** | Yes (domain 10% + partial 10%) | All N rows + self-site | Full N (NA for unlabeled) | `partial_labels_shift_false` |

Set `label_shift = TRUE` in `run_BFL()` to use `partial_labels_shift_true` (unbalanced/shift variant) for Partial or Mix.

---

## Core workflow

```r
library(BFL)

# 1. Fit a local LCVA model per training site, get per-record scores on target
#    (phi: N x C matrix of P(X | cause) probabilities)
local_fit <- run_lcva(X_train, Y_train, ...)
phi       <- predict_lcva(local_fit, X_target, ...)

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
  local_summaries = local_summaries,
  X_target        = X_target,
  Y_target        = Y_target,     # NA for unlabeled records; NULL for Base
  stan_args       = list(iter = 4000, chains = 4, seed = 42)
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

# 6. Site diagnostics
site_similarity_BFL(fit, plot = TRUE)  # pairwise phi correlation heatmap
```

---

## Key functions

| Function | Description |
|----------|-------------|
| `run_BFL()` | Main fitting function. Returns a `"BFL"` object. |
| `score_BFL()` | Scores a `"BFL"` fit against ground truth. Returns a `"BFL_score"` object. |
| `compute_row_hashes()` | Produces a character hash per row for cross-site row alignment. |
| `run_lcva()` | Fits a single-domain LCVA model (wraps `LCVA::LCVA.train`). |
| `predict_lcva()` | Gets per-record posterior phi from a fitted LCVA. |
| `site_similarity_BFL()` | Pairwise Pearson correlation between local site phi matrices. |
| `split_idx_balanced()` | Creates a random balanced (equal-proportion) labeled/unlabeled split. |
| `split_idx_shift_dirichlet()` | Creates a mild-shift split (Dirichlet-drawn prevalences). |
| `split_idx_shift_beta()` | Creates a severe-shift split (Beta(0.2, 0.2) per-cause fractions). |

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
| `stan_idx` | integer | Rows that entered Stan |
| `n_total` | integer | Total N (= `nrow(X_target)`) |
| `nLc` | named int or NULL | Labeled counts outside Stan (for CSMF correction) |
| `has_labels` | logical | Whether partial labels were passed to Stan |
| `label_shift` | logical | Whether the unbalanced shift variant was used |

---

## `score_BFL()` notes

- `Y_eval` must have `length(Y_eval) == fit$n_total`. For Domain/Mix, slice to `test_idx` before passing.
- `eval_idx` is expressed in `Y_eval` (= `test_idx`) space. For Mix, use `match(missing_indices, test_idx)`.
- CSMF correction is applied automatically when `n_stan < n_total` using `fit$nLc`.

---

## Stan models

| File | Variant | When used |
|------|---------|-----------|
| `no_partial_labels.stan` | Base, Domain | `Y_target` is NULL, or all labeled rows fall outside Stan |
| `partial_labels_shift_false.stan` | Partial, Mix | Some Stan rows have known labels; `label_shift = FALSE` |
| `partial_labels_shift_true.stan` | Partial (shift), Mix (shift) | Same but `label_shift = TRUE` |

---

## Notes

- The package was renamed from `BFL2` → `BFL`. Install with `library(BFL)`.
- LCVA fitting requires the `LCVA` package to be installed separately.
- For simulation replication scripts, see the `BFL2_Rerun_Simulation/` project directory.
