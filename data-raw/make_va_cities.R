# ==============================================================
# make_va_cities.R --- generative script for inst/extdata/va_cities.csv
#
# A small, self-contained synthetic verbal-autopsy (VA) dataset used by the
# BFL vignette. Four US cities, four causes of death, ten binary symptoms.
#
#   * Source sites (fully labeled): New York, Seattle, Austin
#   * Target site (labels partly hidden in the vignette): Los Angeles
#
# Deaths are drawn per city from a city-specific cause prevalence (the "twist"
# that makes source weighting meaningful), then symptoms are drawn from a
# cause -> symptom Bernoulli signature (shared across cities up to small
# per-city jitter). No real data; purely illustrative.
#
# Run from the package root to regenerate:
#   Rscript data-raw/make_va_cities.R
# The shipped inst/extdata/va_cities.csv is one frozen realization (seed 42-ish;
# R's RNG differs from the original build, so re-running overwrites with a new
# but equivalent draw --- fine for an illustrative dataset).
# ==============================================================
set.seed(42)

causes   <- c("COVID19", "LungCancer", "HeartAttack", "Stroke")
symptoms <- c("fever", "cough", "breathless", "chronic_cough", "weight_loss",
              "fatigue", "chest_pain", "numbness", "slurred_speech", "confusion")

# base P(symptom | cause): rows = cause, cols = symptom
base <- matrix(c(
  # fev  cou  bre  ccou wl   fat  cp   num  slur conf
  0.85, 0.80, 0.70, 0.15, 0.10, 0.55, 0.10, 0.05, 0.05, 0.10,  # COVID19
  0.15, 0.55, 0.55, 0.85, 0.80, 0.70, 0.20, 0.05, 0.05, 0.10,  # LungCancer
  0.10, 0.20, 0.75, 0.10, 0.10, 0.50, 0.85, 0.10, 0.10, 0.15,  # HeartAttack
  0.10, 0.10, 0.25, 0.05, 0.10, 0.40, 0.15, 0.85, 0.80, 0.75   # Stroke
), nrow = 4, byrow = TRUE, dimnames = list(causes, symptoms))

# per-city cause prevalence (the twist)
prev <- list(
  NewYork    = c(COVID19 = 0.20, LungCancer = 0.20, HeartAttack = 0.45, Stroke = 0.15),  # heart-heavy
  Seattle    = c(COVID19 = 0.45, LungCancer = 0.20, HeartAttack = 0.20, Stroke = 0.15),  # covid-heavy
  Austin     = c(COVID19 = 0.20, LungCancer = 0.45, HeartAttack = 0.20, Stroke = 0.15),  # cancer-heavy
  LosAngeles = c(COVID19 = 0.30, LungCancer = 0.30, HeartAttack = 0.25, Stroke = 0.15)   # TARGET: mixed
)
N <- c(NewYork = 120, Seattle = 120, Austin = 120, LosAngeles = 120)

rows <- list()
for (city in names(prev)) {
  # each city has its own theta = base signature + small jitter
  theta <- pmin(pmax(base + matrix(rnorm(length(base), 0, 0.05), nrow = 4), 0.03), 0.97)
  yv <- sample(causes, N[[city]], replace = TRUE, prob = prev[[city]][causes])
  X  <- t(vapply(yv, function(y) rbinom(length(symptoms), 1, theta[y, ]), numeric(length(symptoms))))
  rows[[city]] <- data.frame(site = city, cause = yv, X, stringsAsFactors = FALSE)
}
va_cities <- do.call(rbind, rows)
colnames(va_cities) <- c("site", "cause", symptoms)
rownames(va_cities) <- NULL

write.csv(va_cities, "data-raw/va_cities.csv", row.names = FALSE)
cat("wrote data-raw/va_cities.csv:", nrow(va_cities), "rows\n")
print(table(va_cities$site, va_cities$cause))
