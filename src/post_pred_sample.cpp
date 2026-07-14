#include <Rcpp.h>
#include <vector>
#include <algorithm>
using namespace Rcpp;

// ---------------------------------------------------------------------------
// Draw a categorical sample from normalized probabilities.
// probs must sum to 1.
// Returns 0-indexed class.
// ---------------------------------------------------------------------------
inline int rcategorical_norm(const double* probs, int C) {
  double u = R::runif(0.0, 1.0);
  double cum = 0.0;

  for (int c = 0; c < C; c++) {
    cum += probs[c];
    if (u <= cum) return c;
  }

  return C - 1;
}

// ---------------------------------------------------------------------------
// Posterior predictive sampler — Zoey-compatible stochastic prediction
//
// For each posterior draw s and observation i:
//   score[c] = pi[s,c] * sum_m phi[i,c,m] * lambda[s,c,m]
//
// score is normalized to sum to 1.
//
// Then:
//   Y[s,i] ~ Categorical(score)
//
// This reproduces Zoey's:
//
//   sample(1:C, 1, prob = x)
//
// behavior exactly.
//
// R column-major layout:
//   phi[i,c,m]    -> phi_vec[i + c*I + m*I*C]
//   lambda[s,c,m] -> lambda_vec[s + c*S + m*S*C]
//   pi[s,c]       -> pi_vec[s + c*S]
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
List post_pred_sample_cpp(
    const NumericVector& phi_vec,
    const NumericVector& lambda_vec,
    const NumericVector& pi_vec,
    int I, int C, int M, int S
) {
  NumericMatrix prob_mean(I, C);
  IntegerMatrix draws_int(S, I);

  std::vector<double> score(C);

  for (int s = 0; s < S; s++) {

    for (int i = 0; i < I; i++) {

      double total = 0.0;

      // ---------------------------------------------------------------
      // Compute unnormalized scores
      // ---------------------------------------------------------------
      for (int c = 0; c < C; c++) {

        double sc = 0.0;

        for (int m = 0; m < M; m++) {

          sc +=
            phi_vec[i + c * I + (long)m * I * C] *
            lambda_vec[s + c * S + (long)m * S * C];
        }

        score[c] = pi_vec[s + c * S] * sc;
        total += score[c];
      }

      // ---------------------------------------------------------------
      // Normalize
      // ---------------------------------------------------------------
      if (total <= 0.0 || !R_finite(total)) {

        for (int c = 0; c < C; c++) {
          score[c] = 1.0 / C;
        }

      } else {

        for (int c = 0; c < C; c++) {
          score[c] /= total;
        }
      }

      // ---------------------------------------------------------------
      // Accumulate posterior mean probabilities
      // ---------------------------------------------------------------
      for (int c = 0; c < C; c++) {
        prob_mean(i, c) += score[c];
      }

      // ---------------------------------------------------------------
      // Zoey-compatible stochastic categorical draw
      // ---------------------------------------------------------------
      draws_int(s, i) = rcategorical_norm(score.data(), C) + 1;
    }
  }

  // -------------------------------------------------------------------
  // Posterior mean probabilities
  // -------------------------------------------------------------------
  for (int i = 0; i < I; i++) {
    for (int c = 0; c < C; c++) {
      prob_mean(i, c) /= S;
    }
  }

  return List::create(
    Named("posterior_pred_Y_prob_mean") = prob_mean,
    Named("posterior_pred_Y") = draws_int
  );
}

// ---------------------------------------------------------------------------
// Modal class per column of an integer draw matrix
// draws: [S x I], values in 1..C
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
IntegerVector modal_class_cpp(const IntegerMatrix& draws, int C) {

  int S = draws.nrow();
  int I = draws.ncol();

  IntegerVector result(I);

  std::vector<int> counts(C);

  for (int i = 0; i < I; i++) {

    std::fill(counts.begin(), counts.end(), 0);

    for (int s = 0; s < S; s++) {

      int cls = draws(s, i) - 1;

      if (cls >= 0 && cls < C) {
        counts[cls]++;
      }
    }

    result[i] =
      static_cast<int>(
        std::max_element(counts.begin(), counts.end()) - counts.begin()
      ) + 1;
  }

  return result;
}