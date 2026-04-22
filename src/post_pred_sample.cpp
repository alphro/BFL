#include <Rcpp.h>
#include <vector>
#include <algorithm>
using namespace Rcpp;

// ---------------------------------------------------------------------------
// Posterior predictive sampler — replaces bfl_post_pred_sampler()
//
// For each posterior draw s = 1..S and observation i = 1..I:
//   score[c] = pi[s,c] * sum_m phi[i,c,m] * lambda[s,c,m]
//   Normalize score to sum to 1.
//   Sample Y[s,i] ~ Categorical(score).
//   Accumulate score into prob_mean.
//
// All arrays are passed as flat vectors in R's column-major order:
//   phi[i,c,m]    -> phi_vec[i + c*I + m*I*C]
//   lambda[s,c,m] -> lambda_vec[s + c*S + m*S*C]
//   pi[s,c]       -> pi_vec[s + c*S]
// ---------------------------------------------------------------------------

// Sample 0-based index from a normalized score array (sums to 1).
static int rcategorical_norm(const double* p, int n) {
  double u = R::unif_rand();
  double cum = 0.0;
  for (int c = 0; c < n - 1; c++) {
    cum += p[c];
    if (u < cum) return c;
  }
  return n - 1;
}

//' Posterior predictive sampler (Rcpp)
//'
//' @param phi_vec    Flattened phi array [I, C, M] column-major.
//' @param lambda_vec Flattened lambda array [S, C, M] column-major.
//' @param pi_vec     Flattened pi matrix [S, C] column-major.
//' @param I,C,M,S   Dimensions.
//'
//' @return List: posterior_pred_Y_prob_mean [I x C], posterior_pred_Y [S x I].
//' @keywords internal
// [[Rcpp::export]]
List post_pred_sample_cpp(
    const NumericVector& phi_vec,
    const NumericVector& lambda_vec,
    const NumericVector& pi_vec,
    int I, int C, int M, int S
) {
  NumericMatrix prob_mean(I, C);  // running sum -> divided by S at end
  IntegerMatrix draws_int(S, I);  // 1-indexed class draws

  std::vector<double> score(C);

  for (int s = 0; s < S; s++) {
    for (int i = 0; i < I; i++) {

      // score[c] = pi[s,c] * sum_m phi[i,c,m] * lambda[s,c,m]
      double total = 0.0;
      for (int c = 0; c < C; c++) {
        double sc = 0.0;
        for (int m = 0; m < M; m++)
          sc += phi_vec[i + c*I + (long)m*I*C] * lambda_vec[s + c*S + (long)m*S*C];
        score[c] = pi_vec[s + c*S] * sc;
        total   += score[c];
      }
      if (total <= 0.0) total = 1.0;

      // Normalize, accumulate, sample
      for (int c = 0; c < C; c++) {
        score[c]       /= total;
        prob_mean(i, c) += score[c];
      }
      draws_int(s, i) = rcategorical_norm(score.data(), C) + 1;  // 1-indexed
    }
  }

  for (int i = 0; i < I; i++)
    for (int c = 0; c < C; c++)
      prob_mean(i, c) /= S;

  return List::create(
    Named("posterior_pred_Y_prob_mean") = prob_mean,
    Named("posterior_pred_Y")           = draws_int
  );
}

// ---------------------------------------------------------------------------
// Modal class per column of an integer draw matrix
//
// Replaces:
//   apply(draws_eval, 2, function(x) which.max(tabulate(x, nbins = C)))
//
// draws: [S x I] integer matrix, values in 1..C.
// Returns: integer vector of length I (1-indexed modal class per column).
// ---------------------------------------------------------------------------

//' Modal class vote across posterior draws (Rcpp)
//'
//' @param draws Integer matrix [S x I] with values in 1..C.
//' @param C     Number of classes.
//' @return Integer vector of length I: 1-indexed modal class per observation.
//' @keywords internal
// [[Rcpp::export]]
IntegerVector modal_class_cpp(const IntegerMatrix& draws, int C) {
  int S = draws.nrow();
  int I = draws.ncol();
  IntegerVector result(I);
  std::vector<int> counts(C);

  for (int i = 0; i < I; i++) {
    std::fill(counts.begin(), counts.end(), 0);
    for (int s = 0; s < S; s++) {
      int cls = draws(s, i) - 1;  // 0-indexed
      if (cls >= 0 && cls < C) counts[cls]++;
    }
    result[i] = (int)(std::max_element(counts.begin(), counts.end())
                      - counts.begin()) + 1;  // 1-indexed
  }
  return result;
}
