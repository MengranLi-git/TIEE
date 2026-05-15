library(Rcpp)

#' Coverage utility
coverage <- function(tq, upper, lower) {
  sum(lower <= tq & upper >= tq, na.rm = TRUE)
}

# --- Rcpp: GPD score function (n x p matrix) ---
cppFunction('
NumericMatrix gpd_cpp(
    const NumericVector& param,
    const NumericVector& exc,
    const NumericMatrix& X,
    const IntegerVector& scale_cols,
    const IntegerVector& shape_cols
){
    int n  = exc.size();
    int p1 = scale_cols.size();
    int p2 = shape_cols.size();
    int p  = p1 + p2;

    // Split parameters: beta_sigma, beta_xi
    NumericVector beta_sigma(p1), beta_xi(p2);
    for(int k=0; k<p1; k++) beta_sigma[k] = param[k];
    for(int k=0; k<p2; k++) beta_xi[k]    = param[p1 + k];

    NumericMatrix score(n, p);

    for (int i = 0; i < n; i++){
        double y = exc[i];

        // Non-exceedance: score row = 0
        if (y <= 0.0) {
            for (int j = 0; j < p; j++) score(i, j) = 0.0;
            continue;
        }

        // Extract covariate rows
        NumericVector Xs(p1), Xx(p2);
        for(int k=0; k<p1; k++) Xs[k] = X(i, scale_cols[k]);
        for(int k=0; k<p2; k++) Xx[k] = X(i, shape_cols[k]);

        // Compute sigma (log-link) and xi
        double eta_s = 0.0, eta_x = 0.0;
        for(int k=0; k<p1; k++) eta_s += Xs[k] * beta_sigma[k];
        for(int k=0; k<p2; k++) eta_x += Xx[k] * beta_xi[k];
        double sg = std::exp(eta_s);
        double xv = eta_x;
        double z  = 1.0 + xv * (y / sg);

        // Partial derivatives of log-likelihood w.r.t. sigma, xi
        double dls = 0.0, dlx = 0.0;
        if (z <= 0) {
            dls = 0.0; dlx = 0.0;
        } else if (std::abs(xv) < 1e-8) {
            // Exponential limit (xi ~ 0)
            dls = -1.0/sg + y/(sg*sg);
            dlx = 0.0;
        } else {
            dls = -1.0/sg + (1.0 + xv) * (y/(sg*sg)) / z;
            dlx = std::log(z)/(xv*xv) - (1.0 + 1.0/xv) * (y/sg) / z;
        }

        // Chain rule: score w.r.t. beta_sigma
        for (int k = 0; k < p1; k++) score(i, k) = dls * sg * Xs[k];
        // Chain rule: score w.r.t. beta_xi
        for (int k = 0; k < p2; k++) score(i, p1 + k) = dlx * Xx[k];
    }
    return score;
}
')

# --- Rcpp: Jacobian of GPD quantile function ---
cppFunction('
NumericMatrix G_gpd_cpp(NumericVector tau,
                        NumericVector diffs,
                        NumericVector xi,
                        NumericVector sigma,
                        double u) {
  int n = xi.size();
  int T = tau.size();
  NumericMatrix out(n, 2);

  for(int i = 0; i < n; i++){
    double xi_i    = xi[i];
    double sigma_i = sigma[i];
    double sum_dxi = 0.0, sum_dsig = 0.0;

    for(int j = 0; j < T; j++){
      double t = tau[j], d = diffs[j];
      double ratio = (1.0 - u) / (1.0 - t);
      double A = pow(ratio, xi_i);
      double dq_dsigma = (A - 1.0) / xi_i;
      double dq_dxi = -sigma_i/(xi_i*xi_i) * (A - 1.0) +
                       (sigma_i/xi_i) * A * log(ratio);
      sum_dsig += d * dq_dsigma;
      sum_dxi  += d * dq_dxi;
    }
    out(i, 0) = sum_dxi;
    out(i, 1) = sum_dsig;
  }
  return out;
}
')

# --- Rcpp: Vectorised GPD density ---
cppFunction('
NumericVector dgpd_cpp(
    const NumericVector& exc,
    const NumericVector& scale,
    const NumericVector& shape
){
    int n = exc.size();
    NumericVector out(n);
    for(int i = 0; i < n; i++){
        double y  = exc[i], sg = scale[i], xi = shape[i];
        if(y <= 0.0 || sg <= 0.0 || ISNAN(y) || ISNAN(sg) || ISNAN(xi)){
            out[i] = 0.0; continue;
        }
        double z = 1.0 + xi * (y / sg);
        if(z <= 0.0){ out[i] = 0.0; continue; }
        if(std::abs(xi) < 1e-8){
            out[i] = (1.0/sg) * std::exp(-y/sg);
        } else {
            out[i] = (1.0/sg) * std::pow(z, -1.0/xi - 1.0);
        }
    }
    return out;
}
')

# --- Rcpp: Vectorised GPD quantile function (n x T matrix) ---
cppFunction('
NumericMatrix qgpd_matrix(NumericVector tau,
                          NumericVector sigma,
                          NumericVector xi) {
  int n = sigma.size(), T = tau.size();
  NumericMatrix Q(n, T);
  for (int i = 0; i < n; i++) {
    double sig = sigma[i], xii = xi[i], sig_over_xi = sig / xii;
    for (int j = 0; j < T; j++) {
      double t = tau[j];
      Q(i, j) = sig_over_xi * (pow(1.0 - t, -xii) - 1.0);
    }
  }
  return Q;
}
')

# --- Rcpp: Objective function for theta estimation ---
cppFunction('
double obj_cpp(double theta,
               const NumericMatrix& Yd,
               const NumericVector& w,
               const NumericVector& diffs,
               double eff,
               double R) {
  int n = Yd.nrow(), T = Yd.ncol();
  double d1 = diffs[0];
  double sum_w = 0.0, sum_diffs = 0.0;
  for (int i = 0; i < w.size(); i++) sum_w += w[i];
  for (int i = 0; i < diffs.size(); i++) sum_diffs += diffs[i];

  double sum1 = 0.0;
  for (int t = 0; t < T; t++)
    for (int i = 0; i < n; i++)
      sum1 += d1 * w[i] * std::abs(Yd(i,t) - theta);

  double sum2 = std::abs(R + theta * sum_w * sum_diffs);
  double sum3 = std::abs(R - theta * (2.0 * n * eff) * sum_diffs);
  return sum1 + sum2 + sum3;
}
')

# --- Rcpp: KDE path via golden-section search ---
cppFunction('
#include <Rcpp.h>
using namespace Rcpp;

inline double obj_fast_inline(double theta,
                              const NumericMatrix& Yd,
                              const NumericVector& w,
                              double d1, double sum_w,
                              double sum_diffs, double eff, double R) {
  int n = Yd.nrow(), T = Yd.ncol();
  double sum1 = 0.0;
  for (int i = 0; i < n; i++) {
    double wi = w[i] * d1;
    for (int t = 0; t < T; t++) sum1 += wi * std::abs(Yd(i, t) - theta);
  }
  double sum2 = std::fabs(R + theta * sum_w * sum_diffs);
  double sum3 = std::fabs(R - theta * (2.0 * n * eff) * sum_diffs);
  return sum1 + sum2 + sum3;
}

double golden_min(const NumericMatrix& Yd, const NumericVector& w,
                  double d1, double sum_w, double sum_diffs,
                  double eff, double R,
                  double a, double b, double tol, int max_iter) {
  const double gr = (std::sqrt(5.0) - 1.0) / 2.0;
  double c = b - gr * (b - a), d = a + gr * (b - a);
  double fc = obj_fast_inline(c, Yd, w, d1, sum_w, sum_diffs, eff, R);
  double fd = obj_fast_inline(d, Yd, w, d1, sum_w, sum_diffs, eff, R);
  int iter = 0;
  while ((b - a) > tol && iter < max_iter) {
    if (fc < fd) {
      b = d; d = c; fd = fc;
      c = b - gr * (b - a);
      fc = obj_fast_inline(c, Yd, w, d1, sum_w, sum_diffs, eff, R);
    } else {
      a = c; c = d; fc = fd;
      d = a + gr * (b - a);
      fd = obj_fast_inline(d, Yd, w, d1, sum_w, sum_diffs, eff, R);
    }
    iter++;
  }
  return 0.5 * (a + b);
}

// [[Rcpp::export]]
NumericVector kde_path_cpp(const NumericVector& tau2,
                           const NumericVector& interval,
                           const NumericMatrix& Yd,
                           const NumericVector& w,
                           double d1, double sum_w,
                           double sum_diffs, double R,
                           double delta, double tol,
                           int max_iter,
                           bool fallback_global = true) {
  double a0 = interval[0], b0 = interval[1];
  int m = tau2.size();
  NumericVector out(m);
  out[0] = golden_min(Yd, w, d1, sum_w, sum_diffs, tau2[0], R, a0, b0, tol, max_iter);
  for (int i = 1; i < m; i++) {
    double center = out[i-1];
    double a = std::max(a0, center - delta);
    double b = std::min(b0, center + delta);
    double th = golden_min(Yd, w, d1, sum_w, sum_diffs, tau2[i], R, a, b, tol, max_iter);
    if (fallback_global) {
      if (std::fabs(th - a) < 1e-12 || std::fabs(th - b) < 1e-12)
        th = golden_min(Yd, w, d1, sum_w, sum_diffs, tau2[i], R, a0, b0, tol, max_iter);
    }
    out[i] = th;
  }
  return out;
}
')
