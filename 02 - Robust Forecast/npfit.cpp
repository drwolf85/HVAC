/*** R
suppressMessages(require(stats))
suppressMessages(require(quantreg))
*/
// #.R_debug_Rcpp = FALSE

// [[Rcpp::plugins(openmp)]]
#include <Rcpp.h>
#include <omp.h>
using namespace Rcpp;

Function subset("[.data.frame");
Function tapply("tapply");
Function cbind("cbind.data.frame");
Function which("which");

Environment Rpkg_stats("package:stats");
Function Rf_kmeans = Rpkg_stats["kmeans"];
Function Rf_modmat = Rpkg_stats["model.matrix"];

Environment Rpkg_quantreg("package:quantreg");
// Function Rf_rq = Rpkg_quantreg["rq"];
Function Rf_rqLasso = Rpkg_quantreg["rq.fit.lasso"]; 
// SEXP getKnots(DataFrame data, int i, int maxit) {
// dst

inline double  epan(double u) {
  return fabs(u) <= 1.0 ? 0.75 * (1.0 - u * u): 0.0;
}

// NumericVector fitXXX(DataFrame data, IntegerVector i, IntegerVector maxit) {
//   List ob = Rf_kmeans(_["x"] = data, _["centers"] = i, _["iter.max"] = maxit);
//   NumericVector invsigma = as<NumericVector>(ob["size"]) / as<NumericVector>(ob["withinss"]);
//   return sqrt(invsigma);
// }

// [[Rcpp::export]]
List getBasis(DataFrame data, List km, double eps) {
  int v, k, bng, bns;
  int nk = as<NumericVector>(km["withinss"]).size();
  NumericVector sigma(nk);
  NumericVector tmp(data.nrows());
  NumericMatrix basis(data.nrows(), nk * data.size());
  NumericMatrix S(nk, data.size());
  NumericMatrix knots = as<NumericMatrix>(km["centers"]);
  List ans;

  bng = 0;
  bns = 0;
  for (v = 0; v < data.size(); v++) {
    // Scale estimation based on the groups
    sigma = tapply(data[v], km["cluster"], "sd");
    sigma = 1.0 / (sigma + eps);
    std::copy(sigma.begin(), sigma.end(), S.begin() + bns);
    bns += sigma.size();
    // Calculation of the basis
    for (k^= k; k < nk; k++) {
      tmp = (as<NumericVector>(data[v]) - knots(k, v)) * sigma[k];
      std::transform(tmp.begin(), tmp.end(), basis.begin() + bng, epan);
      bng += data.nrows();
    }
  }
  ans["sigmas"] = S;
  ans["basis"] = basis;
  return ans;
}

// [[Rcpp::export]]
NumericMatrix getBasisPred(DataFrame data, List km, NumericMatrix S, double eps) {
  int k, bng = 0, v = 0;
  int nk = as<NumericVector>(km["withinss"]).size();
  NumericVector tmp(data.nrows());
  NumericMatrix knots = as<NumericMatrix>(km["centers"]);
  NumericMatrix basis(data.nrows(), nk * data.size());

  for (; v < data.size(); v++) {
    for (k^= k; k < nk; k++) { // Calculation of the basis for prediction
      tmp = (as<NumericVector>(data[v]) - knots(k, v)) * S(k, v);
      std::transform(tmp.begin(), tmp.end(), basis.begin() + bng, epan);
      bng += data.nrows();
    }
  }
  return basis;
}

// [[Rcpp::export]]
List weightCalc(NumericMatrix x, NumericVector y, double scale, int type) {
  int i = 0, b = 0;
  NumericVector wts(y.size(), 1.0);
  NumericVector tmp(y.size());
  List ans;
  if (type) {
    for(; i < y.size(); i++) {
      switch(type) {
        case 2:
          wts[i] = (double) (y.size() - i - 1) / scale;
        case 1:
          wts[i] *= (double) (y.size() - i - 1) / scale;
          wts[i] = exp(-wts[i]);
      }
    }
    wts = wts / sum(wts);
    for (i = 0; i < x.ncol(); i++, b += x.nrow()) {
      std::copy(b + x.begin(), b + x.begin() + x.nrow(), tmp.begin());
      tmp = tmp * wts;
      std::copy(tmp.begin(), tmp.end(), b + x.begin());
    }
  }
  tmp = y * wts;
  ans["x"] = x;
  ans["y"] = tmp;
  return ans;
}

// [[Rcpp::export]]
NumericVector getForErr(Formula formula, DataFrame data, IntegerVector length, IntegerVector nknots, NumericVector scale, IntegerVector type) {
  int i, j;
  int nlen = length(0);
  int nrow = data.nrows();
  double obs;
  double eps = 1e-3;
  NumericVector err(nrow - nlen, 0.0);
  LogicalVector wcol(data.size());
  IntegerVector wo(nlen);
  NumericMatrix D;
  List basis;

  DataFrame subdata;

  // Selection of variables for non-linear dynamics
  for (j ^= j; j < nlen; j++) wo[j] = 1 + j;
  for (i = 1; i < data.size(); i++) {
    if (as<NumericVector>(unique(as<NumericVector>(data[i]))).size() > nknots[0]) {
      wcol[i] = TRUE;
    }
    else {
      wcol[i] = FALSE;
    }
  }
  NumericVector wc = which(wcol);
  // K-means for finding most informative knots as centroids
  subdata = subset(data, wo, wc);
  List ob = Rf_kmeans(_["x"] = subdata, _["centers"] = nknots[0], _["iter.max"] = 200);
  // Definition of the design matrix
  basis = getBasis(subdata, ob, eps);
  D = Rf_modmat(formula, _["data"] = cbind(subset(data, wo, seq_len(data.size())), basis["basis"]), _["subset"] = wo);
  // Weighted robust Lasso regression 
  NumericVector lambdas(D.ncol(), 1.0);
  List W = weightCalc(D, subset(data, wo, 1), scale[0], type[0]);
  List fit = Rf_rqLasso(_["x"] = W["x"], _["y"] = W["y"], _["tau"] = 0.5, _["lambda"] = lambdas);
  // Prediction
  subdata = subset(data, j + 1, wc);
  D = getBasisPred(subdata, ob, basis["sigmas"], eps);
  subdata = subset(data, j + 1, seq_len(data.size()));
  D = Rf_modmat(formula, _["data"] = cbind(subdata, D));
  obs = as<double>(subset(data, j + 1, 1));
  err[0] = as<NumericVector>(fit["coefficients"])(0);
  for (j = 1; j < D.ncol(); j++)
    err[0] += as<NumericVector>(fit["coefficients"])(j) * D(0, j);
  // Prediction errors
  err[0] = obs - err[0];
  for (i = 1; i < err.size(); i++) {
    // Selection of variables for non-linear dynamics
    for (j ^= j; j < nlen; j++) wo[j] = 1 + j + i;
    // K-means for finding most informative knots as centroids
    subdata = subset(data, wo, wc);
    ob = Rf_kmeans(_["x"] = subdata, _["centers"] = ob["centers"], _["iter.max"] = 200);
    // Definition of the design matrix
    basis = getBasis(subdata, ob, eps);
    D = Rf_modmat(formula, _["data"] = cbind(subset(data, wo, seq_len(data.size())), basis["basis"]), _["subset"] = wo);
    // Weighted robust Lasso regression 
    W = weightCalc(D, subset(data, wo, 1), scale[0], type[0]);    
    fit = Rf_rqLasso(_["x"] = W["x"], _["y"] = W["y"], _["tau"] = 0.5, _["lambda"] = lambdas);
    // Prediction
    subdata = subset(data, j + 1 + i, wc);
    D = getBasisPred(subdata, ob, basis["sigmas"], eps);
    subdata = subset(data, j + 1 + i, seq_len(data.size()));
    D = Rf_modmat(formula, _["data"] = cbind(subdata, D));
    obs = as<double>(subset(data, i + j + 1, 1));
    err[i] = as<NumericVector>(fit["coefficients"])(0);
    for (j = 1; j < D.ncol(); j++)
      err[i] += as<NumericVector>(fit["coefficients"])(j) * D(0, j);
    // Prediction errors
    err[i] = obs - err[i];
  }
  return err;
}

// library(Rcpp); sourceCpp("npfit.cpp")
