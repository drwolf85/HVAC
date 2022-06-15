#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>

extern "C"
{
  SEXP updateQR(SEXP, SEXP, SEXP, SEXP, SEXP);
  SEXP singcheckQR(SEXP);
  void F77_NAME(regcf)(int *, int *, double *, double *, double *, double *, double *, int *, int *);
}

// [[Rcpp::depends(biglm)]]

#include <Rcpp.h>

using namespace Rcpp;

Function subset("[.data.frame");
Function tapply("tapply");
Function cbind("cbind.data.frame");
Function which("which");
Function Rf_dotCall(".Call");

Environment Rpkg_stats("package:stats");
Function Rf_kmeans = Rpkg_stats["kmeans"];
Function Rf_modmat = Rpkg_stats["model.matrix"];

Environment Rpkg_quantreg("package:quantreg");
Function Rf_rqLasso = Rpkg_quantreg["rq.fit.lasso"];

inline double  epan(double u) {
  return fabs(u) <= 1.0 ? 0.75 * (1.0 - u * u): 0.0;
}

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
List GlobalFit(Formula formula, LogicalVector regmeths, IntegerVector type, IntegerVector mylength, NumericVector scale, IntegerVector nknots, DataFrame data) {
  int j;
  int nlen = mylength(0);
  int nrow = data.nrows();
  int nobs = nrow - nlen;
  double eps = 1e-3;
  LogicalVector wcol(data.size(), false);
  IntegerVector wo(nlen);
  NumericMatrix D;
  DataFrame subdata;
  List basis;
  List fit;
  // Selection of variables for non-linear dynamics
  for (j ^= j; j < nlen; j++) wo[j] = j + nobs;
  for (j = 1; j < data.size(); j++) {
    if (as<NumericVector>(unique(as<NumericVector>(data[j]))).size() > nknots[0]) {
      wcol[j] = true;
    }
  }
  NumericVector wc = which(wcol);
  // K-means for finding most informative knots as centroids
  subdata = subset(data, wo, wc);
  List ob = Rf_kmeans(_["x"] = subdata, _["centers"] = nknots[0], _["iter.max"] = 200);  
  // Definition of the design matrix
  basis = getBasis(subdata, ob, eps);
  D = Rf_modmat(formula, _["data"] = cbind(subset(data, wo, seq_len(data.size())), basis["basis"]), _["subset"] = wo);
  if (!regmeths(0)) {
    // Weighted robust Lasso regression 
    NumericVector lambdas(D.ncol(), 1.0);
    List W = weightCalc(D, subset(data, wo, 1), scale[0], type[0]);
    fit = Rf_rqLasso(_["x"] = W["x"], _["y"] = W["y"], _["tau"] = 0.5, _["lambda"] = lambdas);
  }
  else {
    int ier = D.nrow() / (D.ncol() + 1);
    int nb22 = D.ncol() * D.ncol() / 2;
    int npar;
    List qr;
    List newqr(ier + 2);
    SEXP mat, resp, weights;
    // qr <- biglm::bigqr.init(NCOL(D))
    SEXP addint = wrap(false);
    qr["D"] = NumericVector(D.ncol(), 0.0);
    qr["rbar"] = NumericVector(Rf_choose(D.ncol(), 2.0));
    qr["thetab"] = NumericVector(D.ncol(), 0.0);
    qr["ss"] = 0.0;
    qr["checked"] = addint;
    qr["tol"] = NumericVector(D.ncol(), 0.0);
    newqr(0) = wrap(qr);
    for (j ^= j; j <= ier; j++) {
      nlen = j * (D.ncol() + 1);
      npar = nlen + D.ncol();
      if (npar >= D.nrow()) npar = D.nrow() - 1;
      PROTECT(weights = wrap(rep(1.0, npar - nlen + 1)));
      PROTECT(resp = subset(data, seq(nlen + 1, npar + 1), 1));
      PROTECT(mat = wrap(D(Range(nlen, npar), _)));
      newqr(j + 1) = updateQR(mat, resp, weights, newqr(j), addint);
      UNPROTECT(3);
    }
    qr = as<List>(singcheckQR(newqr(ier + 1)));
    ier = 1;
    NumericVector coef(D.ncol(), 0.0);
    npar = D.ncol();
    F77_CALL(regcf)(&npar, &nb22, REAL(wrap(qr["D"])), REAL(wrap(qr["rbar"])), REAL(wrap(qr["thetab"])), REAL(wrap(qr["tol"])), coef.begin(), &npar, &ier);
    fit["coefficients"] = coef;
  }
  fit["wc"] = wc;
  fit["sigmas"] = basis["sigmas"];
  fit["centers"] = ob["centers"];
  fit["withinss"] = ob["withinss"];
  fit["formula"] = formula;
  return fit;
}

// [[Rcpp::export]]
NumericVector GlobalPredict(List fit, LogicalVector regmeths, DataFrame actions, NumericMatrix odata) {
  int i, j;
  double eps = 1e-3;
  NumericVector fitted(actions.nrows(), 0.0);
  // Organize the data and action for prediction
  SEXP newdata = cbind(actions, odata);
  DataFrame subdata = subset(newdata, 1, fit["wc"]);
  // Construction of the design matrix
  NumericMatrix D = getBasisPred(subdata, fit, fit["sigmas"], eps);
  D = Rf_modmat(fit["formula"], _["data"] = cbind(newdata, D));
  for (j ^= j; j < D.nrow(); j++) {
    for (i ^= i; i < D.ncol(); i++) {
      fitted[j] += D(j, i) * as<NumericVector>(fit["coefficients"])(i);
    }
    if (regmeths[0]) {
      if (fitted[j] < 0.0) fitted[j] = 0.0;
    }
  }
  return fitted;
}

inline bool isDominant(NumericVector F1, NumericVector F2) {
  if (F1.size() != F2.size()) ::Rf_error("multi-objective vectors must have the same length");
  bool dom = true;
  bool som = false;
  int i;
  for (i ^= i; i < F1.size(); i++) {
    if (F1(i) > F2(i)) dom = false;
    if (F1(i) < F2(i)) som = true;
  }
  return dom & som;
}

// [[Rcpp::export]]
LogicalVector ParetoFront(NumericMatrix Fs, LogicalVector cs) {
  int j, k;
  bool dom, sdom;
  LogicalVector ans(Fs.nrow(), true);

  for (j ^= j; j < Fs.nrow(); j++) {
    if (!cs[j]) {
      ans[j] = false;
      continue;
    }
    dom = true;
    for (k ^= k; k < j; k++) {
      if (!cs[k]) continue;
      sdom = isDominant(Fs(j, _), Fs(k, _));
      dom &= sdom;
      ans[k] &= !sdom;
    }
    for (k = j + 1; k < Fs.nrow(); k++) {
      if (!cs[k]) continue;
      sdom = isDominant(Fs(j, _), Fs(k, _));
      dom &= sdom;
      ans[k] &= !sdom;
    }
    if (dom) {
      for (++j; j < Fs.nrow(); j++) {
        ans[j] = false;
      }
      return ans;
    }
  }
  return ans;
}

/*** R

FuzzyFront <- function(preds, errors, Lfun, feasible = NULL, B = 100) {
  if (is.null(feasible)) feasible <- rep(TRUE, nrow(preds))
  res <- replicate(B, {
    epsB <- sample(nrow(errors), nrow(preds), replace = TRUE)
    tmp <- preds + errors[epsB, ]
    tmp <- sapply(seq_along(Lfun), function(i) Lfun[[i]](tmp[, i]))
    return(ParetoFront(tmp, feasible))
  })
  return(rowMeans(res))
}

*/