/*** R
suppressMessages(require(stats))
*/
// #.R_debug_Rcpp = FALSE

// [[Rcpp::plugins(openmp)]]
#include <Rcpp.h>
#include <omp.h>
using namespace Rcpp;

Function subset("[.data.frame");
Function cbind("cbind");

Environment Rpkg_stats("package:stats");
Function Rf_lm = Rpkg_stats["lm"];
Function Rf_arima = Rpkg_stats["arima"];
Function Rf_predict = Rpkg_stats["predict"];
Function Rf_predictLM = Rpkg_stats["predict.lm"];

// [[Rcpp::export]]
NumericVector getARIMAXForErr(Formula formula, IntegerVector order, SEXP seasonal, DataFrame xreg, IntegerVector length) {
  int i, j;
  int nlen = length(0);
  int nrow = xreg.nrows();
  NumericVector err(nrow - nlen, 0.0);
  NumericVector init(order(0) + order(2), 0.0);
  NumericVector pred(1);
  NumericVector obs(1);
  IntegerVector wo(nlen);
  List fitLM;
  SEXP fitTS;

  for (i ^= i; i < err.size(); i++) {
    // Selection of variables for the linear dynamics and definition of the objects
Rprintf("%d ", i);
    for (j ^= j; j < nlen; j++) wo[j] = 1 + j + i;
    // ARIMAX regression
    fitLM = Rf_lm(formula, _["data"] = subset(xreg, wo, seq_len(xreg.size())));
Rprintf("LM ");
    if (Rf_isNull(seasonal)) {
      fitTS = Rf_arima(_["x"] = as<List>(fitLM)["residuals"], _["order"] = order, _["method"] = "CSS");
    } else {
      fitTS = Rf_arima(_["x"] = as<List>(fitLM)["residuals"], _["order"] = order, _["seasonal"] = seasonal, _["method"] = "CSS");
    }
Rprintf("TS ");
    // Prediction
    obs = subset(xreg, i + j + 1, 1);
    err[i] = obs[0];
    pred = Rf_predict(_["object"] = fitTS, _["n.ahead"] = 1, _["se.fit"] = 0);
Rprintf("ts.pred ");
    err[i] -= pred[0];
    pred = Rf_predictLM(_["object"] = fitLM, _["newdata"] = subset(xreg, i + j + 1, -1));
Rprintf("lm.pred ");
    err[i] -= pred[0];
Rprintf("EOL\n");
// if (i > 3) break;
  }
  return err;
}

// library(Rcpp); sourceCpp("arimax.cpp")
