// [[Rcpp::plugins(openmp)]]
#include <Rcpp.h>
#include <omp.h>
using namespace Rcpp;

Environment Rpkg_stats("package:stats");
Function Rf_IQR = Rpkg_stats["IQR"];

/*********************************
 * STAND-ALONE ACCURACY MEASURES *
 *********************************/

// [[Rcpp::export]]
double MAE(NumericVector e) {
  return mean(abs(e));
}

// [[Rcpp::export]]
double MAPE(NumericVector e, NumericVector o) {
  int len = o.size() - e.size();
  if (len < 0) ::Rf_error("observations cannot be less than prediction errors");
  int i = 0, rm = 0;
  double ans = 0.0;
  for (; i < e.size(); i++) {
    if (o[i + len] != 0.0) {
      ans += fabs(e[i] / o[i + len]);
    } 
    else {
      ++rm;
    }
  }
  ans /= (double) (e.size() - rm);
  return ans;
}

// [[Rcpp::export]]
double sMAPE(NumericVector e, NumericVector o) {
  int len = o.size() - e.size();
  if (len < 0) ::Rf_error("observations cannot be less than prediction errors");
  int i = 0, rm = 0;
  double ans = 0.0, deno;
  for (; i < e.size(); i++) {
    deno = fabs(o[i + len]);
    deno += fabs(o[i + len] - e[i]);
    if (deno != 0.0) {
      ans += fabs(e[i]) / (deno);
    } 
    else {
      ++rm;
    }
  }
  ans /= 0.5 * (double) (e.size() - rm);
  return ans;
}

// [[Rcpp::export]]
double msMAPE(NumericVector e, NumericVector o) {
  int len = o.size() - e.size();
  if (len < 0) ::Rf_error("observations cannot be less than prediction errors");
  int i = 0, rm = 0;
  double ans = 0.0, deno, si = 0.0, ybar = 0.0;
  for (; i < e.size(); i++) {
    ybar *= (double) i;
    ybar += o[i + len];
    ybar /= (double) (i + 1);
    si *= (double) i;
    si += fabs(o[i + len] - ybar);
    si /= (double) (i + 1);
    deno = fabs(o[i + len]);
    deno += fabs(o[i + len] - e[i]);
    deno *= 2;
    deno += si;
    if (deno != 0.0) {
      ans += fabs(e[i] / deno);
    } 
    else {
      ++rm;
    }
  }
  ans /= (double) (e.size() - rm);
  return ans;
}

// [[Rcpp::export]]
double RMSE(NumericVector e) {
  return sqrt(mean(e * e));
}

// [[Rcpp::export]]
double NMSE(NumericVector e, NumericVector o) {
  int len = o.size() - e.size();
  if (len < 0) ::Rf_error("observations cannot be less than prediction errors");
  NumericVector eo(e.size());
  std::copy(o.begin() + len, o.end(), eo.begin());
  eo = eo - mean(eo);
  return sqrt(mean(e * e)/mean(eo * eo));
}

// [[Rcpp::export]]
double KLN(NumericVector e, NumericVector o) {
  int len = o.size() - e.size();
  if (len < 0) ::Rf_error("observations cannot be less than prediction errors");
  int i = 0, rm = 0;
  double ans = 0.0, deno, si = 0.0, ybar = 0.0;
  for (; i < e.size(); i++) {
    ybar *= (double) i;
    ybar += o[i + len];
    ybar /= (double) (i + 1);
    si *= (double) i;
    deno = o[i + len] - ybar;
    si += deno * deno;
    si /= (double) (i + 1);
    if (si != 0.0) {
      ans += e[i] * e[i] / si;
    }
    else {
      ++rm;
    }
  }
  ans /= (double) (e.size() - rm);
  return sqrt(ans);
}

// [[Rcpp::export]]
double KLDE1(NumericVector e, NumericVector o) {
  int len = o.size() - e.size();
  if (len < 0) ::Rf_error("observations cannot be less than prediction errors");
  NumericVector ae = abs(e);
  int i = 0, rm = 0;
  double ans = 0.0, deno, si = 0.0, ybar = 0.0;
  for (; i < e.size(); i++) {
    ybar *= (double) i;
    ybar += o[i + len];
    ybar /= (double) (i + 1);
    si *= (double) i;
    deno = o[i + len] - ybar;
    si += fabs(deno);
    si /= (double) (i + 1);
    if (si != 0.0) {
      deno = -ae[i] / si;
      ans += exp(deno) - deno - 1.0;
    }
    else {
      ++rm;
    }
  }
  ans /= (double) (e.size() - rm);
  return ans;  
}

// [[Rcpp::export]]
double KLDE2(NumericVector e, NumericVector o) {
  int len = o.size() - e.size();
  if (len < 0) ::Rf_error("observations cannot be less than prediction errors");
  NumericVector ae = abs(e);
  int i = 0, rm = 0;
  double ans = 0.0, deno, si = 0.0, ybar = 0.0;
  for (; i < e.size(); i++) {
    ybar *= (double) i;
    ybar += o[i + len];
    ybar /= (double) (i + 1);
    si *= (double) i;
    deno = o[i + len] - ybar;
    si += deno * deno;
    si /= (double) (i + 1);
    if (si != 0.0) {
      deno = -ae[i] / sqrt(si);
      ans += exp(deno) - deno - 1.0;
    }
    else {
      ++rm;
    }
  }
  ans /= (double) (e.size() - rm);
  return ans;  
}

// [[Rcpp::export]]
double KLIQR(NumericVector e, NumericVector o) {
  int len = o.size() - e.size();
  if (len < 0) ::Rf_error("observations cannot be less than prediction errors");
  NumericVector eo(e.size());
  std::copy(o.begin() + len, o.end(), eo.begin());
  return RMSE(e) / as<double>(Rf_IQR(eo));
}

/******************************
 * RELATIVA ACCURACY MEASURES *
 ******************************/

// [[Rcpp::export]]
double RSE(NumericVector e, NumericVector o) {
  int len = o.size() - e.size();
  if (len < 1) ::Rf_error("length(observations) >=  length('prediction errors') + 1");
  int i = 0, rm = 0;
  double ans = 0.0, si = 0.0;
  for (; i < e.size(); i++) {
    si = o[i + len] - o[i + len -1];
    si *= si;
    if (si != 0.0) {
      ans += e[i] * e[i] / si;
    }
    else {
      ++rm;
    }
  }
  ans /= (double) (e.size() - rm);
  return sqrt(ans);
}

// [[Rcpp::export]]
double m1RSE(NumericVector e, NumericVector o) {
  int len = o.size() - e.size();
  if (len < 1) ::Rf_error("length(observations) >= length('prediction errors') + 1");
  int i = 0, rm = 0;
  double ans = 0.0, ar, ybar = 0.0, si = 0.0;
  for (; i < e.size(); i++) {
    ybar *= (double) i;
    ybar += o[i + len];
    ybar /= (double) (i + 1);
    si *= (double) i;
    ar = o[i + len] - ybar;
    si += ar * ar;
    si /= (double) (i + 1);
    ar = o[i + len] - o[i + len -1];
    ar *= ar;
    ar += si;
    if (ar != 0.0) {
      ans += e[i] * e[i] / ar;
    }
    else {
      ++rm;
    }
  }
  ans /= (double) (e.size() - rm);
  return sqrt(ans);
}

// [[Rcpp::export]]
double m2RSE(NumericVector e, NumericVector o) {
  int len = o.size() - e.size();
  if (len < 1) ::Rf_error("length(observations) >= length('prediction errors') + 1");
  int i = 0, rm = 0;
  double ans = 0.0, ar, si = 0.0;
  for (; i < e.size(); i++) {
    si *= (double) i;
    ar = o[i + len] - o[i + len -1];
    si += ar * ar;
    si /= (double) (i + 1);
    ar *= ar;
    ar += si;
    if (ar != 0.0) {
      ans += e[i] * e[i] / ar;
    }
    else {
      ++rm;
    }
  }
  ans /= (double) (e.size() - rm);
  return sqrt(ans);
}

// [[Rcpp::export]]
double TU2(NumericVector e, NumericVector o) {
  NumericVector od = diff(o);
  int len = od.size() - e.size();
  if (len < 1) ::Rf_error("length(observations) >= length('prediction errors') + 1");
  NumericVector eo(e.size());
  std::copy(od.begin() + len, od.end(), eo.begin());
  return sqrt(sum(e * e) / sum(eo * eo));
}

// [[Rcpp::export]]
double RAE(NumericVector e, NumericVector o) {
  NumericVector od = diff(o);
  int len = od.size() - e.size();
  if (len < 1) ::Rf_error("length(observations) >= length('prediction errors') + 1");
  NumericVector eo(e.size());
  std::copy(od.begin() + len, od.end(), eo.begin());
  return sqrt(sum(abs(e)) / sum(abs(eo)));
}

// [[Rcpp::export]]
double MSEr1(NumericVector e) {
  double ans = 0.0, ms = 0.0, e2;
  int i = 0, rm = 0;
  for (; i < e.size(); i++) {
    e2 = e[i] * e[i];
    ms *= (double) i;
    ms += e2;
    ms /= (double) (i + 1);
    if (ms != 0.0) {
      ans += e2 / ms; 
    }
    else {
      ++rm;
    }
  }
  ans /= (double) (e.size() - rm);  
  return sqrt(ans);
}

// [[Rcpp::export]]
double MSEr2(NumericVector e, NumericVector o) {
  int len = o.size() - e.size();
  if (len < 0) ::Rf_error("observations cannot be less than prediction errors");
  double ans = 0.0, ms = 0.0, ybar = 0.0, e2;
  int i = 0, rm = 0;
  for (; i < e.size(); i++) {
    e2 = o[i + len] - e[i];
    ybar *= (double) i;
    ybar += e2;
    ybar /= (double) (i + 1);
    e2 = o[i + len] - ybar;
    e2 *= e2;
    ms *= (double) i;
    ms += e2;
    ms /= (double) (i + 1);
    if (ms != 0.0) {
      ans += e[i] * e[i] / ms; 
    }
    else {
      ++rm;
    }
  }
  ans /= (double) (e.size() - rm);  
  return sqrt(ans);
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
LogicalVector ParetoFrontUnconstr(NumericMatrix Fs) {
  int j, k;
  bool dom, sdom;
  LogicalVector ans(Fs.nrow(), true);

  for (j ^= j; j < Fs.nrow(); j++) {
    dom = true;
    for (k ^= k; k < j; k++) {
      sdom = isDominant(Fs(j, _), Fs(k, _));
      dom &= sdom;
      ans[k] &= !sdom;
    }
    for (k = j + 1; k < Fs.nrow(); k++) {
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
