#include <R.h>
#include <Rinternals.h>

SEXP cumSumOnZeros(SEXP x) {
  int i, n;
  double val, *px;

  PROTECT(x = coerceVector(x, REALSXP));
  px = REAL(x);
  n = length(x);
  val = px[0];
  for (i = 1; i < n; i++) {
    if (px[i] > 0.0) {
      val = px[i];
    }
    else {
      px[i] = val;
    }
  }
  UNPROTECT(1);
  return x;
}
