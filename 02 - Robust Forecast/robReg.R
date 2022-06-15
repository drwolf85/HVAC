## load data and functions required
{ if (wb < 3.5) { load("pmv.RData") } else { load("dgi.RData") } }
suppressMessages(require(stats))
suppressMessages(require(quantreg))
suppressMessages(require(Rcpp))
sourceCpp("npfit.cpp")

#1 erf1o <- getForErr("x1_0 ~ .", dat, 288L * 7L * 3L, 10L, 1, 0L)
#2 erf1e <- getForErr("x1_0 ~ .", dat, 288L * 7L * 3L, 10L, 288, 1L)
#3 erf1n <- getForErr("x1_0 ~ .", dat, 288L * 7L * 3L, 10L, 288, 2L)
#4 erf2o <- getForErr("x2_0 ~ .:.", dat, 288L * 7L * 3L, 10L, 1, 0L)
#5 erf2e <- getForErr("x2_0 ~ .:.", dat, 288L * 7L * 3L, 10L, 288, 1L)
#6 erf2n <- getForErr("x2_0 ~ .:.", dat, 288L * 7L * 3L, 10L, 288, 2L)

#1 save(erf1o, file = "erf1o.RData")
#2 save(erf1e, file = "erf1e.RData")
#3 save(erf1n, file = "erf1n.RData")
#4 save(erf2o, file = "erf2o.RData")
#5 save(erf2e, file = "erf2e.RData")
#6 save(erf2n, file = "erf2n.RData")


