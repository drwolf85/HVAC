## load data and functions required
{ if (wb < 1.5) { load("pmv.RData") } else { load("dgi.RData") } }
suppressMessages(require(stats))
suppressMessages(require(Rcpp))
sourceCpp("arimax.cpp")

ndata <- 288L * 7L * 3L

default.warn <- getOption("warn")
options(warn = -1)

arimax1 <- getARIMAXForErr(x1_0 ~ ., c(0, 1, 2), NULL, dat, ndata)
save(arimax1, file = "arimax1.RData")
arimax2 <- getARIMAXForErr(x2_0 ~ ., c(0, 1, 1), NULL, dat, ndata)
save(arimax2, file = "arimax2.RData")
options(warn = default.warn)

