library(Rcpp)
sourceCpp("benchforex.cpp")

robust <- new.env()
load("errbustforex.RData", robust)
arimax <- new.env()
load("errARIMAX.RData", arimax)

files <- grep("[dgpv]", dir(, ".RData"), value = TRUE)
gofoR <- sapply(ls(robust), function(x) {
  myenv <- new.env()
  assign("e", get(x, envir = robust), envir = myenv)
  load(files[grepl("2", x) + 1])
  assign("o", dat[, 1L], envir = myenv)
  ans <- c(eval(quote(MAE(e)), envir = myenv),
           eval(quote(MAPE(e, o)), envir = myenv),
           eval(quote(sMAPE(e, o)), envir = myenv),
           eval(quote(msMAPE(e, o)), envir = myenv),
           eval(quote(RMSE(e)), envir = myenv),
           eval(quote(NMSE(e, o)), envir = myenv),
           eval(quote(KLN(e, o)), envir = myenv),
           eval(quote(KLDE1(e, o)), envir = myenv),
           eval(quote(KLDE2(e, o)), envir = myenv),
           eval(quote(KLIQR(e, o)), envir = myenv),
           eval(quote(RSE(e, o)), envir = myenv),
           eval(quote(m1RSE(e, o)), envir = myenv),
           eval(quote(m2RSE(e, o)), envir = myenv),
           eval(quote(TU2(e, o)), envir = myenv),
           eval(quote(RAE(e, o)), envir = myenv),
           eval(quote(MSEr1(e)), envir = myenv),
           eval(quote(MSEr2(e, o)), envir = myenv))
  return(ans)})
rownames(gofoR) <- c("MAE", "MAPE", "sMAPE", "msMAPE", "RMSE", "NMSE", "KLN", "KLDE1", "KLDE2", "KLIQR", "RSE", "m1RSE", "m2RSE", "TU2", "RAE", "MSEr1", "MSEr2")

gofoA <- sapply(ls(arimax), function(x) {
  myenv <- new.env()
  assign("e", get(x, envir = arimax), envir = myenv)
  load(files[grepl("2", x) + 1])
  assign("o", dat[, 1L], envir = myenv)
  ans <- c(eval(quote(MAE(e)), envir = myenv),
           eval(quote(MAPE(e, o)), envir = myenv),
           eval(quote(sMAPE(e, o)), envir = myenv),
           eval(quote(msMAPE(e, o)), envir = myenv),
           eval(quote(RMSE(e)), envir = myenv),
           eval(quote(NMSE(e, o)), envir = myenv),
           eval(quote(KLN(e, o)), envir = myenv),
           eval(quote(KLDE1(e, o)), envir = myenv),
           eval(quote(KLDE2(e, o)), envir = myenv),
           eval(quote(KLIQR(e, o)), envir = myenv),
           eval(quote(RSE(e, o)), envir = myenv),
           eval(quote(m1RSE(e, o)), envir = myenv),
           eval(quote(m2RSE(e, o)), envir = myenv),
           eval(quote(TU2(e, o)), envir = myenv),
           eval(quote(RAE(e, o)), envir = myenv),
           eval(quote(MSEr1(e)), envir = myenv),
           eval(quote(MSEr2(e, o)), envir = myenv))
  return(ans)})
rownames(gofoA) <- rownames(gofoR)

save(gofoA, gofoR, file = "benchmarks.RData")


alpha = .05
vR <- sapply(ls(robust), function(x) quantile(abs(get(x, envir = robust)), 1 - alpha))
bR <- sapply(ls(robust), function(x) abs(median(get(x, envir = robust))))
vA <- sapply(ls(arimax), function(x) quantile(abs(get(x, envir = arimax)), 1 - alpha))
bA <- sapply(ls(arimax), function(x) abs(median(get(x, envir = arimax))))

m1 <- cbind(c(vR[1:3], vA[1]), c(bR[1:3], bA[1]))
m2 <- cbind(c(vR[-1:-3], vA[-1]), c(bR[-1:-3], bA[-1]))
m1 <- m1[order(m1[, 1]), ]
m2 <- m2[order(m2[, 1]), ]

w1 <- ParetoFrontUnconstr(m1)
w2 <- ParetoFrontUnconstr(m2)

## Pareto Front for model selection
pdf("wtwFront.pdf")
plot(rbind(m1, m2), type = "n", xlab = "Variability", ylab = "Bias", main = "Pareto Front: weighted-time-window selection")
points(as.data.frame(m1)[!w1, ], col = 8)
points(as.data.frame(m1)[w1, ], type = "b")
points(as.data.frame(m2)[!w2, ], pch = 3, col = 8)
points(as.data.frame(m2)[w2, ], pch = 3, type = "b")
legend("center", c("Dominant point (PMV model)", "Dominated point (PMV model)", "Dominant point (DGI-model)", "Dominated point (DGI-model)"), pch = c(1, 1, 3, 3), col = c(1, 8, 1, 8))
dev.off()

## Stochastic dominance??
E1 <- cbind(sapply(ls(robust)[1:3], get, envir = robust), arimax$arimax1)
E2 <- cbind(sapply(ls(robust)[-1:-3], get, envir = robust), arimax$arimax2)

prop.table(table(c(ls(robust)[1:3], "arimax1")[apply(abs(E1), 1, which.min)]))
prop.table(table(c(ls(robust)[-1:-3], "arimax2")[apply(abs(E2), 1, which.min)]))

