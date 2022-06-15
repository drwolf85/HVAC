## Loading packages and tools
suppressMessages(require(stats))
suppressMessages(require(quantreg))
suppressMessages(require(Rcpp))
suppressMessages(require(biglm))

## Compile the C++ code via RCpp
biglmPath <- system.file("libs", package = "biglm", mustWork = TRUE)
biglmPath <- paste0("-L", biglmPath, " -l:biglm.so `$(R_HOME)/bin/Rscript -e \"Rcpp:::LdFlags()\"`")
Sys.setenv("PKG_LIBS" = biglmPath)
sourceCpp("pareto.cpp")

## Loading the data
load("data.RData")
load("dnow.RData")
load("errs.RData")

## Possible actions of the system
actions <- expand.grid(fc_0 = c(0, .33, .66, 1), d1_0 = c(0, 100), d2_0 = c(0, 100), b_0 = c(0, 50, 100))

## Model settings
formulas <- c("x1_0 ~ .", "x2_0 ~ .:.", "yth ~ .:.")
HMMmeth <- c(FALSE, FALSE, TRUE) # specify if the response is a latent variable
tmwintyp <- c(0L, 0L, 0L) # specify the weighting-time-window
lengths <- c(288L * 7L * 3L, 288L * 7L * 3L, 288L * 7L * 3L)#29430L)
scales <- c(1, 1, 1)
knots <- c(10L, 10L, 10L)
datasets <- list(pmv, dgi, ec)
errs <- list(erf1o, erf2o, resEC)
errs <- sapply(errs, function(xx) tail(xx, min(lengths)))

## Estimation
estTime <- system.time(estim <- mapply(GlobalFit, formulas, HMMmeth, tmwintyp, lengths, scales, knots, datasets))
cat("Elapsed time for estimation:\n")
print(estTime)

## Prediction settings
nowdata <- list(opmv, odgi, t(as.matrix(oec[, -1])))
decisions <- list(data.frame(x1_0 = 0, fc_0 = actions[, 1]), cbind.data.frame(x2_0 = 0, actions[, -1]), cbind.data.frame(yth = 0 , actions[, -4], bl = abs(actions[, 4] - oec[, 1])))

## Prediction
predTime <- system.time(preds <- mapply(GlobalPredict, estim, HMMmeth, decisions, nowdata))
cat("Elapsed time for Prediction:\n")
print(predTime)

Lfun <- c(abs, function(x) abs(x - 11), function(x) ifelse(x < 0, 0, x))

losses <- preds
losses[, 1] <- abs(losses[, 1])
losses[, 2] <- abs(losses[, 2] - 11)

## Optimization based on Prediction (classical results with constraints)
allowFC <- FALSE
# allowFC <- (0.155 + 1.05 * cos(0.35 + 0.218 * data$hour)) < 0
allowDMX <- odgi[, "v8_1"] == 1 
lamps <- (actions$d1_0 > 0 | actions$d2_0 >0)
feasible <- xor(!(lamps | allowDMX), allowDMX) & xor(!(actions$fc_0 > 0 | allowFC), allowFC)
optimTime <- system.time(wpf <- ParetoFront(losses, feasible))
cat("Elapsed time for Optimization:\n")
print(predTime)

preds[wpf, ]
losses[wpf, ]
actions[wpf, ]

## Optimization based on Bootstrap without constraints

domprobs <- FuzzyFront(preds, errs, Lfun)

preds[which.max(domprobs), ]
losses[which.max(domprobs), ]
actions[which.max(domprobs), ]

## Optimization based on Bootstrap with constraints
domprobs <- FuzzyFront(preds, errs, Lfun, feasible)

preds[which.max(domprobs), ]
losses[which.max(domprobs), ]
actions[which.max(domprobs), ]
