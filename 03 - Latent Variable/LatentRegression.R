rm(list = ls(all = TRUE))
graphics.off()

dataset <- "dati_Uff2_stratega2.csv"
# dataset <- "dati_Uff2_stratega.csv"
source("readData.R")
.missing <- c("v3", "v9", "v10", "v12")
source("missing.R")

.determ_fit <- FALSE

## Data ajdustments
data$d1[data$d1 < 0] <- 0
data$d2[data$d2 < 0] <- 0
data$fc[data$fc < 0] <- 0
data$d1[data$d1 > 0] <- 100
data$d2[data$d2 > 0] <- 100
data$fc[data$fc > 1] <- 1
bl <- c(0, abs(diff(data$b))) # blinder
# data$y_load[data$y_load < 0] <- 0

## E-step
wh <- data$y > 0
rt <- diff(cumsum(data$y)[unique(c(1L, which(wh)))]) / diff(seq_along(wh)[unique(c(1L, which(wh)))])
tmp <- data$y + 0
tmp[unique(c(1L, which(wh)))] <- c(rt, 0)
dyn.load("EMalgo.so")
invisible(.Call("cumSumOnZeros", tmp))
thresh <- optimize(function(x) (sum(tmp[tmp > x]) - sum(data$y))**2, c(0, .1))
tmp[tmp <= thresh$minimum] <- 0

## M-step
formula <- formula(tmp ~ bl + bl:. + .:.)
suppressMessages(require(stats))
suppressMessages(require(quantreg))
suppressMessages(require(Rcpp))
suppressMessages(require(biglm))
sourceCpp("npfit.cpp")
# Selection of useful data
datsel <- cbind(data[-1L, c("fc", "d1", "d2")], bl = bl[-1L], data[-NROW(data), c("x1", "x2", "y_load", "v5", "v6", "v8")])
# Non-parametric set of basis function got with k-means
set.seed(0)
km <- kmeans(datsel[, c("x1", "x2", "y_load", "v5", "v6")], 10, iter.max = 200)
mybasis <- getBasis(datsel[, c("x1", "x2", "y_load", "v5", "v6")], km, 0.01)
X <- model.matrix(yth ~ .:., cbind(yth = tmp[-1L], datsel, mybasis$basis))
# Estimate the model coefficients by successive update of the QR-matrix
mybounds <- seq(1, NROW(X), by = NCOL(X) + 1)
mybounds[length(mybounds)] <- NROW(X)
wh <- apply(embed(mybounds, 2)[, 2:1], 1, function(x) x[1]:x[2])
qr <- biglm:::bigqr.init(NCOL(X))
for (i in seq_along(wh)) {
  qr <- .Call("updateQR", X[wh[[i]], ], tmp[1 + wh[[i]]], rep(1, length(wh[[i]])), qr, FALSE, PACKAGE = "biglm")
}
qr <- .Call("singcheckQR", qr, PACKAGE = "biglm")
coef <- .Fortran("regcf", as.integer(length(qr$D)), as.integer(length(qr$D)^2 * 0.5), qr$D, qr$rbar, qr$thetab, qr$tol, beta = numeric(length(qr$D)), nreq = as.integer(length(qr$D)), ier = integer(1), PACKAGE = "biglm")$beta
# Calculation of fitted values
fitted <- X %*% coef
fitted[fitted < 0] <- 0

best <- function(xxx) mean(abs(data$y[-1] - diff((floor((xxx + cumsum(fitted)) * 10))) / 10))
rEnergyCons <- optimize(best, 0:1 * 0.1)

# rqfit <- rq(tmp ~ .:., data = cbind(data[, c("fc", "d1", "d2", "x1", "x2", "y_load", "v5", "v6", "v8")], bl = bl, mybasis$basis), method="lasso")
