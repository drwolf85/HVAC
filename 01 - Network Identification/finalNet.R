rm(list = ls(all = TRUE))
graphics.off()

library(Rcpp)

## Network analysis
source("Entropy.R")

dataset <- "dati_Uff2_stratega2.csv"
source("readData.R")
.missing <- c("v3", "v9", "v10", "v12")
source("missing.R")
sourceCpp("fit.cpp")


## Data ajdustments
data$d1[data$d1 < 0] <- 0
data$d2[data$d2 < 0] <- 0
data$fc[data$fc < 0] <- 0
data$d1[data$d1 > 100] <- 100
data$d2[data$d2 > 100] <- 100
data$fc[data$fc > 1] <- 1

## Paramters for the analysis
steps <- 1

## Selection of light-related variables
nomivar <- c("x2", "v5", "v6", "v8", "v11", "d1", "d2", "b")
dat <- data[, nomivar]
datname <- names(dat)
dat <- as.matrix(dat)
colnames(dat) <- datname
datcolo <- c("#00ccff", rep("#ffcc00", 4), rep("#00ff00", 3))
dat <- apply(dat, 2, quantize)
dat <- embed(dat, steps + 1)[, -2L:-5L]

res <- combn(ncol(dat), 2, function(xx) mi_qv(dat[, xx[1]], dat[, xx[2]]), simplify = TRUE)
attr(res, "Size") <- ncol(dat)
attr(res, "call") <- ""
attr(res, "Diag") <- FALSE
attr(res, "Upper") <- FALSE
class(res) <- "dist"

## Organize the variables in a network
library(igraph)
g <- graph.full(ncol(dat), directed = FALSE)
{if (steps == 0) {
  V(g)$label <- datname
  V(g)$color <- datcolo
}
else {
  V(g)$label <- paste(rep(datname, times = steps + 1), rep(0:steps, each = length(datname)), sep = "_")[-2L:-5L]
  V(g)$color <- rep(datcolo, times = steps + 1)[-2L:-5L]
}}
E(g)$color <- rgb(0, 0, 0, 1 - res)
#thresh <- getThresh(dim(dat), probs = 0.05, nsim = 100L) # critical value with significance level at 0.05 
thresh <- 0.979201112359260772422 
g2 <- delete.edges(g, E(g)[res > thresh])
g2 <- delete.vertices(g2, which(degree(g2) < 1))

## Modello
dat <- data[, nomivar]
dat <- as.matrix(dat)
dat <- embed(dat, steps + 1)[, -2L:-5L]
colnames(dat) <- V(g)$label
dat <- as.data.frame(dat)

## Fitting via lm #2 MODIFY WITHOUT THE TIME T-0
alf <- anova(linfit <- lm(x2_0 ~ .:.:., data = dat))
pe <- rev(exp(-seq_len(nrow(dat)) / 288))
pe <- pe / sum(pe)
alfw <- anova(linfitW <- lm(x2_0 ~ .:.:., data = dat, weights = pe))

#erfo <- getForErrUnif("x2_0 ~ .:.:.", dat, 288L * 3L)
#erfe <- getForErrExp("x2_0 ~ .:.:.", dat, 288L * 3L, scale = 288)
#erfn <- getForErrNorm("x2_0 ~ .:.:.", dat, 288L * 3L, scale = 288)

# ## Fitting via GA
# system.time(fits <- replicate(10L, nnet(x2_0 ~ ., data = dat, linout = TRUE, size = 4, maxit = 1000, trace = FALSE), simplify = FALSE))
# pos <- order(sapply(fits, "[[", "value"))
# fits[[pos[10L]]]$wts <- 0.5 * (fits[[pos[1L]]]$wts + fits[[pos[2L]]]$wts)
# sd <- sd(fits[[pos[1L]]]$wts)
# mu <- fits[[pos[1L]]]$wts
# for (i in 3L:9L) fits[[pos[i]]] <- nnet(x2_0 ~ ., data = dat, wts = rnorm(length(mu), mu, sd), linout = TRUE, size = 4, maxit = 1000, trace = FALSE)
# fits[[pos[10L]]] <- nnet(x2_0 ~ ., data = dat, wts = fits[[pos[10L]]]$wts, linout = TRUE, size = 4, maxit = 1000, trace = FALSE)
# pos <- order(sapply(fits, "[[", "value"))
# ## GA # fast implementation
# fit <- fits[[pos[1L]]]
# bestR2_x2 <- 1 - fit$value / (var(dat[,"x2_0"]) * (nrow(dat) - 1))


## Paramters for the analysis
steps <- 4

## Selection of temperature-related variables
nomivar <- c("x1", paste0("v", c(1:2, 4, 7, 13:15)), "y_load", "fc")
dat <- data[, nomivar]
datname <- names(dat)
dat <- as.matrix(dat)
colnames(dat) <- datname
datcolo <- c("#00ccff", rep("#ffcc00", 8), rep("#00ff00", 1))
dat <- apply(dat, 2, quantize)
dat <- embed(dat, steps + 1)[, -2L:-9L]

res <- combn(ncol(dat), 2, function(xx) mi_qv(dat[, xx[1]], dat[, xx[2]]), simplify = TRUE)
attr(res, "Size") <- ncol(dat)
attr(res, "call") <- ""
attr(res, "Diag") <- FALSE
attr(res, "Upper") <- FALSE
class(res) <- "dist"

## Organize the variables in a network
library(igraph)
g <- graph.full(ncol(dat), directed = FALSE)
{if (steps == 0) {
  V(g)$label <- datname
  V(g)$color <- datcolo
}
else {
  V(g)$label <- paste(rep(datname, times = steps + 1), rep(0:steps, each = length(datname)), sep = "_")[-2L:-9L]
  V(g)$color <- rep(datcolo, times = steps + 1)[-2L:-9L]
}}
E(g)$color <- rgb(0, 0, 0, 1 - res)
# thresh <- getThresh(dim(dat), probs = 0.05, nsim = 100L) # critical value with significance level at 0.05 
thresh <- 0.9854157114599839673019
g2 <- delete.edges(g, E(g)[res > thresh])
g2 <- delete.vertices(g2, which(degree(g2) < 1))



## Modello
dat <- data[, nomivar]
dat <- as.matrix(dat)
dat <- embed(dat, steps + 1)[, -2L:-9L]
colnames(dat) <- V(g)$label
wo <- as.logical(g2[1, ])
wo[1] <- TRUE
dat <- as.data.frame(dat[, wo])

## Fitting via lm
alf <- anova(linfit <- lm(x1_0 ~ ., data = dat))
pe <- rev(exp(-seq_len(nrow(dat)) / 288))
pe <- pe / sum(pe)
alfw <- anova(linfitW <- lm(x1_0 ~ ., data = dat, weights = pe))

# Outlier-robust penalties??
# L12 <- function(x) 2/(1+exp(-abs(x)^2))-1
# curve(L12(x), -5,5)
# L11 <- function(x) 2/(1+exp(-abs(x)))-1
# curve(L11(x), -5,5, add=TRUE, col = 2)
# L21 <- function(x) (2/(1+exp(-abs(x)))-1)^2
# curve(L21(x), -5,5, add=TRUE, col = 3)

