rm(list = ls(all = TRUE))
graphics.off()

## Network analysis
source("Entropy.R")

dataset <- "dati_Uff2_stratega2.csv"
source("readData.R")
.missing <- c("v3", "v9", "v10", "v12")
source("missing.R")

library(igraph)

## Data ajdustments
data$d1[data$d1 < 0] <- 0
data$d2[data$d2 < 0] <- 0
data$fc[data$fc < 0] <- 0
data$d1[data$d1 > 100] <- 100
data$d2[data$d2 > 100] <- 100
data$fc[data$fc > 1] <- 1

## Paramters for the analysis
steps <- 4

## Analysis of dependences
dat <- data[, -1:-7]
datname <- names(data[, -1:-7])
datcolo <- c(rep("#00ff00", 4), rep("#ffcc00", 8), rep("#00ccff", 4), rep("#ffcc00", 3))
dat <- apply(dat, 2, quantize)
dat <- embed(dat, steps + 1)

idx <- matrix(seq_len(nrow(dat) %/% 8L * 8L), nrow(dat) %/% 8L, 8L)

#system.time(thresh <- getThresh(c(nrow(dat) %/% 8L, ncol(dat)), nsim = 100L))
#thresh <- 0.8458930618700987169589 # 5% with steps = 0
thresh <- 0.844962627685229894503 # 5% with steps = 4

g <- list()
g2 <- list()
for (i in seq_len(8L)) {
  res <- combn(ncol(dat), 2, function(xx) mi_qv(dat[idx[, i], xx[1]], dat[idx[, i], xx[2]]), simplify = TRUE)
  res[res < 0] <- 0 
  attr(res, "Size") <- ncol(dat)
  attr(res, "call") <- ""
  attr(res, "Diag") <- FALSE
  attr(res, "Upper") <- FALSE
  class(res) <- "dist"

  # library(igraph)
  g[[i]] <- graph.full(ncol(dat), directed = FALSE)
  {if (steps == 0) {
    V(g[[i]])$label <- datname
    V(g[[i]])$color <- datcolo
  }
  else {
    V(g[[i]])$label <- paste(rep(datname, times = steps + 1), rep(0:steps, each = length(datname)), sep="_")
    V(g[[i]])$color <- rep(datcolo, times = steps + 1)
  }}
  E(g[[i]])$color <- rgb(0, 0, 0, 1 - res)
  E(g[[i]])$mi <- res ## MUTUAL INFORMATION for the test
  g2[[i]] <- delete.edges(g[[i]], E(g[[i]])[res > thresh])
  g2[[i]] <- delete.vertices(g2[[i]], which(degree(g2[[i]]) < 1))
}

{if (steps < .5) {
  tb1 <- sapply(g2, function(x) V(x)$label[unlist(x[[which(V(x)$label == "x1"), ]])])
  tb1 <- table(unlist(tb1)) / 8L
  tb2 <- sapply(g2, function(x) V(x)$label[unlist(x[[which(V(x)$label == "x2"), ]])])
  tb2 <- table(unlist(tb2)) / 8L
} else {
  tb1 <- sapply(g2, function(x) V(x)$label[unlist(x[[which(V(x)$label == "x1_0"), ]])])
  tb1 <- table(unlist(tb1)) / 8L
  tb2 <- sapply(g2, function(x) V(x)$label[unlist(x[[which(V(x)$label == "x2_0"), ]])])
  tb2 <- table(unlist(tb2)) / 8L
}}


