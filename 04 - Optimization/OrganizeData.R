rm(list = ls(all = TRUE))
graphics.off()

## Loading data
dataset <- "dati_Uff2_stratega2.csv"
# dataset <- "dati_Uff2_stratega.csv"
source("readData.R")
.missing <- c("v3", "v9", "v10", "v12")
source("missing.R")

## Data ajdustments
data$d1[data$d1 < 0] <- 0
data$d2[data$d2 < 0] <- 0
data$fc[data$fc < 0] <- 0
data$d1[data$d1 > 0] <- 100
data$d2[data$d2 > 0] <- 100
data$fc[data$fc > 1] <- 1
bl <- c(0, abs(diff(data$b))) # blinder absolute movements
wh <- data$y > 0
rt <- diff(cumsum(data$y)[unique(c(1L, which(wh)))]) / diff(seq_along(wh)[unique(c(1L, which(wh)))])
tmp <- data$y + 0
tmp[unique(c(1L, which(wh)))] <- c(rt, 0)
dyn.load("EMalgo.so")
invisible(.Call("cumSumOnZeros", tmp))
thresh <- optimize(function(x) (sum(tmp[tmp > x]) - sum(data$y))**2, c(0, .1))
tmp[tmp <= thresh$minimum] <- 0


## ENERGY CONSUMPTION - learning data and prediction set
dat <- cbind(data[-1L, c("fc", "d1", "d2")], bl = bl[-1L], data[-NROW(data), c("x1", "x2", "y_load", "v5", "v6", "v8")])
ec <- cbind(yth = tmp[-1L], dat)
names(ec)[2:5] <- paste(names(ec)[2:5], 0, sep= "_")
names(ec)[6:11] <- paste(names(ec)[6:11], 1, sep= "_")
oec <- data[NROW(data), c("b", "x1", "x2", "y_load", "v5", "v6", "v8")]
names(oec) <- paste(c("b", "x1", "x2", "y_load", "v5", "v6", "v8"), 1, sep = "_")
oec <- as.matrix(oec)


## PREDICTED MEAN VOTE - learning data and prediction set
steps <- 4
nomivar <- c("x1", paste0("v", c(1:2, 4, 7, 14:15)), "y_load", "fc")
dat <- data[, nomivar]
dat <- as.matrix(dat)
dat <- embed(dat, steps + 1)[, -2L:-8L]
colnames(dat) <- paste(rep(nomivar, times = steps + 1), rep(0:steps, each = length(nomivar)), sep = "_")[-2L:-8L]
pmv <- as.data.frame(dat)
opmv <- embed(as.matrix(tail(data[, nomivar], steps)), steps)
colnames(opmv) <- paste(rep(nomivar, times = steps), rep(1:steps, each = length(nomivar)), sep = "_")


## DAYLIGHT GLARE INDEX - learning data and prediction set
steps <- 1
nomivar <- c("x2", "v5", "v6", "v8", "v11", "d1", "d2", "b")
dat <- data[, nomivar]
dat <- as.matrix(dat)
dat <- embed(dat, steps + 1)[, -2L:-5L]
colnames(dat) <- paste(rep(nomivar, times = steps + 1), rep(0:steps, each = length(nomivar)), sep = "_")[-2L:-5L]
dgi <- as.data.frame(dat)
odgi <- embed(as.matrix(tail(data[, nomivar], steps)), steps)
colnames(odgi) <- paste(rep(nomivar, times = steps), rep(1:steps, each = length(nomivar)), sep = "_")


## Saving data
save(pmv, dgi, ec, file = "data.RData")
save(opmv, odgi, oec, file = "dnow.RData")
