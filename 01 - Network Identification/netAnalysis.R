rm(list = ls(all = TRUE))
graphics.off()

## Network analysis
source("Entropy.R")

dataset <- "dati_Uff2_stratega2.csv"
source("readData.R")
.missing <- c("v3", "v9", "v10", "v12")
source("missing.R")

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

res <- combn(ncol(dat), 2, function(xx) mi_qv(dat[, xx[1]], dat[, xx[2]]), simplify = TRUE)
attr(res, "Size") <- ncol(dat)
attr(res, "call") <- ""
attr(res, "Diag") <- FALSE
attr(res, "Upper") <- FALSE
class(res) <- "dist"

thresh <- 0.9854157114599839673019 # 5% quantile obtained by simulation of i.i.d. normals (n.var. and n.oss. are the same of this dataset)
#thresh <- 0.9851305041418656482932 # minimum of the same simulations
cat("Threshold: ", thresh, "\n", sep = "")
cat("Ratio of connections: ", round(sum(res < thresh) / choose(ncol(dat), 2) * 100, 2), "%\n", sep = "")


library(igraph)
g <- graph.full(ncol(dat), directed = FALSE)
{if (steps == 0) {
  V(g)$label <- datname
  V(g)$color <- datcolo
}
else {
  V(g)$label <- paste(rep(datname, times = steps + 1), rep(0:steps, each = length(datname)), sep="_")
  V(g)$color <- rep(datcolo, times = steps + 1)
}}
E(g)$color <- rgb(0, 0, 0, 1 - res)
g2 <- delete.edges(g, E(g)[res > thresh])
g2 <- delete.vertices(g2, which(degree(g2) < 1))


## Plotting the networks
# coords <- layout.circle(g)
plot(g, vertex.label.color= "black")#, layout= coords)
plot(g, edge.color = rgb(0, 0, 0, res <= thresh), vertex.label.color= "black")#, layout= coords)

# coords <- layout.circle(g2)
plot(g2, edge.color = rgb(0.5, 0.5, 0.5), vertex.label.color= "black")#, layout = coords)

## Extracting the nodes names related to energy consumption at time $t = 0$
# V(g2)$label[ g2[[which(V(g2)$label == "y_0"), ]][[1]] ]
