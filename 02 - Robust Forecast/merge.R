## Merge the robust results ...
load("erf1o.RData")
load("erf1e.RData")
load("erf1n.RData")
load("erf2o.RData")
load("erf2e.RData")
load("erf2n.RData")

## ... in a unique file.
save.image("errbustforex.RData")

rm(list = ls(all = TRUE))

## Merge the ARIMAX results ...
load("arimax1.RData")
load("arimax2.RData")

## ... in a unique file.
save.image("errARIMAX.RData")

