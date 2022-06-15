## R Script for controlling missing data
# This script require a dataset called "data" and
# a vector ".missing" with the names of the variables

## user missing varibles will be removed
dataNames <- names(data)
notMissing <- !(dataNames %in% .missing)
data <- data[, notMissing]

## constant variables will be shown to the user
posConst <- c()
for (i in seq_len(dim(data)[2])) {
  if (var(data[, i]) == 0) posConst <- c(posConst, i)
}

cat("Constant Variables in the dataset:\n", paste0("[", names(data)[posConst], "]"), "\n")



# posVariab <- seq_len(dim(data)[2])[-posConst]

# par(ask = TRUE)
# for (i in posVariab) {
#   plot(data[, i], main = names(data)[i])
# }

# pairs(data[, 12:18], pch = ".")
# plot(data$v12, type = "l")
