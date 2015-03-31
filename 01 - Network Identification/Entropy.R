## Variable quantization
quantize <- function(x, length.out = 100) {
  # Function for quantization of variable x
  # number of quantization levels will be "length.out + 2"
  if(is.factor(x) || is.character(x)) return(factor(x))
  if(is.numeric(x)) {
    if (diff(range(x)) != 0) {
      cut(x, c(-Inf, seq(min(x), max(x), length.out = length.out), Inf))
    }
    else {
      cut(x, c(-Inf, x[1], Inf))
    }
  }
  else{
    stop("argument \"x\" is not a factor nor a numeric object")
  }
}

bicut <- function(x) cut(x, c(-Inf, mean(x), Inf))

## Dissimilarity based on mutual entropy of quantized variables
mi_qv <- function(q1, q2) {
  q1 <- quantize(q1)
  q2 <- quantize(q2)
  # Joint probability distribution
  pj <- table(q1, q2)
  pj <- pj / sum(pj)
  # Independent probability distribution
  pi <- outer(table(q1), table(q2))
  pi <- pi/ sum(pi)
  # Entropies calculations
  ej <- -sum(pj * log(pj), na.rm = TRUE)
  ei <- -sum(pj * log(pi), na.rm = TRUE)
#   me <- log(nlevels(q1) * nlevels(q2))
  # Entropic dependence dissimilarity
  dst <- 2 - ei / ej
  return(ifelse(is.na(dst), 1, dst))
}

## Compoutation of the critical value via sumulations
getThresh <- function(dim, probs = 0.05, nsim = 100L) {
  stats <- replicate(nsim, {
    AA <- array(runif(prod(dim)), dim)
    min(combn(ncol(AA), 2, function(xx) mi_qv(AA[, xx[1]], AA[, xx[2]]), simplify = TRUE))
    }, simplify = "array")
  cv <- quantile(stats, probs)
  return(cv)
}

