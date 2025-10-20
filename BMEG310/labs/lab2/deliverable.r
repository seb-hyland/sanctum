# First, define helper function to sum the elements of a vector
sumFun <- function(x) {
  sum <- 0
  for (item in x) {
    sum <- sum + item
  }
  sum
}

varFun <- function(x) {
  len <- length(x)
  mean <- sumFun(x) / len

  square_error <- sumFun((x - mean)**2)
  variance <- (1 / (len - 1)) * square_error

  print(variance)
}

varFun(6:36)
