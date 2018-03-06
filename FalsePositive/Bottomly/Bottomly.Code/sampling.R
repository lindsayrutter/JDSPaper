
sub.sample <- function(count.data, n1=3, n2=3) {
  number <- ncol(count.data)
  n <- n1 + n2
  sampling <- sample(number, n)
  subsample <- count.data[, sampling]
  return(subsample)
}

sampling.groups <- function (sub.group1, sub.group2, n1=3, n2=3) {
  number1 <- ncol(sub.group1)
  number2 <- ncol(sub.group2)
  sampling1 <- sample(number1, n1)
  sampling2 <- sample(number2, n2)
  subsample <- cbind(sub.group1[, sampling1], sub.group2[, sampling2])
  return(subsample)
}


    