#' norm
#'
#' This function is used to normalize a vector to have unit length
#'
#' @param x a numeric vector
#'
#' @return a normalized vector with length 1
#'
#'
#'
#'
#'
#'

norm = function(x){
  y=x/sqrt(sum(x^2))
  return(y)
}
