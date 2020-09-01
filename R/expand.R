#' expand
#'
#' This function is used to expand the Y(category) to a indicator vector
#'
#' @param y the current state
#' @param ml the length of the indicator vector
#'
#' @return a indicator vector
#'
#' @export
#'
#' @examples
#'
#'
#'
#'
#'

expand=function(y,ml){    #expand it into indicate variable
  ry=rep(0,ml)
  ry[y-1]=1
  return(ry)
}
