#' @title Absolute Distance Comparison Function
#'
#' @author Adam Struzik
#'
#' @description
#' Creates a function that calculates the absolute distance between two values.
#'
#' @return
#' Returns a function taking two arguments, `x` and `y`, and returning their absolute difference.
#'
#' @examples
#' cmp <- abs_distance()
#' cmp(1, 5) # returns 4
#'
#' @export
abs_distance <- function() {
  function(x, y) {
    abs(x - y)
  }
}
