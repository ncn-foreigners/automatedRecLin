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
#' @export
abs_distance <- function() {
  function(x, y) {
    abs(x - y)
  }
}

#' @importFrom reclin2 cmp_jarowinkler
#'
#' @title Jaro-Winkler Distance Complement
#'
#' @author Adam Struzik
#'
#' @description
#' Creates a function that calculates the complement of the Jaro-Winkler
#' distance between two strings (i.e.,
#' \eqn{1 - \text{Jaro-Winkler distance}}).
#'
#' @return
#' Returns a function taking two string arguments, `x` and `y`,
#' and returning the complement of the Jaro-Winkler distance.
#'
#' @export
jarowinkler_complement <- function() {
  temp <- reclin2::cmp_jarowinkler()
  function(x, y) {
    1 - temp(x, y)
  }
}
