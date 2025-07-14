#' @title Calculate absolute distance
#'
#' @export
abs_distance <- function() {
  function(x, y) {
    abs(x - y)
  }
}
