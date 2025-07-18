#' @title Control parameters for the kliep function
#'
#' @export
control_kliep <- function(scale = NULL,
                          ...) {
  append(list(scale = scale),
         list(...))
}
