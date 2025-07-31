#' @title Controls for the kliep Function
#'
#' @author Adam Struzik
#'
#' @description
#' Controls for the \link[densityratio]{kliep} function used in the package.
#'
#' @param scale `"numerator"`, `"denominator"` or `NULL`, indicating whether to standardize each numeric variable according
#' to the numerator means and standard deviations, the denominator means and standard deviations,
#' or apply no standardization at all.
#' @param ... Additional arguments.
#'
#' @return
#' Returns a list with parameters.
#'
#' @export
control_kliep <- function(scale = NULL,
                          ...) {
  append(list(scale = scale),
         list(...))
}
