#' @title Controls for the `kliep` Function
#'
#' @author Adam Struzik
#'
#' @description
#' Controls for the \link[densityratio]{kliep} function used in the package.
#'
#' @param scale `"numerator"`, `"denominator"` or `NULL`, indicating whether to standardize each numeric variable according
#' to the numerator means and standard deviations, the denominator means and standard deviations,
#' or apply no standardization at all.
#' @param progressbar Logical indicating whether or not to display a progressbar.
#' @param nfold Number of cross-validation folds used in order to calculate the optimal sigma value (default is 2-fold cv).
#' @param ... Additional arguments.
#'
#' @return
#' Returns a list with parameters.
#'
#' @export
control_kliep <- function(scale = NULL,
                          progressbar = FALSE,
                          nfold = 2,
                          ...) {
  append(list(scale = scale, progressbar = progressbar, nfold = nfold),
         list(...))
}
