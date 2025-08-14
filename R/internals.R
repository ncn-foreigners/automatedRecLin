binary_formula <- function(df) {
  par <- apply(df, 2, mean)
  par
}

p_0_formula <- function(df) {
  par <- apply(df, 2, function(x) {
    x <- x[x == 0]
    length(x) / NROW(df)
  })
  par
}

f_alpha <- function(alpha, gamma) {
  gamma <- gamma[gamma > 0]
  gamma_mean <- mean(gamma)
  sum(log(gamma) - log(gamma_mean) - digamma(alpha) + log(alpha))
}

f_alpha_iterative <- function(alpha, beta, gamma) {
  gamma <- gamma[gamma > 0]
  sum(log(gamma) + log(beta) - digamma(alpha))
}

gamma_plus_formula <- function(df) {
  par <- apply(df, 2, function(x) {
    x <- x[x > 0]
    mean(x)
  })
  par
}

alpha_formula <- function(df, fun) {
  par <- apply(df, 2, function(y) {
    fun(x = 1, fn = f_alpha, gamma = y, method = "Newton")$x
  })
  par
}

alpha_formula_iterative <- function(df, fun, beta) {
  K <- length(beta)
  par <- sapply(1:K, function(x) {
    fun(x = 1, fn = f_alpha_iterative, gamma = df[, as.numeric(x)], beta = beta[x], method = "Newton")$x
  })
  par
}

#' @importFrom stats dgamma
hurdle_gamma_density <- function(x, p_0, alpha, beta) {
  ifelse(x == 0, p_0, 1) *
    ifelse(x > 0, (1 - p_0) * stats::dgamma(x = x, shape = alpha, scale = 1 / beta), 1)
}

fixed_n_M <- function(n, ratio_gamma) {
  function(n_M) {
    # sum(n_M * ratio_gamma / (n_M * (ratio_gamma - 1) + n))
    sum(pmin(n_M * ratio_gamma / (n_M * (ratio_gamma - 1) + n), 1))
  }
}
