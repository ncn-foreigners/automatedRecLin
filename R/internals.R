binary_formula <- function(df) {
  par <- apply(df, 2, mean)
  par
}

p_0_formula <- function(df) {
  par <- apply(df, 2, function(x) {
    x <- x[x = 0]
    length(x) / NROW(df)
  })
}

f_alpha <- function(alpha, gamma) {
  gamma <- gamma[gamma > 0]
  gamma_mean <- mean(gamma)
  sum(log(gamma) - log(gamma_mean) - digamma(alpha) + log(alpha))
}

gamma_plus_formula <- function(df) {
  par <- apply(df, 2, function(x) {
    x <- x[x > 0]
    mean(x)
  })
}


alpha_formula <- function(df, fun) {
  par <- apply(df, 2, function(y) {
    fun(x = 1, fn = f_alpha, gamma = y, method = "Newton")$x
  })
  par
}
