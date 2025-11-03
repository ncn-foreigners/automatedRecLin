library(automatedRecLin)
library(data.table)

data("A_example")
data("B_example")

variables <- c("name", "surname")
comparators <- list("name" = jarowinkler_complement(),
                    "surname" = jarowinkler_complement())
methods_cpar <- list("name" = "continuous_parametric",
                     "surname" = "continuous_parametric")
methods_cnonpar <- list("name" = "continuous_nonparametric",
                        "surname" = "continuous_nonparametric")


expect_silent(
  mec(A = A_example, B = B_example, variables = variables)
)

expect_equal(
  mec(A = A_example, B = B_example, variables = variables)$M_est,
  data.table("a" = 1:4, "b" = 1:4, "ratio" = rep(720, 4))
)

expect_equal(
  mec(A = A_example, B = B_example, variables = variables,
      comparators = comparators, methods = methods_cpar)$M_est[, c("a", "b")],
  data.table("a" = c(7, 6, 8, 5, 1, 2, 3, 4), "b" = c(7, 6, 8, 5, 1, 2, 3, 4))
)

expect_warning(
  mec(A = A_example, B = B_example, variables = variables,
      comparators = comparators, methods = methods_cnonpar)
)
