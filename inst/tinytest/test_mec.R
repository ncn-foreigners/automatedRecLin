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
methods_hm <- list("name" = "hit_miss",
                   "surname" = "hit_miss")
true_matches <- data.table("a" = 1:8, "b" = 1:8)
confusion_b <- matrix(c(4, 0, 4, 112), nrow = 2, ncol = 2)
rownames(confusion_b) <- c("Actual Positive", "Actual Negative")
colnames(confusion_b) <- c("Predicted Positive", "Predicted Negative")

set.seed(1)

expect_silent(
  mec(A = A_example, B = B_example, variables = variables)
)

expect_equal(
  mec(A = A_example, B = B_example, variables = variables)$M_est,
  data.table("a" = 1:4, "b" = 1:4, "ratio" = rep(720, 4))
)

expect_equal(
  mec(A = A_example, B = B_example, variables = variables)$b_params,
  data.table("variable" = c("gamma_name", "gamma_surname"),
             "theta" = c(1, 1), "eta" = c(0.04166666666666667823149, 0.03333333333333331205406))
)

expect_equal(
  mec(A = A_example, B = B_example, variables = variables, true_matches = true_matches)$eval_metrics,
  c("FLR" = 0, "MMR" = 0.5)
)

expect_equal(
  mec(A = A_example, B = B_example, variables = variables, true_matches = true_matches)$confusion,
  confusion_b
)

expect_equal(
  mec(A = A_example, B = B_example, variables = variables, methods = methods_hm)$M_est,
  data.table("a" = 1:4, "b" = 1:4, "ratio" = rep(573.0984617692614619955, 4))
)

expect_equal(
  mec(A = A_example, B = B_example, variables = variables, methods = methods_hm)$hm_params,
  data.table("variable" = c("gamma_name", "gamma_surname"),
             "theta" = c(1, 1), "eta" = c(0.04616298284003408219922, 0.03847198038861038282832))
)

expect_equal(
  mec(A = A_example, B = B_example, variables = variables,
      comparators = comparators, methods = methods_cpar)$M_est[, c("a", "b")],
  data.table("a" = c(7, 6, 8, 5, 1, 2, 3, 4), "b" = c(7, 6, 8, 5, 1, 2, 3, 4))
)

expect_equal(
  mec(A = A_example, B = B_example, variables = variables,
      comparators = comparators, methods = methods_cpar)$cpar_params,
  data.table("variable" = c("gamma_name", "gamma_surname"),
             "p_0_M" = c(0.625, 0.5),
             "alpha_M" = c(138.4622794465536514963, 120.6657058981201657843),
             "beta_M" = c(2199.106791209967013856, 1974.529732878328104562),
             "p_0_U" = c(0.04166666666666666435370, 0.03333333333333333287074),
             "alpha_U" = c(6.516735740295027667912, 4.622775398311523176176),
             "beta_U" = c(11.173089162681424824086, 7.167260899080440061937))
)

expect_warning(
  mec(A = A_example, B = B_example, variables = variables,
      comparators = comparators, methods = methods_cnonpar)
)

expect_equal(
  mec(A = A_example, B = B_example, variables = variables,
      comparators = comparators, methods = methods_cnonpar)$M_est[, c("a", "b")],
  data.table("a" = 1:4, "b" = 1:4)
)
