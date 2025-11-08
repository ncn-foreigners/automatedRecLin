library(automatedRecLin)
library(data.table)
library(xgboost)

df_1 <- data.frame(
  "name" = c("James", "Emma", "William", "Olivia", "Thomas",
             "Sophie", "Harry", "Amelia", "George", "Isabella"),
  "surname" = c("Smith", "Johnson", "Brown", "Taylor", "Wilson",
                "Davis", "Clark", "Harris", "Lewis", "Walker")
)
df_2 <- data.frame(
  "name" = c("James", "Ema", "Wimliam", "Olivia", "Charlotte",
             "Henry", "Lucy", "Edward", "Alice", "Jack"),
  "surname" = c("Smith", "Johnson", "Bron", "Tailor", "Moore",
                "Evans", "Hall", "Wright", "Green", "King")
)
matches <- data.frame("a" = 1:4, "b" = 1:4)
df_new_1 <- data.frame(
  "name" = c("John", "Emily", "Mark", "Anna", "David"),
  "surname" = c("Smith", "Johnson", "Taylor", "Williams", "Brown")
)
df_new_2 <- data.frame(
  "name" = c("John", "Emely", "Mark", "Michael"),
  "surname" = c("Smitth", "Johnson", "Tailor", "Henders")
)
true_matches <- data.frame("a" = 1:3, "b" = 1:3)

comparators = list("name" = jarowinkler_complement(),
                   "surname" = jarowinkler_complement())
methods_cpar <- list("name" = "continuous_parametric",
                     "surname" = "continuous_parametric")
methods_cnonpar <- list("name" = "continuous_nonparametric",
                        "surname" = "continuous_nonparametric")

confusion_cpar <- matrix(c(2, 0, 1, 17), nrow = 2, ncol = 2)
rownames(confusion_cpar) <- c("Actual Positive", "Actual Negative")
colnames(confusion_cpar) <- c("Predicted Positive", "Predicted Negative")

expect_silent(
  train_rec_lin(A = df_1, B = df_2,
                variables = c("name", "surname"),
                matches = matches)
)

expect_warning(
  train_rec_lin(A = df_1, B = df_2,
                variables = c("name", "surname"),
                matches = matches,
                prob_ratio = "2")
)

expect_equal(
  train_rec_lin(A = df_1, B = df_2,
                variables = c("name", "surname"),
                matches = matches)$b_params,
  data.table("variable" = c("gamma_name", "gamma_surname"),
             "theta" = c(0.5, 0.5),
             "eta" = c(0.02, 0.02))
)

expect_equal(
  train_rec_lin(A = df_1, B = df_2,
                variables = c("name", "surname"),
                matches = matches)$pi_est,
  0.04
)

expect_equal(
  train_rec_lin(A = df_1, B = df_2,
                variables = c("name", "surname"),
                matches = matches)$match_prop,
  0.4
)

expect_silent(
  train_rec_lin(A = df_1, B = df_2,
                variables = c("name", "surname"),
                matches = matches,
                comparators = comparators,
                methods = methods_cpar)
)

expect_warning(
  train_rec_lin(A = df_1, B = df_2,
                variables = c("name", "surname"),
                matches = matches,
                comparators = comparators,
                methods = methods_cpar,
                prob_ratio = "2")
)

expect_equal(
  train_rec_lin(A = df_1, B = df_2,
                variables = c("name", "surname"),
                matches = matches,
                comparators = comparators,
                methods = methods_cpar)$cpar_params,
  data.table("variable" = c("gamma_name", "gamma_surname"),
             "p_0_M" = c(0.5, 0.5),
             "p_0_Omega" = c(0.02000000000000000041633, 0.02000000000000000041633),
             "alpha_M" = c(224.66614858286149569722, 15.65946723172959487158),
             "alpha_Omega" = c(5.407152306652986517577, 6.105019794944631428280),
             "beta_M" = c(2516.2608641280485244351, 176.1690063569579081104),
             "beta_Omega" = c(7.975415092480718115553, 9.363290038376906210260))
)

expect_warning(
  train_rec_lin(A = df_1, B = df_2,
                variables = c("name", "surname"),
                matches = matches,
                comparators = comparators,
                methods = methods_cnonpar)
)

expect_equal(
  train_rec_lin(A = df_1, B = df_2,
                variables = c("name", "surname"),
                matches = matches,
                comparators = comparators,
                methods = methods_cnonpar)$cnonpar_params,
  data.table("variable" = c("gamma_name", "gamma_surname"),
             "p_0_M_cnonpar" = c(0.5, 0.5),
             "p_0_U_cnonpar" = c(0.02, 0.02))
)

model_b <- train_rec_lin(A = df_1, B = df_2,
                         variables = c("name", "surname"),
                         matches = matches)

expect_equal(
  predict(model_b, df_new_1, df_new_2)$M_est,
  data.table("a" = 1:3, "b" = 1:3, "ratio" = rep(12.75510204081633602868, 3))
)

model_cpar <- train_rec_lin(A = df_1, B = df_2,
                            variables = c("name", "surname"),
                            matches = matches,
                            comparators = comparators,
                            methods = methods_cpar)

expect_equal(
  predict(model_cpar, df_new_1, df_new_2)$M_est,
  data.table("a" = c(1, 3), "b" = c(1, 3),
             "ratio" = c(59126.609831455622042995, 4201.052699951693284675))
)

expect_equal(
  predict(model_cpar, df_new_1, df_new_2, true_matches = true_matches)$eval_metrics,
  c("FLR" = 0, "MMR" = 1/3)
)

expect_equal(
  predict(model_cpar, df_new_1, df_new_2, true_matches = true_matches)$confusion,
  confusion_cpar
)

model_cnonpar <- train_rec_lin(A = df_1, B = df_2,
                               variables = c("name", "surname"),
                               matches = matches,
                               comparators = comparators,
                               methods = methods_cnonpar)

expect_equal(
  predict(model_cnonpar, df_new_1, df_new_2)$M_est,
  data.table("a" = c(2, 3, 1), "b" = c(2, 3, 1),
             "ratio" = c(45.72851745592178218658, 28.13814113792743043518, 25.60155559185211160411))
)

# custom model

vectors <- comparison_vectors(A = df_1, B = df_2, variables = c("name", "surname"),
                              comparators = comparators, matches = matches)
train_data <- xgb.DMatrix(
  data = as.matrix(vectors$Omega[, c("gamma_name", "gamma_surname")]),
  label = vectors$Omega$match
)
params <- list(objective = "binary:logistic",
               eval_metric = "logloss")
model_xgb <- xgboost::xgboost(data = train_data, params = params,
                              nrounds = 100, verbose = 0)

expect_silent(
  custom_rec_lin_model(model_xgb, vectors)
)

expect_equal(
  custom_rec_lin_model(model_xgb, vectors)$pi_est,
  0.04
)

expect_equal(
  custom_rec_lin_model(model_xgb, vectors)$match_prop,
  0.4
)

custom_xgb_model <- custom_rec_lin_model(model_xgb, vectors)

expect_equal(
  predict(custom_xgb_model, df_new_1, df_new_2,
          data_type = "matrix", type = "prob")$M_est,
  data.table("a" = 1:3, "b" = 1:3,
             "ratio" = rep(37.02279418163015378695, 3))
)

expect_equal(
  predict(custom_xgb_model, df_new_1, df_new_2,
          data_type = "matrix", type = "prob")$flr_est,
  0.13274158431759
)

expect_equal(
  predict(custom_xgb_model, df_new_1, df_new_2,
          data_type = "matrix", type = "prob")$mmr_est,
  0.13274158431759
)
