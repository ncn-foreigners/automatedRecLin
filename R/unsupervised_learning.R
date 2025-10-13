#' @import data.table
#' @importFrom stats runif
#' @importFrom stats predict
#' @importFrom purrr partial
#' @importFrom utils head
#'
#' @title Unsupervised Maximum Entropy Classifier for Record Linkage
#'
#' @author Adam Struzik
#'
#' @description
#' Implements several extensions to the maximum entropy classification (MEC) algorithm for record linkage
#' (see [Lee et al. (2022)](https://www150.statcan.gc.ca/n1/pub/12-001-x/2022001/article/00007-eng.htm)),
#' iteratively estimating probability/density ratios to classify record pairs into matches and non-matches
#' based on comparison vectors.
#'
#' @param A A duplicate-free `data.frame` or `data.table`.
#' @param B A duplicate-free `data.frame` or `data.table`.
#' @param variables A character vector of key variables used to create comparison vectors.
#' @param comparators A named list of functions for comparing pairs of records.
#' @param methods A named list of methods used for estimation (`"binary"`, `"continuous_parametric"`, `"continuous_nonparametric"` or `"hit_miss"`).
#' @param duplicates_in_A Logical indicating whether to allow `A` to have duplicate records.
#' @param start_params Start parameters for the `"binary"`, `"continuous_parametric"` and `"hit_miss"` methods.
#' @param nonpar_hurdle Logical indicating whether to use a hurdle model or not
#' (used only if the `"continuous_nonparametric"` method has been chosen for at least one variable).
#' @param set_construction A method for constructing the predicted set of matches (`"size"` or `"flr"`).
#' @param target_flr A target false link rate (FLR) (used only if `set_construction == "flr"`).
#' @param max_iter_bisection A maximum number of iterations for the bisection procedure (used only if `set_construction == "flr"`).
#' @param tol Error tolerance in the bisection procedure (used only if `set_construction == "flr"`).
#' @param delta A numeric value specifying the tolerance for the change in the estimated number of matches between iterations.
#' @param eps A numeric value specifying the tolerance for the change in model parameters between iterations.
#' @param max_iter_em A maximum number of iterations for the EM algorithm
#' (used only if the `"hit_miss"` method has been chosen for at least one variable).
#' @param tol_em Error tolerance in the EM algorithm
#' (used only if the `"hit_miss"` method has been chosen for at least one variable).
#' @param controls_nleqslv Controls passed to the \link[nleqslv]{nleqslv} function
#' (only if the `"continuous_parametric"` method has been chosen for at least one variable).
#' @param controls_kliep Controls passed to the \link[densityratio]{kliep} function
#' (only if the `"continuous_nonparametric"` method has been chosen for at least one variable).
#' @param true_matches A `data.frame` or `data.table` indicating known matches.
#'
#' @details
#' Consider two datasets without duplicates: \eqn{A} and \eqn{B}.
#' Let the bipartite comparison space \eqn{\Omega = A \times B} consist of
#' matches \eqn{M} and non-matches \eqn{U} between the records in files
#' \eqn{A} and \eqn{B}. For any pair of records \eqn{(a,b) \in \Omega},
#' let \eqn{\pmb{\gamma}_{ab} = (\gamma_{ab}^1,\gamma_{ab}^2,
#' \ldots,\gamma_{ab}^K)'} be the comparison vector between
#' a set of key variables. The original MEC algorithm uses the binary
#' comparison function to evaluate record pairs across two datasets.
#' However, this approach may be insufficient when handling datasets
#' with frequent errors across multiple variables.
#'
#' We propose the use of continuous comparison functions to address
#' the limitations of binary comparison methods. We consider every
#' semi-metric, i.e., a function \eqn{d: A \times B \to \mathbb{R}},
#' satisfying the following conditions:\cr
#' \enumerate{
#' \item{\eqn{d(x,y) \geq 0},}
#' \item{\eqn{d(x,y) = 0} if and only if \eqn{x = y},}
#' \item{\eqn{d(x,y) = d(y,x)}.}
#' }
#' For example, we can use \eqn{1 - \text{Jaro-Winkler distance}} for character variables
#' (which is implemented in the `automatedRecLin` package as the `jarowinkler_complement` function)
#' or the Euclidean distance for numerical variables. The `automatedRecLin` package allows the use of
#' a different comparison function for each key variable (which should be specified
#' as a list in the `comparators` argument). The default function
#' for each key variable is \link[reclin2]{cmp_identical}
#' (the binary comparison function).
#'
#' The `mec` function offers different approaches to estimate the
#' probability/density ratio between matches and non-matches,
#' which should be specified as a list in the `methods` argument.
#' The available methods suitable for the binary comparison function
#' are `"binary"` and `"hit_miss"`. Both assume that \eqn{\gamma_{ab}^k|M}
#' and \eqn{\gamma_{ab}^k|U} follow Bernoulli distributions.
#' `"binary"` and `"hit_miss"` both estimate the parameters for the matches iteratively,
#' but `"binary"` estimates the parameters for the non-matches
#' only at the start, while `"hit_miss"` does
#' so iteratively using a hit-miss model (for details see
#' [Lee et al. (2022)](https://www150.statcan.gc.ca/n1/pub/12-001-x/2022001/article/00007-eng.htm)).
#' `"binary"` is the default method for each variable.
#'
#' For the continuous semi-metrics we suggest the usage
#' of `"continuous_parametric"` or `"continuous_nonparametric"`
#' method. The `"continuous_parametric"` method assumes that
#' \eqn{\gamma_{ab}^k|M} and \eqn{\gamma_{ab}^k|U} follow
#' hurdle Gamma distributions. The density function of a hurdle
#' Gamma distribution is characterized by three parameters
#' \eqn{p_0 \in (0,1)} and \eqn{\alpha, \beta > 0} as follows:
#' \deqn{
#' f(x;p_0,\alpha,\beta) = p_0^{\mathbb{I}(x = 0)}[(1 - p_0) v(x;\alpha,\beta)]^{\mathbb{I}(x > 0)},
#' }
#' where
#' \deqn{
#' v(x;\alpha,\beta) = \frac{\beta^{\alpha} x^{\alpha - 1} \exp(-\beta x)}
#' {\Gamma(\alpha)}
#' }
#' is the density function of a Gamma distribution
#' (for details see [Vo et al. (2023)](https://ideas.repec.org/a/eee/csdana/v179y2023ics0167947322002365.html)).
#' At the beginning, the algorithm estimates the parameters for the non-matches
#' and then does it iteratively for the matches.
#' The `"continuous_nonparametric"` method does not assume anything about
#' the distributions of the comparison vectors. It iteratively directly
#' estimates the density ratio between the matches and the non-matches, using
#' the Kullback-Leibler Importance Estimation Procedure (KLIEP).
#' For details see [Sugiyama et al. (2008)](https://link.springer.com/article/10.1007/s10463-008-0197-x).
#'
#' The `mec` function allows the construction of the predicted set
#' of matches using its estimated size or the bisection procedure,
#' described in [Lee et al. (2022)](https://www150.statcan.gc.ca/n1/pub/12-001-x/2022001/article/00007-eng.htm),
#' based on a target False Link Rate (FLR). To use the second option, set `set_construction = "flr"` and
#' specify a target FLR using the `target_flr` argument.
#'
#' The assumption that \eqn{A} and \eqn{B} contain no duplicate records
#' might be relaxed by allowing \eqn{A} to have duplicates. To do so,
#' set `"duplicates_in_A = TRUE`.
#'
#' @return
#' Returns a list containing:\cr
#' \itemize{
#' \item{`M_est` -- a `data.table` with predicted matches,}
#' \item{`n_M_est` -- estimated classification set size,}
#' \item{`flr_est` -- estimated false link rate (FLR),}
#' \item{`mmr_est` -- estimated missing match rate (MMR),}
#' \item{`iter_bisection` -- the number of iterations in the bisection procedure,}
#' \item{`b_vars` -- a character vector of variables used for the `"binary"` method (with the prefix `"gamma_"`),}
#' \item{`cpar_vars` -- a character vector of variables used for the `"continuous_parametric"` method (with the prefix `"gamma_"`),}
#' \item{`cnonpar_vars` -- a character vector of variables used for the `"continuous_nonparametric"` method (with the prefix `"gamma_"`),}
#' \item{`hm_vars` -- a character vector of variables used for the `"hit_miss"` method (with the prefix `"gamma_"`),}
#' \item{`b_params` -- parameters estimated using the `"binary"` method,}
#' \item{`cpar_params` -- parameters estimated using the `"continuous_parametric"` method,}
#' \item{`hm_params` -- parameters estimated using the `"hit_miss"` method,}
#' \item{`ratio_kliep` -- a result of the \link[densityratio]{kliep} function,}
#' \item{`variables` -- a character vector of key variables used for comparison,}
#' \item{`set_construction` -- a method for constructing the predicted set of matches,}
#' \item{`eval_metrics` -- standard metrics for quality assessment (if `true_matches` is provided),}
#' \item{`confusion` -- confusion matrix (if `true_matches` is provided).}
#' }
#'
#' @references
#' Lee, D., Zhang, L.-C. and Kim, J. K. (2022). Maximum entropy classification for record linkage.
#' Survey Methodology, Statistics Canada, Catalogue No. 12-001-X, Vol. 48, No. 1.
#'
#' Vo, T. H., Chauvet, G., Happe, A., Oger, E., Paquelet, S., and Garès, V. (2023).
#' Extending the Fellegi-Sunter record linkage model for mixed-type data with application to the French national health data system.
#' Computational Statistics & Data Analysis, 179, 107656.
#'
#' Sugiyama, M., Suzuki, T., Nakajima, S. et al. Direct importance estimation for covariate shift adaptation.
#' Ann Inst Stat Math 60, 699–746 (2008). \doi{10.1007/s10463-008-0197-x}
#'
#' @examples
#' df_1 <- data.frame(
#'   name = c("Emma", "Liam", "Olivia", "Noah", "Ava",
#'            "Ethan", "Sophia", "Mason", "Isabella", "James"),
#'   surname = c("Smith", "Johnson", "Williams", "Brown", "Jones",
#'               "Garcia", "Miller", "Davis", "Rodriguez", "Wilson"),
#'   city = c("New York", "Los Angeles", "Chicago", "Houston", "Phoenix",
#'            "Philadelphia", "San Antonio", "San Diego", "Dallas", "San Jose")
#' )
#'
#' df_2 <- data.frame(
#'   name = c(
#'     "Emma", "Liam", "Olivia", "Noah",
#'     "Ava", "Ehtan", "Sopia", "Mson",
#'     "Charlotte", "Benjamin", "Amelia", "Lucas"
#'   ),
#'   surname = c(
#'      "Smith", "Johnson", "Williams", "Brown",
#'     "Jnes", "Garca", "Miler", "Dvis",
#'     "Martinez", "Lee", "Hernandez", "Clark"
#'   ),
#'   city = c(
#'     "New York", "Los Angeles", "Chicago", "Houston",
#'     "Phonix", "Philadelpia", "San Antnio", "San Dieg",
#'     "Seattle", "Miami", "Boston", "Denver"
#'   )
#' )
#' true_matches <- data.frame(
#'   "a" = 1:8,
#'   "b" = 1:8
#' )
#'
#' variables <- c("name", "surname", "city")
#' comparators <- list(
#'   "name" = jarowinkler_complement(),
#'   "surname" = jarowinkler_complement(),
#'   "city" = jarowinkler_complement()
#' )
#' methods <- list(
#'   "name" = "continuous_parametric",
#'   "surname" = "continuous_parametric",
#'   "city" = "continuous_parametric"
#' )
#'
#' set.seed(1)
#' result <- mec(A = df_1, B = df_2,
#'               variables = variables,
#'               comparators = comparators,
#'               methods = methods,
#'               true_matches = true_matches)
#' result
#' @export
mec <- function(A,
                B,
                variables,
                comparators = NULL,
                methods = NULL,
                duplicates_in_A = FALSE,
                start_params = NULL,
                nonpar_hurdle = FALSE,
                set_construction = NULL,
                target_flr = 0.03,
                max_iter_bisection = 100,
                tol = 0.005,
                delta = 0.5,
                eps = 0.05,
                max_iter_em = 10,
                tol_em = 1,
                controls_nleqslv = list(),
                controls_kliep = control_kliep(),
                true_matches = NULL) {

  if (!is.null(methods)) {
    stopifnot("`methods` should be a list." =
                is.list(methods))
  }

  if (!is.null(start_params)) {
    stopifnot("`start_params` should be a list." =
                is.list(start_params))
  }

  if (!is.null(true_matches)) {

    if (!(is.data.frame(true_matches) || is.data.table(true_matches))) {
      warning("`true_matches` should be a data.frame or a data.table. Setting `true_matches` to `NULL`.")
      true_matches <- NULL
    } else if (!(length(colnames(true_matches)) == 2 && all(colnames(true_matches) == c("a", "b")))) {
      warning("`true_matches` should consist of two columns: a, b. Setting `true_matches` to `NULL`.")
    }

  }

  if (is.null(set_construction)) {
    set_construction <- "size"
  }

  data.table::setDT(A)
  data.table::setDT(B)
  data.table::set(A, j = "a", value = seq_len(nrow(A)))
  data.table::set(B, j = "b", value = seq_len(nrow(B)))
  M <- merge(A, B, by = variables)

  stopifnot("There are no records with perfect agreement on the key variables.
            Please provide relevant datasets." =
              NROW(M) > 0)

  unique_values <- lapply(variables, function(col) {
    length(unique(c(A[[col]], B[[col]])))
  })
  names(unique_values) <- variables

  for (var in variables) {
    if (unique_values[[var]] == 1) {
      data.table::set(A, j = substitute(var), value = NULL)
      data.table::set(B, j = substitute(var), value = NULL)
      variables <- variables[variables != var]
      comparators[[var]] <- NULL
      methods[[var]] <- NULL
      warning(paste("The variable", var, "has only one unique value and has been removed."))
    }
  }

  vectors <- comparison_vectors(A = A,
                                B = B,
                                variables = variables,
                                comparators = comparators)
  Omega <- vectors$Omega
  comparators <- vectors$comparators

  missing_variables <- variables[!(variables %in% names(methods))]
  methods[missing_variables] <- "binary"
  methods <- methods[variables]

  b_vars <- NULL
  cpar_vars <- NULL
  cnonpar_vars <- NULL
  hm_vars <- NULL

  if (any(methods == "binary")) {
    b_vars <- paste0("gamma_", names(which(methods == "binary")))
  }

  if (any(methods == "continuous_parametric")) {
    cpar_vars <- paste0("gamma_", names(which(methods == "continuous_parametric")))
  }

  if (any(methods == "continuous_nonparametric")) {
    cnonpar_vars <- paste0("gamma_", names(which(methods == "continuous_nonparametric")))
  }

  if (any(methods == "hit_miss")) {
    hm_vars <- paste0("gamma_", names(which(methods == "hit_miss")))
  }

  ratio_kliep <- NULL
  B_values <- NULL

  M <- merge(M, Omega, by = c("a", "b"), all = FALSE)
  M <- M[, colnames(Omega), with = FALSE]
  U <- data.table::fsetdiff(Omega, M)
  n <- NROW(Omega)
  n_M <- NROW(M)
  data.table::set(Omega, j = "ratio", value = 1)

  if (is.null(start_params)) {

    start_params <- list()

    if (length(b_vars) > 0) {
      start_params[["binary"]] <- data.table(
        variable = b_vars,
        theta = runif(length(b_vars), min = 0.9)
      )

    }

    if (length(cpar_vars) > 0) {
      start_params[["continuous_parametric"]] <- data.table(
        variable = cpar_vars,
        p_0_M = runif(length(cpar_vars), min = 0.8, max = 0.9),
        # p_0_M = n_M / rep(min(NROW(A), NROW(B)), length(cpar_vars)),
        alpha_M = runif(length(cpar_vars), min = 0.1, max = 1),
        beta_M = runif(length(cpar_vars), min = 10, max = 20)
      )
    }

    if (length(hm_vars) > 0) {
      start_params[["hit_miss"]] <- data.table(
        variable = hm_vars,
        theta = runif(length(hm_vars), min = 0.9)
      )

    }

  }

  if (length(b_vars) > 0) {

    b_params <- start_params$binary
    Omega_b <- Omega[, b_vars, with = FALSE]
    M_b <- M[, b_vars, with = FALSE]
    U_b <- U[, b_vars, with = FALSE]

    eta_b <- binary_formula(Omega_b)
    b_params$eta <- eta_b

    b_numerator_list <- lapply(b_vars,
                                    function(col) {
                                      stats::dbinom(x = Omega_b[[col]],
                                                    size = 1,
                                                    prob = as.numeric(b_params[b_params[["variable"]] == col, "theta"]))
                                    })
    b_numerator <- Reduce(`*`, b_numerator_list)
    b_denominator_list <- lapply(b_vars,
                                      function(col) {
                                        stats::dbinom(x = Omega_b[[col]],
                                                      size = 1,
                                                      prob = as.numeric(b_params[b_params[["variable"]] == col, "eta"]))
                                      })
    b_denominator <- Reduce(`*`, b_denominator_list)
    data.table::set(Omega, j = "ratio", value = Omega[["ratio"]] * b_numerator / b_denominator)
    data.table::set(Omega, j = "b_denominator", value = b_denominator)

  }

  if (length(cpar_vars) > 0) {

    cpar_params <- start_params$continuous_parametric
    Omega_cpar <- Omega[, cpar_vars, with = FALSE]
    M_cpar <- M[, cpar_vars, with = FALSE]
    U_cpar <- U[, cpar_vars, with = FALSE]

    p_0_U <- p_0_formula(Omega_cpar)
    gamma_plus_U <- gamma_plus_formula(Omega_cpar)
    modified_nleqslv <- purrr::partial(nleqslv::nleqslv, control = controls_nleqslv)
    alpha_U <- alpha_formula(Omega_cpar, modified_nleqslv)
    beta_U <- alpha_U / gamma_plus_U
    cpar_params$p_0_U <- p_0_U
    cpar_params$alpha_U <- alpha_U
    cpar_params$beta_U <- beta_U

    cpar_numerator_list <- lapply(cpar_vars,
                                                   function(col) {
                                                     hurdle_gamma_density(x = Omega_cpar[[col]],
                                                                          p_0 = as.numeric(cpar_params[cpar_params[["variable"]] == col, "p_0_M"]),
                                                                          alpha = as.numeric(cpar_params[cpar_params[["variable"]] == col, "alpha_M"]),
                                                                          beta = as.numeric(cpar_params[cpar_params[["variable"]] == col, "beta_M"]))
                                                   })
    cpar_numerator <- Reduce(`*`, cpar_numerator_list)
    cpar_denominator_list <- lapply(cpar_vars,
                                                     function(col) {
                                                       hurdle_gamma_density(x = Omega_cpar[[col]],
                                                                            p_0 = as.numeric(cpar_params[cpar_params[["variable"]] == col, "p_0_U"]),
                                                                            alpha = as.numeric(cpar_params[cpar_params[["variable"]] == col, "alpha_U"]),
                                                                            beta = as.numeric(cpar_params[cpar_params[["variable"]] == col, "beta_U"]))
                                                     })
    cpar_denominator <- Reduce(`*`, cpar_denominator_list)
    data.table::set(Omega, j = "ratio", value = Omega[["ratio"]] * cpar_numerator / cpar_denominator)
    data.table::set(Omega, j = "cpar_denominator", value = cpar_denominator)

  }

  if (length(cnonpar_vars) > 0) {

    Omega_cnonpar <- Omega[, cnonpar_vars, with = FALSE]
    M_cnonpar <- M[, cnonpar_vars, with = FALSE]
    U_cnonpar <- U[, cnonpar_vars, with = FALSE]

    Omega_indexes <- paste0(Omega[["a"]], "_", Omega[["b"]])
    M_indexes <- paste0(M[["a"]], "_", M[["b"]])

    if (nonpar_hurdle) {

      p_0_M_cnonpar <- runif(length(cnonpar_vars), min = 0.5)
      p_0_U_cnonpar <- p_0_formula(Omega_cnonpar)
      names(p_0_M_cnonpar) <- cnonpar_vars

      ratio_temp <- lapply(cnonpar_vars, function(x) {
        r <- numeric(n)
        r[which(Omega_indexes %in% M_indexes)] <- stats::runif(length(which(Omega_indexes %in% M_indexes)),
                                                               min = 5, max = 10)
        r[setdiff(1:n, which(Omega_indexes %in% M_indexes))] <- stats::runif(n - length(which(Omega_indexes %in% M_indexes)),
                                                                             min = 0.1, max = 1)
        r
      })
      names(ratio_temp) <- cnonpar_vars
      cnonpar_ratio_list <- lapply(cnonpar_vars, function(x) {
        gamma_vec <- Omega_cnonpar[[x]]
        ifelse(gamma_vec == 0, p_0_M_cnonpar[x] / p_0_U_cnonpar[x], 1) *
          ifelse(gamma_vec > 0, (1 - p_0_M_cnonpar[x]) / (1 - p_0_U_cnonpar[x]) * ratio_temp[[x]], 1)
      })
      ratio_kliep <- Reduce(`*`, cnonpar_ratio_list)
      data.table::set(Omega, j = "ratio", value = Omega[["ratio"]] * ratio_kliep)

    } else {

      ratio_temp <- as.numeric(Omega$ratio)
      ratio_temp[which(Omega_indexes %in% M_indexes)] <- (Omega$ratio)[which(Omega_indexes %in% M_indexes)] * stats::runif(length(which(Omega_indexes %in% M_indexes)),
                                                                                                                           min = 5, max = 10)
      ratio_temp[setdiff(1:n, which(Omega_indexes %in% M_indexes))] <- (Omega$ratio)[setdiff(1:n, which(Omega_indexes %in% M_indexes))] * stats::runif(n - length(which(Omega_indexes %in% M_indexes)),
                                                                                                                                                       min = 0.1, max = 5)
      data.table::set(Omega, j = "ratio", value = ratio_temp)

    }

  }

  if (length(hm_vars) > 0) {

    hm_params <- start_params$hit_miss
    Omega_hm <- Omega[, hm_vars, with = FALSE]
    M_hm <- M[, hm_vars, with = FALSE]
    U_hm <- U[, hm_vars, with = FALSE]

    eta_hm <- binary_formula(Omega_hm)
    hm_params$eta <- eta_hm

    hm_numerator_list <- lapply(hm_vars,
                               function(col) {
                                 stats::dbinom(x = Omega_hm[[col]],
                                               size = 1,
                                               prob = as.numeric(hm_params[hm_params[["variable"]] == col, "theta"]))
                               })
    hm_numerator <- Reduce(`*`, hm_numerator_list)
    hm_denominator_list <- lapply(hm_vars,
                                 function(col) {
                                   stats::dbinom(x = Omega_hm[[col]],
                                                 size = 1,
                                                 prob = as.numeric(hm_params[hm_params[["variable"]] == col, "eta"]))
                                 })
    hm_denominator <- Reduce(`*`, hm_denominator_list)
    data.table::set(Omega, j = "ratio", value = Omega[["ratio"]] * hm_numerator / hm_denominator)
    data.table::set(Omega, j = "hm_denominator", value = hm_denominator)

    values_list_m <- lapply(hm_vars, function(x) {
      name <- substr(x, 7, nchar(x))
      values <- unique(c(A[[name]], B[[name]]))
      m_est <- sapply(values, function(y) sum(A[[name]] == y) / NROW(A))
      data.table::data.table(
        "value" = values,
        "m_est" = m_est
      )
    })
    names(values_list_m) <- hm_vars
  }

  iter <- 1

  repeat {

    g_est <- pmin(NROW(M) * Omega$ratio / (NROW(M) * (Omega$ratio - 1) + n), 1)
    n_M_old <- n_M
    n_M <- sum(g_est)

    if (n_M > min(NROW(A), NROW(B))) {
      n_M <- min(NROW(A), NROW(B))
    }

    data.table::set(Omega, j = "g_est", value = g_est)
    Omega <- Omega[order(-get("ratio")), ]

    M <- data.table("a" = numeric(), "b" = numeric())
    for (var in colnames(Omega)) {
      M[[var]] <- numeric()
    }
    M[["ratio"]] <- numeric()

    if (!duplicates_in_A) {

      used_a <- c()
      used_b <- c()

      for (i in 1:NROW(Omega)) {

        current_a <- Omega$a[i]
        current_b <- Omega$b[i]
        if (!(current_a %in% used_a) && !(current_b %in% used_b)) {
          # M <- rbind(M, Omega[i, c("a", "b", "ratio")])
          M <- rbind(M, Omega[i, ])
          used_a <- c(used_a, current_a)
          used_b <- c(used_b, current_b)
        }
        if (NROW(M) >= n_M) {
          break
        }

      }

    } else {

      used_a <- c()

      for (i in 1:NROW(Omega)) {

        current_a <- Omega$a[i]
        if (!(current_a %in% used_a)) {
          M <- rbind(M, Omega[i, ])
          used_a <- c(used_a, current_a)
        }
        if (NROW(M) >= n_M) {
          break
        }

      }

    }

    M <- head(M, round(n_M))
    U <- data.table::fsetdiff(Omega, M)

    if (NROW(M) == 0) {
      break
    }

    if (iter >= 2) {

      old_params <- c()
      params <- c()
      if (length(b_vars) > 0) {
        old_params <- c(old_params, theta_b_old)
        params <- c(params, theta_b)
      }
      if (length(cpar_vars) > 0) {
        old_params <- c(old_params, p_0_M_old, alpha_M_old, beta_M_old)
        params <- c(params, p_0_M, alpha_M, beta_M)
      }
      if (length(hm_vars) > 0) {
        old_params <- c(old_params, theta_hm_old)
        params <- c(params, theta_hm)
      }

      if (length(cnonpar_vars) == 0) {

        if ((abs(n_M_old - n_M) < delta) || norm(old_params - params, type = "2") < eps) {
          break
        }

      } else {

        if ((abs(n_M_old - n_M) < delta)) {
          break
        }

      }

    }

    data.table::set(Omega, j = "ratio", value = 1)

    if (length(b_vars) > 0) {

      Omega_b <- Omega[, b_vars, with = FALSE]
      M_b <- M[, b_vars, with = FALSE]
      U_b <- U[, b_vars, with = FALSE]

      theta_b_old <- b_params$theta
      theta_b <- binary_formula(M_b)
      b_params$theta <- theta_b

      b_numerator_list <- lapply(b_vars,
                                      function(col) {
                                        stats::dbinom(x = Omega_b[[col]],
                                                      size = 1,
                                                      prob = as.numeric(b_params[b_params[["variable"]] == col, "theta"]))
                                      })
      b_numerator <- Reduce(`*`, b_numerator_list)
      data.table::set(Omega, j = "ratio", value = Omega[["ratio"]] * b_numerator / Omega[["b_denominator"]])
    }

    if (length(cpar_vars) > 0) {
      Omega_cpar <- Omega[, cpar_vars, with = FALSE]
      M_cpar <- M[, cpar_vars, with = FALSE]
      U_cpar <- U[, cpar_vars, with = FALSE]
      p_0_M_old <- cpar_params$p_0_M
      alpha_M_old <- cpar_params$alpha_M
      beta_M_old <- cpar_params$beta_M
      p_0_M <- p_0_formula(M_cpar)
      gamma_plus_M <- gamma_plus_formula(M_cpar)
      # alpha_M <- alpha_formula_iterative(M_cpar, modified_nleqslv, beta_M_old)
      alpha_M <- alpha_formula(M_cpar, modified_nleqslv)
      beta_M <- alpha_M / gamma_plus_M
      beta_M[is.nan(beta_M)] <- beta_M_old[is.nan(beta_M)]
      cpar_params$p_0_M <- p_0_M
      cpar_params$alpha_M <- alpha_M
      cpar_params$beta_M <- beta_M
      cpar_numerator_list <- lapply(cpar_vars,
                                                     function(col) {
                                                       hurdle_gamma_density(x = Omega_cpar[[col]],
                                                                            p_0 = as.numeric(cpar_params[cpar_params[["variable"]] == col, "p_0_M"]),
                                                                            alpha = as.numeric(cpar_params[cpar_params[["variable"]] == col, "alpha_M"]),
                                                                            beta = as.numeric(cpar_params[cpar_params[["variable"]] == col, "beta_M"]))
                                                     })

      cpar_numerator <- Reduce(`*`, cpar_numerator_list)
      data.table::set(Omega, j = "ratio", value = Omega[["ratio"]] * cpar_numerator / Omega[["cpar_denominator"]])

    }

    if (length(cnonpar_vars) > 0) {

      Omega_cnonpar <- Omega[, cnonpar_vars, with = FALSE]
      M_cnonpar <- M[, cnonpar_vars, with = FALSE]
      U_cnonpar <- U[, cnonpar_vars, with = FALSE]

      if (nonpar_hurdle) {

        p_0_M_cnonpar <- p_0_formula(M_cnonpar)

        ratio_kliep_list <- lapply(cnonpar_vars, function (x) {
          gamma_M <- M_cnonpar[[x]]
          gamma_M <- data.table::data.table(gamma_M[gamma_M > 0])
          names(gamma_M) <- x
          gamma_U <- U_cnonpar[[x]]
          gamma_U <- data.table::data.table(gamma_U[gamma_U > 0])
          names(gamma_U) <- x
          if (length(gamma_M) > 0) {
            ratio_plus <- do.call(
              densityratio::kliep,
              c(list(
                df_numerator = gamma_M,
                df_denominator = gamma_U
              ),
              controls_kliep)
            )
            gamma_vec <- Omega_cnonpar[[x]]
            gamma_df <- Omega_cnonpar[, x, with = FALSE]
            kliep_pred <- as.vector(stats::predict(ratio_plus, gamma_df))
            ifelse(gamma_vec == 0, p_0_M_cnonpar[x] / p_0_U_cnonpar[x], 1) *
              ifelse(gamma_vec > 0, (1 - p_0_M_cnonpar[x]) * (1 - p_0_U_cnonpar[x]) * kliep_pred, 1)
          } else {
            gamma_vec <- Omega_cnonpar[[x]]
            ifelse(gamma_vec == 0, p_0_M_cnonpar[x] / p_0_U_cnonpar[x], 1)
          }

        })
        ratio_kliep <- Reduce(`*`, ratio_kliep_list)
        data.table::set(Omega, j = "ratio", value = Omega[["ratio"]] * ratio_kliep)

      } else {

        ratio_kliep <- do.call(
          densityratio::kliep,
          c(list(
            df_numerator = M_cnonpar,
            df_denominator = U_cnonpar
          ),
          controls_kliep)
        )

        data.table::set(Omega, j = "ratio", value = Omega[["ratio"]] * stats::predict(ratio_kliep, Omega_cnonpar))

      }

    }

    if (length(hm_vars) > 0) {

      Omega_hm <- Omega[, hm_vars, with = FALSE]
      M_hm <- M[, hm_vars, with = FALSE]
      U_hm <- U[, hm_vars, with = FALSE]

      theta_hm_old <- hm_params$theta
      theta_hm <- binary_formula(M_hm)
      hm_params$theta <- theta_hm

      p_est <- n_M / max(NROW(A), NROW(B))

      u_est_list <- lapply(hm_vars, function(x) {
        u_est <- runif(NROW(values_list_m[[x]]), 0, 1)
        u_est <- u_est / sum(u_est)
      })
      names(u_est_list) <- hm_vars

      m_bk_prod <- sapply(1:NROW(B), function(x) {
        m_bk_list <- lapply(hm_vars, function(y) {
          var_name <- substr(y, 7, nchar(y))
          d_value <- (B[[var_name]])[x]
          m_kd <- as.vector((values_list_m[[y]])[(values_list_m[[y]])[["value"]] == d_value, ][["m_est"]])
        })
        names(m_bk_list) <- hm_vars
        Reduce(`*`, m_bk_list)
      })

      iter_em <- 1
      while (iter_em <= max_iter_em + 1) {

        if (iter_em >= 2) {
          delta_b_old <- em_params[["delta_b"]]
        }

        delta_m_u <- lapply(1:NROW(B), function(x) {

          u_bk_list <- lapply(hm_vars, function(y) {
            var_name <- substr(y, 7, nchar(y))
            values <- cbind(values_list_m[[y]], u_est_list[[y]])
            colnames(values) <- c("value", "m_est", "u_est")
            d_value <- (B[[var_name]])[x]
            u_kd <- as.vector(values[values[["value"]] == d_value, ][["u_est"]])
          })
          names(u_bk_list) <- hm_vars
          u_bk_prod <- Reduce(`*`, u_bk_list)

          data.table::data.table(delta_b = p_est * m_bk_prod[x] / (p_est * m_bk_prod[x] + (1 - p_est) * u_bk_prod),
                                 m_bk_prod = m_bk_prod[x], u_bk_prod = u_bk_prod)

        })

        em_params <- data.table::rbindlist(delta_m_u)

        if (iter_em >= 3) {
          log_lik_old <- log_lik
        }

        if (iter_em >= 2) {
          log_lik <- sum(
            ifelse(delta_b_old == 0,
                   0,
                   delta_b_old * log(p_est * em_params[["m_bk_prod"]]))
          ) +
            sum(ifelse(em_params[["u_bk_prod"]] == 0,
                       0,
                       (1 - delta_b_old) * log((1 - p_est) * em_params[["u_bk_prod"]]))
                )
        }

        if (iter_em >= 3) {
          if (abs(log_lik - log_lik_old) <= tol_em) break
        }

        u_est_list <- lapply(hm_vars, function(x) {

          var_name <- substr(x, 7, nchar(x))
          z_bk <- B[[var_name]]

          sapply((values_list_m[[x]])[["value"]], function(y) {
            sum((1 - em_params[["delta_b"]]) * (z_bk == y)) / sum(1 - em_params[["delta_b"]])
          })

        })
        names(u_est_list) <- hm_vars

        iter_em <- iter_em + 1

      }

      hm_numerator_list <- lapply(hm_vars,
                                 function(col) {
                                   stats::dbinom(x = Omega_hm[[col]],
                                                 size = 1,
                                                 prob = as.numeric(hm_params[hm_params[["variable"]] == col, "theta"]))
                                 })
      hm_numerator <- Reduce(`*`, hm_numerator_list)
      eta_hm <- sapply(hm_vars, function(col) {
        ((1 - p_est) * sum(u_est_list[[col]] * (values_list_m[[col]])[["m_est"]]) +
           p_est * (1 - 1 / NROW(A)) * sum(((values_list_m[[col]])[["m_est"]]) ^ 2)) /
          (1 - p_est / NROW(A))
      })
      hm_params$eta <- eta_hm

      hm_denominator_list <- lapply(hm_vars,
                                    function(col) {
                                      stats::dbinom(x = Omega_hm[[col]],
                                                    size = 1,
                                                    prob = as.numeric(hm_params[hm_params[["variable"]] == col, "eta"]))
                                    })
      hm_denominator <- Reduce(`*`, hm_denominator_list)
      data.table::set(Omega, j = "hm_denominator", value = hm_denominator)

      data.table::set(Omega, j = "ratio", value = Omega[["ratio"]] * hm_numerator / Omega[["hm_denominator"]])
    }

    iter <- iter + 1

  }

  n_M_est <- n_M
  iter_bisection <- NULL

  if (set_construction == "flr") {

    min_treshold <- min(Omega$ratio)
    max_treshold <- max(Omega$ratio)
    treshold <- (min_treshold + max_treshold) / 2

    iter_bisection <- 0

    while (iter_bisection < max_iter_bisection) {

      M <- Omega[get("ratio") >= treshold, ]
      flr_est <- 1 / NROW(M) * sum(1 - M[["g_est"]])

      if (abs(flr_est - target_flr) <= tol) {

        iter_bisection <- iter_bisection + 1
        break

      } else if (flr_est < target_flr) {

        max_treshold <- treshold
        treshold <- (min_treshold + max_treshold) / 2

      } else {

        min_treshold <- treshold
        treshold <- (min_treshold + max_treshold) / 2

      }

      iter_bisection <- iter_bisection + 1

    }

  } else if (set_construction == "size") {

    flr_est <- 1 / NROW(M) * sum(1 - M[["g_est"]])

  }

  mmr_est <- 1 - sum(M[["g_est"]] / n_M_est)
  if (mmr_est < 0) {
    mmr_est <- NULL
  }

  if (!is.null(true_matches)) {

    data.table::setDT(true_matches)

    eval <- evaluation(M, true_matches, n)
    eval_metrics <- unlist(get_metrics(
      TP = eval$TP,
      FP = eval$FP,
      FN = eval$FN,
      TN = eval$TN
    ))
    confusion <- get_confusion(
      TP = eval$TP,
      FP = eval$FP,
      FN = eval$FN,
      TN = eval$TN
    )

  }

  structure(
    list(
      M_est = M[, c("a", "b", "ratio")],
      n_M_est = n_M_est,
      flr_est = flr_est,
      mmr_est = if (is.null(mmr_est)) NULL else mmr_est,
      iter_bisection = iter_bisection,
      b_vars = if (length(b_vars) == 0) NULL else b_vars,
      cpar_vars = if (length(cpar_vars) == 0) NULL else cpar_vars,
      cnonpar_vars = if (length(cnonpar_vars) == 0) NULL else cnonpar_vars,
      hm_vars = if (length(hm_vars) == 0) NULL else hm_vars,
      b_params = if (length(b_vars) == 0) NULL else b_params,
      cpar_params = if (length(cpar_vars) == 0) NULL else cpar_params,
      hm_params = if (length(hm_vars) == 0) NULL else hm_params,
      ratio_kliep = if (is.null(ratio_kliep)) NULL else ratio_kliep,
      variables = variables,
      set_construction = set_construction,
      eval_metrics = if (is.null(true_matches)) NULL else eval_metrics,
      confusion = if (is.null(true_matches)) NULL else confusion
    ),
    class = "mec_rec_lin"
  )

}
