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
#' @param set_construction A method for constructing the predicted set of matches (`"size"`, `"flr"` or `"mmr"`).
#' @param target_rate A target false link rate (FLR) or missing match rate
#' (MMR) (used only if `set_construction == "flr"` or `set_construction == "mmr"`).
#' @param max_iter_bisection A maximum number of iterations for the bisection procedure
#' (used only if `set_construction == "flr"` or `set_construction == "mmr"`).
#' @param tol Error tolerance in the bisection procedure
#' (used only if `set_construction == "flr"` or `set_construction == "mmr"`).
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
#' based on a target False Link Rate (FLR)
#' or missing match rate (MMR). To use the second option, set `set_construction = "flr"`
#' or `set_construction = "mmr"` and
#' specify a target error rate using the `target_rate` argument.
#'
#' The assumption that \eqn{A} and \eqn{B} contain no duplicate records
#' might be relaxed by allowing \eqn{A} to have duplicates. To do so,
#' set `duplicates_in_A = TRUE`.
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
                nonpar_hurdle = TRUE,
                set_construction = NULL,
                target_rate = 0.03,
                max_iter_bisection = 100,
                tol = 0.005,
                delta = 0.5,
                eps = 0.05,
                max_iter_em = 10,
                tol_em = 1,
                controls_nleqslv = list(),
                controls_kliep = control_kliep(),
                true_matches = NULL) {

  if (!is.null(start_params)) {
    stopifnot("`start_params` should be a list." =
                is.list(start_params))
  }

  if (is.null(set_construction)) {
    set_construction <- "size"
  }
  validate_choice(set_construction, c("size", "flr", "mmr"), "set_construction")
  methods <- validate_methods(
    methods = methods,
    variables = variables,
    allowed_methods = c("binary", "continuous_parametric", "continuous_nonparametric", "hit_miss")
  )

  data.table::setDT(A)
  data.table::setDT(B)
  true_matches <- sanitize_true_matches(
    true_matches = true_matches,
    n_A = nrow(A),
    n_B = nrow(B),
    arg_name = "true_matches"
  )
  stopifnot("There are no records with perfect agreement on the key variables.
            Please provide relevant datasets." =
              has_perfect_agreement(A, B, variables))
  data.table::set(A, j = "a", value = seq_len(nrow(A)))
  data.table::set(B, j = "b", value = seq_len(nrow(B)))
  M <- merge(
    A[, c("a", variables), with = FALSE],
    B[, c("b", variables), with = FALSE],
    by = variables
  )

  preprocessing <- drop_constant_key_variables(
    A = A,
    B = B,
    variables = variables,
    comparators = comparators,
    methods = methods
  )
  A <- preprocessing$A
  B <- preprocessing$B
  variables <- preprocessing$variables
  comparators <- preprocessing$comparators
  methods <- preprocessing$methods
  constant_vars <- preprocessing$constant_vars

  for (var in constant_vars) {
    warning(paste("The variable", var, "has only one unique value and has been removed."))
  }

  vectors <- comparison_vectors(A = A,
                                B = B,
                                variables = variables,
                                comparators = comparators)
  Omega <- vectors$Omega
  comparators <- vectors$comparators
  method_variables <- extract_method_variables(methods, include_hit_miss = TRUE)
  b_vars <- method_variables$b_vars
  cpar_vars <- method_variables$cpar_vars
  cnonpar_vars <- method_variables$cnonpar_vars
  hm_vars <- method_variables$hm_vars

  ratio_kliep <- NULL
  kliep_warning_emitted <- FALSE

  warn_kliep_once <- function(variables, message) {
    if (!kliep_warning_emitted) {
      warn_kliep_issue("mec()", variables, message)
      kliep_warning_emitted <<- TRUE
    }
  }

  exact_match_idx <- data.table(
    omega_idx = seq_len(NROW(Omega)),
    a = Omega[["a"]],
    b = Omega[["b"]]
  )[M[, c("a", "b"), with = FALSE], on = c("a", "b"), nomatch = 0L][["omega_idx"]]
  n <- NROW(Omega)
  exact_match_mask <- rep(FALSE, n)
  exact_match_mask[exact_match_idx] <- TRUE
  M_idx <- exact_match_idx
  M <- Omega[M_idx]
  n_M <- length(M_idx)
  all_idx <- seq_len(n)
  selection_mask <- rep(FALSE, n)
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

    b_params <- align_parameter_table(start_params$binary, b_vars)
    Omega_b <- Omega[, b_vars, with = FALSE]

    eta_b <- binary_formula(Omega_b)
    b_params$eta <- eta_b

    b_denominator <- bernoulli_product(Omega_b, b_params$eta)
    data.table::set(
      Omega,
      j = "ratio",
      value = Omega[["ratio"]] * bernoulli_ratio(Omega_b, b_params$theta, b_params$eta)
    )
    data.table::set(Omega, j = "b_denominator", value = b_denominator)

  }

  if (length(cpar_vars) > 0) {

    cpar_params <- align_parameter_table(start_params$continuous_parametric, cpar_vars)
    Omega_cpar <- Omega[, cpar_vars, with = FALSE]

    p_0_U <- p_0_formula(Omega_cpar)
    gamma_plus_U <- gamma_plus_formula(Omega_cpar)
    modified_nleqslv <- purrr::partial(nleqslv::nleqslv, control = controls_nleqslv)
    alpha_U <- alpha_formula(Omega_cpar, modified_nleqslv)
    beta_U <- alpha_U / gamma_plus_U
    cpar_params$p_0_U <- p_0_U
    cpar_params$alpha_U <- alpha_U
    cpar_params$beta_U <- beta_U

    cpar_denominator <- hurdle_gamma_product(
      Omega_cpar,
      cpar_params$p_0_U,
      cpar_params$alpha_U,
      cpar_params$beta_U
    )
    data.table::set(
      Omega,
      j = "ratio",
      value = Omega[["ratio"]] * hurdle_gamma_ratio(
        Omega_cpar,
        cpar_params$p_0_M,
        cpar_params$alpha_M,
        cpar_params$beta_M,
        cpar_params$p_0_U,
        cpar_params$alpha_U,
        cpar_params$beta_U
      )
    )
    data.table::set(Omega, j = "cpar_denominator", value = cpar_denominator)

  }

  if (length(cnonpar_vars) > 0) {

    Omega_cnonpar <- Omega[, cnonpar_vars, with = FALSE]
    n_exact_matches <- sum(exact_match_mask)

    if (nonpar_hurdle) {

      p_0_M_cnonpar <- runif(length(cnonpar_vars), min = 0.5)
      p_0_U_cnonpar <- p_0_formula(Omega_cnonpar)
      names(p_0_M_cnonpar) <- cnonpar_vars

      ratio_temp <- lapply(cnonpar_vars, function(x) {
        r <- numeric(n)
        r[exact_match_mask] <- stats::runif(n_exact_matches, min = 5, max = 10)
        r[!exact_match_mask] <- stats::runif(n - n_exact_matches, min = 0.1, max = 1)
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
      ratio_temp[exact_match_mask] <- ratio_temp[exact_match_mask] * stats::runif(n_exact_matches,
                                                                                  min = 5, max = 10)
      ratio_temp[!exact_match_mask] <- ratio_temp[!exact_match_mask] * stats::runif(n - n_exact_matches,
                                                                                    min = 0.1, max = 5)
      data.table::set(Omega, j = "ratio", value = ratio_temp)

    }

  }

  if (length(hm_vars) > 0) {

    hm_params <- align_parameter_table(start_params$hit_miss, hm_vars)
    Omega_hm <- Omega[, hm_vars, with = FALSE]

    eta_hm <- binary_formula(Omega_hm)
    hm_params$eta <- eta_hm

    hm_denominator <- bernoulli_product(Omega_hm, hm_params$eta)
    data.table::set(
      Omega,
      j = "ratio",
      value = Omega[["ratio"]] * bernoulli_ratio(Omega_hm, hm_params$theta, hm_params$eta)
    )
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

    g_est <- pmin(length(M_idx) * Omega$ratio / (length(M_idx) * (Omega$ratio - 1) + n), 1)
    n_M_old <- n_M
    n_M <- sum(g_est)

    if (n_M > min(NROW(A), NROW(B))) {
      n_M <- min(NROW(A), NROW(B))
    }

    data.table::set(Omega, j = "g_est", value = g_est)
    M_idx <- select_mec_indices(
      a = Omega[["a"]],
      b = Omega[["b"]],
      ratio = Omega[["ratio"]],
      n_M = n_M,
      duplicates_in_A = duplicates_in_A
    )
    M <- Omega[M_idx]

    if (length(M_idx) == 0L) {
      break
    }

    selection_mask[M_idx] <- TRUE
    U_idx <- all_idx[!selection_mask]
    selection_mask[M_idx] <- FALSE

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

      M_b <- Omega_b[M_idx]
      theta_b_old <- b_params$theta
      theta_b <- binary_formula(M_b)
      b_params$theta <- theta_b

      b_numerator <- bernoulli_product(Omega_b, b_params$theta)
      data.table::set(Omega, j = "ratio", value = Omega[["ratio"]] * b_numerator / Omega[["b_denominator"]])
    }

    if (length(cpar_vars) > 0) {
      M_cpar <- Omega_cpar[M_idx]
      p_0_M_old <- cpar_params$p_0_M
      alpha_M_old <- cpar_params$alpha_M
      beta_M_old <- cpar_params$beta_M
      p_0_M <- p_0_formula(M_cpar)
      gamma_plus_M <- gamma_plus_formula(M_cpar)
      alpha_M <- alpha_formula(M_cpar, modified_nleqslv)
      beta_M <- alpha_M / gamma_plus_M
      beta_M[is.nan(beta_M)] <- beta_M_old[is.nan(beta_M)]
      cpar_params$p_0_M <- p_0_M
      cpar_params$alpha_M <- alpha_M
      cpar_params$beta_M <- beta_M
      cpar_numerator <- hurdle_gamma_product(
        Omega_cpar,
        cpar_params$p_0_M,
        cpar_params$alpha_M,
        cpar_params$beta_M
      )
      data.table::set(Omega, j = "ratio", value = Omega[["ratio"]] * cpar_numerator / Omega[["cpar_denominator"]])

    }

    if (length(cnonpar_vars) > 0) {

      M_cnonpar <- Omega_cnonpar[M_idx]
      U_cnonpar <- Omega_cnonpar[U_idx]

      if (nonpar_hurdle) {

        # When the current split is too sparse for KLIEP, keep the last ratio update.
        p_0_M_cnonpar <- p_0_formula(M_cnonpar)
        p_0_U_cnonpar <- p_0_formula(U_cnonpar)
        ratio_kliep_old <- ratio_kliep

        tryCatch({
          ratio_kliep_models <- fit_kliep_hurdle_models(
            df_numerator = M_cnonpar,
            df_denominator = U_cnonpar,
            variables = cnonpar_vars,
            controls_kliep = controls_kliep
          )
          missing_models <- missing_kliep_models(ratio_kliep_models)

          if (length(missing_models) < length(cnonpar_vars)) {
            if (length(missing_models) > 0L) {
              warn_kliep_once(missing_models, "using only the hurdle mass term for those variables in the current iteration.")
            }

            ratio_kliep <- kliep_hurdle_ratio(
              df = Omega_cnonpar,
              variables = cnonpar_vars,
              p_0_numerator = p_0_M_cnonpar,
              p_0_denominator = p_0_U_cnonpar,
              ratio_kliep_list = ratio_kliep_models
            )
          } else {
            ratio_kliep <- ratio_kliep_old
            warn_kliep_once(cnonpar_vars, "could not be fitted in the current iteration; using the previous ratio estimate.")
          }
          data.table::set(Omega, j = "ratio", value = Omega[["ratio"]] * ratio_kliep)
        },
          error = function(e) {
            warn_kliep_once(cnonpar_vars, paste("failed with error:", e$message))
          })

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

      M_hm <- Omega_hm[M_idx]
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

      hm_numerator <- bernoulli_product(Omega_hm, hm_params$theta)
      eta_hm <- sapply(hm_vars, function(col) {
        ((1 - p_est) * sum(u_est_list[[col]] * (values_list_m[[col]])[["m_est"]]) +
           p_est * (1 - 1 / NROW(A)) * sum(((values_list_m[[col]])[["m_est"]]) ^ 2)) /
          (1 - p_est / NROW(A))
      })
      hm_params$eta <- eta_hm

      hm_denominator <- bernoulli_product(Omega_hm, hm_params$eta)
      data.table::set(Omega, j = "hm_denominator", value = hm_denominator)

      data.table::set(Omega, j = "ratio", value = Omega[["ratio"]] * hm_numerator / Omega[["hm_denominator"]])
    }

    iter <- iter + 1

  }

  n_M_est <- n_M
  selection_summary <- summarize_mec_selection(
    a = Omega[["a"]],
    b = Omega[["b"]],
    ratio = Omega[["ratio"]],
    g_est = Omega[["g_est"]],
    n_M_est = n_M_est,
    duplicates_in_A = duplicates_in_A,
    set_construction = set_construction,
    target_rate = target_rate,
    tol = tol,
    max_iter = max_iter_bisection
  )
  M <- Omega[selection_summary$selected_idx]
  flr_est <- selection_summary$flr_est
  mmr_est <- selection_summary$mmr_est
  iter_bisection <- selection_summary$iter

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
      n = n,
      b_vars = if (length(b_vars) == 0) NULL else b_vars,
      cpar_vars = if (length(cpar_vars) == 0) NULL else cpar_vars,
      cnonpar_vars = if (length(cnonpar_vars) == 0) NULL else cnonpar_vars,
      hm_vars = if (length(hm_vars) == 0) NULL else hm_vars,
      b_params = if (length(b_vars) == 0) NULL else b_params,
      cpar_params = if (length(cpar_vars) == 0) NULL else cpar_params,
      hm_params = if (length(hm_vars) == 0) NULL else hm_params,
      ratio_kliep = if (is.null(ratio_kliep)) NULL else ratio_kliep,
      delta = delta,
      eps = eps,
      controls_nleqslv = controls_nleqslv,
      variables = variables,
      set_construction = set_construction,
      eval_metrics = if (is.null(true_matches)) NULL else eval_metrics,
      confusion = if (is.null(true_matches)) NULL else confusion
    ),
    class = "mec_rec_lin"
  )

}

#' @noRd
fit_mec_unsupervised_omega <- function(A,
                                       B,
                                       variables,
                                       comparators,
                                       methods,
                                       Omega,
                                       exact_matches,
                                       start_params = NULL,
                                       nonmatch_params = NULL,
                                       nonpar_hurdle = TRUE,
                                       delta = 0.5,
                                       eps = 0.05,
                                       max_iter_em = 10,
                                       tol_em = 1,
                                       controls_nleqslv = list(),
                                       controls_kliep = control_kliep(),
                                       context = "mec_blocking()") {
  if (!is.null(start_params)) {
    stopifnot("`start_params` should be a list." =
                is.list(start_params))
  }
  if (!is.null(nonmatch_params)) {
    stopifnot("`nonmatch_params` should be a list." =
                is.list(nonmatch_params))
  }

  method_variables <- extract_method_variables(methods, include_hit_miss = TRUE)
  b_vars <- method_variables$b_vars
  cpar_vars <- method_variables$cpar_vars
  cnonpar_vars <- method_variables$cnonpar_vars
  hm_vars <- method_variables$hm_vars

  exact_match_idx <- data.table(
    omega_idx = seq_len(NROW(Omega)),
    a = Omega[["a"]],
    b = Omega[["b"]]
  )[exact_matches[, c("a", "b"), with = FALSE], on = c("a", "b"), nomatch = 0L][["omega_idx"]]

  if (length(exact_match_idx) == 0L) {
    stop("The MEC training space contains no exact agreements on the key variables.")
  }

  n <- NROW(Omega)
  exact_match_mask <- rep(FALSE, n)
  exact_match_mask[exact_match_idx] <- TRUE
  M_idx <- exact_match_idx
  n_M <- length(M_idx)
  all_idx <- seq_len(n)
  selection_mask <- rep(FALSE, n)

  b_params <- NULL
  cpar_params <- NULL
  cnonpar_params <- NULL
  hm_params <- NULL
  ratio_kliep <- NULL
  ratio_kliep_list <- NULL
  kliep_warning_emitted <- FALSE

  warn_kliep_once <- function(variables, message) {
    if (!kliep_warning_emitted) {
      warn_kliep_issue(context, variables, message)
      kliep_warning_emitted <<- TRUE
    }
  }

  data.table::set(Omega, j = "ratio", value = 1)

  if (is.null(start_params)) {
    start_params <- list()

    if (length(b_vars) > 0L) {
      start_params[["binary"]] <- data.table(
        variable = b_vars,
        theta = runif(length(b_vars), min = 0.9)
      )
    }

    if (length(cpar_vars) > 0L) {
      start_params[["continuous_parametric"]] <- data.table(
        variable = cpar_vars,
        p_0_M = runif(length(cpar_vars), min = 0.8, max = 0.9),
        alpha_M = runif(length(cpar_vars), min = 0.1, max = 1),
        beta_M = runif(length(cpar_vars), min = 10, max = 20)
      )
    }

    if (length(hm_vars) > 0L) {
      start_params[["hit_miss"]] <- data.table(
        variable = hm_vars,
        theta = runif(length(hm_vars), min = 0.9)
      )
    }
  }

  if (length(b_vars) > 0L) {
    b_params <- align_parameter_table(start_params$binary, b_vars)
    Omega_b <- Omega[, b_vars, with = FALSE]
    if (!is.null(nonmatch_params) && !is.null(nonmatch_params$binary)) {
      eta_b <- align_parameter_table(nonmatch_params$binary, b_vars)$eta
    } else {
      eta_b <- binary_formula(Omega_b)
    }
    b_params$eta <- eta_b
    b_denominator <- bernoulli_product(Omega_b, b_params$eta)
    data.table::set(
      Omega,
      j = "ratio",
      value = Omega[["ratio"]] * bernoulli_ratio(Omega_b, b_params$theta, b_params$eta)
    )
    data.table::set(Omega, j = "b_denominator", value = b_denominator)
  }

  if (length(cpar_vars) > 0L) {
    cpar_params <- align_parameter_table(start_params$continuous_parametric, cpar_vars)
    Omega_cpar <- Omega[, cpar_vars, with = FALSE]
    modified_nleqslv <- purrr::partial(nleqslv::nleqslv, control = controls_nleqslv)
    if (!is.null(nonmatch_params) && !is.null(nonmatch_params$continuous_parametric)) {
      cpar_nonmatch_params <- align_parameter_table(nonmatch_params$continuous_parametric, cpar_vars)
      p_0_U <- cpar_nonmatch_params$p_0_U
      alpha_U <- cpar_nonmatch_params$alpha_U
      beta_U <- cpar_nonmatch_params$beta_U
    } else {
      p_0_U <- p_0_formula(Omega_cpar)
      gamma_plus_U <- gamma_plus_formula(Omega_cpar)
      alpha_U <- alpha_formula(Omega_cpar, modified_nleqslv)
      beta_U <- alpha_U / gamma_plus_U
    }
    cpar_params$p_0_U <- p_0_U
    cpar_params$alpha_U <- alpha_U
    cpar_params$beta_U <- beta_U
    cpar_denominator <- hurdle_gamma_product(
      Omega_cpar,
      cpar_params$p_0_U,
      cpar_params$alpha_U,
      cpar_params$beta_U
    )
    data.table::set(
      Omega,
      j = "ratio",
      value = Omega[["ratio"]] * hurdle_gamma_ratio(
        Omega_cpar,
        cpar_params$p_0_M,
        cpar_params$alpha_M,
        cpar_params$beta_M,
        cpar_params$p_0_U,
        cpar_params$alpha_U,
        cpar_params$beta_U
      )
    )
    data.table::set(Omega, j = "cpar_denominator", value = cpar_denominator)
  }

  if (length(cnonpar_vars) > 0L) {
    Omega_cnonpar <- Omega[, cnonpar_vars, with = FALSE]
    n_exact_matches <- sum(exact_match_mask)

    if (nonpar_hurdle) {
      p_0_M_cnonpar <- runif(length(cnonpar_vars), min = 0.5)
      p_0_U_cnonpar <- p_0_formula(Omega_cnonpar)
      names(p_0_M_cnonpar) <- cnonpar_vars
      cnonpar_params <- data.table(
        variable = cnonpar_vars,
        p_0_M_cnonpar = p_0_M_cnonpar,
        p_0_U_cnonpar = p_0_U_cnonpar
      )

      ratio_temp <- lapply(cnonpar_vars, function(x) {
        r <- numeric(n)
        r[exact_match_mask] <- stats::runif(n_exact_matches, min = 5, max = 10)
        r[!exact_match_mask] <- stats::runif(n - n_exact_matches, min = 0.1, max = 1)
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
      ratio_temp[exact_match_mask] <- ratio_temp[exact_match_mask] * stats::runif(
        n_exact_matches,
        min = 5,
        max = 10
      )
      ratio_temp[!exact_match_mask] <- ratio_temp[!exact_match_mask] * stats::runif(
        n - n_exact_matches,
        min = 0.1,
        max = 5
      )
      data.table::set(Omega, j = "ratio", value = ratio_temp)
    }
  }

  if (length(hm_vars) > 0L) {
    hm_params <- align_parameter_table(start_params$hit_miss, hm_vars)
    Omega_hm <- Omega[, hm_vars, with = FALSE]
    eta_hm <- binary_formula(Omega_hm)
    hm_params$eta <- eta_hm
    hm_denominator <- bernoulli_product(Omega_hm, hm_params$eta)
    data.table::set(
      Omega,
      j = "ratio",
      value = Omega[["ratio"]] * bernoulli_ratio(Omega_hm, hm_params$theta, hm_params$eta)
    )
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
    g_est <- pmin(length(M_idx) * Omega$ratio / (length(M_idx) * (Omega$ratio - 1) + n), 1)
    n_M_old <- n_M
    n_M <- sum(g_est)

    if (n_M > min(NROW(A), NROW(B))) {
      n_M <- min(NROW(A), NROW(B))
    }

    data.table::set(Omega, j = "g_est", value = g_est)
    M_idx <- select_mec_indices(
      a = Omega[["a"]],
      b = Omega[["b"]],
      ratio = Omega[["ratio"]],
      n_M = n_M,
      duplicates_in_A = FALSE
    )

    if (length(M_idx) == 0L) {
      break
    }

    selection_mask[M_idx] <- TRUE
    U_idx <- all_idx[!selection_mask]
    selection_mask[M_idx] <- FALSE

    if (iter >= 2) {
      old_params <- c()
      params <- c()
      if (length(b_vars) > 0L) {
        old_params <- c(old_params, theta_b_old)
        params <- c(params, theta_b)
      }
      if (length(cpar_vars) > 0L) {
        old_params <- c(old_params, p_0_M_old, alpha_M_old, beta_M_old)
        params <- c(params, p_0_M, alpha_M, beta_M)
      }
      if (length(hm_vars) > 0L) {
        old_params <- c(old_params, theta_hm_old)
        params <- c(params, theta_hm)
      }

      if (length(cnonpar_vars) == 0L) {
        if ((abs(n_M_old - n_M) < delta) || norm(old_params - params, type = "2") < eps) {
          break
        }
      } else if ((abs(n_M_old - n_M) < delta)) {
        break
      }
    }

    data.table::set(Omega, j = "ratio", value = 1)

    if (length(b_vars) > 0L) {
      M_b <- Omega_b[M_idx]
      theta_b_old <- b_params$theta
      theta_b <- binary_formula(M_b)
      b_params$theta <- theta_b
      b_numerator <- bernoulli_product(Omega_b, b_params$theta)
      data.table::set(Omega, j = "ratio", value = Omega[["ratio"]] * b_numerator / Omega[["b_denominator"]])
    }

    if (length(cpar_vars) > 0L) {
      M_cpar <- Omega_cpar[M_idx]
      p_0_M_old <- cpar_params$p_0_M
      alpha_M_old <- cpar_params$alpha_M
      beta_M_old <- cpar_params$beta_M
      p_0_M <- p_0_formula(M_cpar)
      gamma_plus_M <- gamma_plus_formula(M_cpar)
      alpha_M <- alpha_formula(M_cpar, modified_nleqslv)
      beta_M <- alpha_M / gamma_plus_M
      beta_M[is.nan(beta_M)] <- beta_M_old[is.nan(beta_M)]
      cpar_params$p_0_M <- p_0_M
      cpar_params$alpha_M <- alpha_M
      cpar_params$beta_M <- beta_M
      cpar_numerator <- hurdle_gamma_product(
        Omega_cpar,
        cpar_params$p_0_M,
        cpar_params$alpha_M,
        cpar_params$beta_M
      )
      data.table::set(Omega, j = "ratio", value = Omega[["ratio"]] * cpar_numerator / Omega[["cpar_denominator"]])
    }

    if (length(cnonpar_vars) > 0L) {
      M_cnonpar <- Omega_cnonpar[M_idx]
      U_cnonpar <- Omega_cnonpar[U_idx]

      if (nonpar_hurdle) {
        p_0_M_cnonpar <- p_0_formula(M_cnonpar)
        p_0_U_cnonpar <- p_0_formula(U_cnonpar)
        cnonpar_params <- data.table(
          variable = cnonpar_vars,
          p_0_M_cnonpar = p_0_M_cnonpar,
          p_0_U_cnonpar = p_0_U_cnonpar
        )
        ratio_kliep_old <- ratio_kliep

        tryCatch({
          ratio_kliep_models <- fit_kliep_hurdle_models(
            df_numerator = M_cnonpar,
            df_denominator = U_cnonpar,
            variables = cnonpar_vars,
            controls_kliep = controls_kliep
          )
          missing_models <- missing_kliep_models(ratio_kliep_models)

          if (length(missing_models) < length(cnonpar_vars)) {
            if (length(missing_models) > 0L) {
              warn_kliep_once(missing_models, "using only the hurdle mass term for those variables in the current iteration.")
            }

            ratio_kliep_list <- ratio_kliep_models
            ratio_kliep <- kliep_hurdle_ratio(
              df = Omega_cnonpar,
              variables = cnonpar_vars,
              p_0_numerator = p_0_M_cnonpar,
              p_0_denominator = p_0_U_cnonpar,
              ratio_kliep_list = ratio_kliep_models
            )
          } else {
            ratio_kliep <- ratio_kliep_old
            warn_kliep_once(cnonpar_vars, "could not be fitted in the current iteration; using the previous ratio estimate.")
          }
          data.table::set(Omega, j = "ratio", value = Omega[["ratio"]] * ratio_kliep)
        },
          error = function(e) {
            warn_kliep_once(cnonpar_vars, paste("failed with error:", e$message))
          })
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

    if (length(hm_vars) > 0L) {
      M_hm <- Omega_hm[M_idx]
      theta_hm_old <- hm_params$theta
      theta_hm <- binary_formula(M_hm)
      hm_params$theta <- theta_hm

      p_est <- n_M / max(NROW(A), NROW(B))
      u_est_list <- lapply(hm_vars, function(x) {
        u_est <- runif(NROW(values_list_m[[x]]), 0, 1)
        u_est <- u_est / sum(u_est)
      })
      names(u_est_list) <- hm_vars

      m_bk_prod <- sapply(seq_len(NROW(B)), function(x) {
        m_bk_list <- lapply(hm_vars, function(y) {
          var_name <- substr(y, 7, nchar(y))
          d_value <- (B[[var_name]])[x]
          as.vector((values_list_m[[y]])[(values_list_m[[y]])[["value"]] == d_value, ][["m_est"]])
        })
        names(m_bk_list) <- hm_vars
        Reduce(`*`, m_bk_list)
      })

      iter_em <- 1
      while (iter_em <= max_iter_em + 1) {
        if (iter_em >= 2) {
          delta_b_old <- em_params[["delta_b"]]
        }

        delta_m_u <- lapply(seq_len(NROW(B)), function(x) {
          u_bk_list <- lapply(hm_vars, function(y) {
            var_name <- substr(y, 7, nchar(y))
            values <- cbind(values_list_m[[y]], u_est_list[[y]])
            colnames(values) <- c("value", "m_est", "u_est")
            d_value <- (B[[var_name]])[x]
            as.vector(values[values[["value"]] == d_value, ][["u_est"]])
          })
          names(u_bk_list) <- hm_vars
          u_bk_prod <- Reduce(`*`, u_bk_list)
          data.table::data.table(
            delta_b = p_est * m_bk_prod[x] / (p_est * m_bk_prod[x] + (1 - p_est) * u_bk_prod),
            m_bk_prod = m_bk_prod[x],
            u_bk_prod = u_bk_prod
          )
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

      hm_numerator <- bernoulli_product(Omega_hm, hm_params$theta)
      eta_hm <- sapply(hm_vars, function(col) {
        ((1 - p_est) * sum(u_est_list[[col]] * (values_list_m[[col]])[["m_est"]]) +
           p_est * (1 - 1 / NROW(A)) * sum(((values_list_m[[col]])[["m_est"]]) ^ 2)) /
          (1 - p_est / NROW(A))
      })
      hm_params$eta <- eta_hm
      hm_denominator <- bernoulli_product(Omega_hm, hm_params$eta)
      data.table::set(Omega, j = "hm_denominator", value = hm_denominator)
      data.table::set(Omega, j = "ratio", value = Omega[["ratio"]] * hm_numerator / Omega[["hm_denominator"]])
    }

    iter <- iter + 1
  }

  n_M_est <- n_M
  selection_summary <- summarize_mec_selection(
    a = Omega[["a"]],
    b = Omega[["b"]],
    ratio = Omega[["ratio"]],
    g_est = Omega[["g_est"]],
    n_M_est = n_M_est,
    duplicates_in_A = FALSE,
    set_construction = "size"
  )
  M <- Omega[selection_summary$selected_idx]

  model <- list(
    b_vars = if (length(b_vars) == 0L) NULL else b_vars,
    cpar_vars = if (length(cpar_vars) == 0L) NULL else cpar_vars,
    cnonpar_vars = if (length(cnonpar_vars) == 0L) NULL else cnonpar_vars,
    hm_vars = if (length(hm_vars) == 0L) NULL else hm_vars,
    b_params = if (is.null(b_params)) NULL else b_params,
    cpar_params = if (is.null(cpar_params)) NULL else cpar_params,
    cnonpar_params = if (is.null(cnonpar_params)) NULL else cnonpar_params,
    hm_params = if (is.null(hm_params)) NULL else hm_params,
    ratio_kliep = if (is.null(ratio_kliep) || is.numeric(ratio_kliep)) NULL else ratio_kliep,
    ratio_kliep_list = if (is.null(ratio_kliep_list)) NULL else ratio_kliep_list,
    variables = variables,
    comparators = comparators,
    methods = methods,
    n_M_est = n_M_est,
    prob_est = n_M_est / n
  )

  list(
    model = model,
    Omega = Omega,
    M_est = M[, c("a", "b", "ratio"), with = FALSE],
    n_M_est = n_M_est,
    flr_est = selection_summary$flr_est,
    mmr_est = selection_summary$mmr_est
  )
}

#' @title Blocked Unsupervised Maximum Entropy Classifier for Record Linkage
#'
#' @author Adam Struzik
#'
#' @description
#' Runs graph-based blocking using \link[blocking:blocking]{blocking::blocking()},
#' fits one pooled unsupervised maximum entropy classifier (MEC) on selected
#' within-block pairs, and applies the fitted density-ratio model blockwise.
#'
#' @param A A duplicate-free `data.frame` or `data.table`.
#' @param B A duplicate-free `data.frame` or `data.table`.
#' @param variables A character vector of key variables used to create MEC comparison vectors.
#' @param comparators A named list of functions for comparing pairs of records.
#' @param methods A named list of methods used for estimation (`"binary"` or
#' `"continuous_parametric"`). Other unsupervised MEC methods are not
#' supported by `mec_blocking()` at this stage.
#' @param blocking_x Optional input passed as `x` to \link[blocking:blocking]{blocking::blocking()}.
#' @param blocking_y Optional input passed as `y` to \link[blocking:blocking]{blocking::blocking()}.
#' @param blocking_variables Variables used to create blocking strings when
#' `blocking_x` and `blocking_y` are not supplied.
#' @param blocking_sep Separator used when concatenating `blocking_variables`.
#' @param controls_blocking A list of additional arguments passed to
#' \link[blocking:blocking]{blocking::blocking()}, except `x` and `y`.
#' @param min_training_pairs Minimum number of within-block training pairs.
#' If `NULL`, `min_training_nonmatches` should also be `NULL`, and all blocks
#' are used for training.
#' @param min_training_nonmatches Minimum lower bound on within-block nonmatches.
#' If `NULL`, `min_training_pairs` should also be `NULL`, and all blocks are
#' used for training.
#' @param block_sampling_seed Optional seed for random training-block sampling.
#' @param nonmatch_sample_size Number of pairs sampled from the full
#' Cartesian product of `A` and `B` to estimate nonmatch distribution
#' parameters. If `NULL`, all pairs are used.
#' @param nonmatch_sampling_seed Optional seed for nonmatch pair sampling.
#' @param prob_ratio Probability/density ratio type (`"1"` or `"2"`). The
#' default `"2"` uses the blockwise fixed-point equation.
#' @param start_params Start parameters for the `"binary"` and
#' `"continuous_parametric"` methods.
#' @param nonpar_hurdle Currently unused in `mec_blocking()`.
#' @param fixed_method A method for solving blockwise fixed-point equations using
#' the \link[FixedPoint]{FixedPoint} function.
#' @param delta A numeric value specifying the tolerance for the change in the
#' estimated number of matches between MEC iterations.
#' @param eps A numeric value specifying the tolerance for the change in model
#' parameters between MEC iterations.
#' @param max_iter_em Currently unused in `mec_blocking()`.
#' @param tol_em Currently unused in `mec_blocking()`.
#' @param controls_nleqslv Controls passed to the \link[nleqslv]{nleqslv} function
#' (only if the `"continuous_parametric"` method has been chosen for at least one
#' variable).
#' @param controls_kliep Currently unused in `mec_blocking()`.
#' @param true_matches A `data.frame` or `data.table` indicating known matches.
#' @param keep_blocking_result Logical indicating whether to store the raw object
#' returned by \link[blocking:blocking]{blocking::blocking()}.
#' @param keep_training_data Logical indicating whether to store pooled training
#' comparison vectors.
#' @param verbose Logical indicating whether to print progress messages.
#'
#' @details
#' The function assumes one-to-one linkage. The blocking stage defines disjoint
#' bipartite blocks. MEC is trained once on the pooled union of selected
#' within-block Cartesian products and is then applied separately in each final
#' block. The ANN distance returned by `blocking::blocking()` is not used.
#'
#' If both `min_training_pairs` and `min_training_nonmatches` are `NULL`, all
#' final blocks are used for pooled training. If both are supplied, blocks are
#' sampled without replacement until both thresholds are reached. Supplying only
#' one threshold is an error.
#'
#' Nonmatch distribution parameters are estimated from a simple random sample
#' from the full Cartesian product of `A` and `B`, not from the selected
#' training blocks.
#'
#' @return
#' Returns an object of class `"mec_blocking"` containing predicted matches,
#' estimated match counts, block summaries, selected training blocks, fitted
#' density-ratio components, and optional evaluation diagnostics.
#'
#' @examples
#' df_1 <- data.frame(
#'   name = c("Emma", "Liam", "Olivia", "Noah", "Ava"),
#'   surname = c("Smith", "Jones", "Brown", "Davis", "Miller"),
#'   city = c("Boston", "Boston", "Austin", "Austin", "Denver")
#' )
#' df_2 <- data.frame(
#'   name = c("Emma", "Liam", "Olivia", "Noah", "Ava"),
#'   surname = c("Smith", "Jones", "Brown", "Davis", "Miller"),
#'   city = c("Boston", "Boston", "Austin", "Austin", "Denver")
#' )
#'
#' blocking_x <- matrix(
#'   c(1, 0, 0, 1, 1, 1, 2, 0, 0, 2),
#'   ncol = 2,
#'   byrow = TRUE
#' )
#' blocking_y <- blocking_x
#'
#' result <- mec_blocking(
#'   A = df_1,
#'   B = df_2,
#'   variables = c("name", "surname", "city"),
#'   blocking_x = blocking_x,
#'   blocking_y = blocking_y,
#'   controls_blocking = list(
#'     representation = "custom_matrix",
#'     ann = "kd",
#'     distance = "euclidean",
#'     seed = 1
#'   ),
#'   true_matches = data.frame(a = 1:5, b = 1:5)
#' )
#' result
#' @import data.table
#' @importFrom stats predict runif
#' @importFrom purrr partial
#' @importFrom FixedPoint FixedPoint
#' @importFrom blocking blocking
#' @export
mec_blocking <- function(
    A,
    B,
    variables,
    comparators = NULL,
    methods = NULL,
    blocking_x = NULL,
    blocking_y = NULL,
    blocking_variables = variables,
    blocking_sep = " ",
    controls_blocking = list(),
    min_training_pairs = NULL,
    min_training_nonmatches = NULL,
    block_sampling_seed = NULL,
    nonmatch_sample_size = NULL,
    nonmatch_sampling_seed = NULL,
    prob_ratio = "2",
    start_params = NULL,
    nonpar_hurdle = TRUE,
    fixed_method = "Newton",
    delta = 0.5,
    eps = 0.05,
    max_iter_em = 10,
    tol_em = 1,
    controls_nleqslv = list(),
    controls_kliep = control_kliep(),
    true_matches = NULL,
    keep_blocking_result = FALSE,
    keep_training_data = FALSE,
    verbose = FALSE) {

  stopifnot("`A` should be a data.frame or a data.table." =
              is.data.frame(A) | is.data.table(A))
  stopifnot("`B` should be a data.frame or a data.table." =
              is.data.frame(B) | is.data.table(B))
  stopifnot("`variables` should be a character vector." =
              is.character(variables))
  stopifnot("`variables` should contain at least one variable." =
              length(variables) > 0L)

  methods <- validate_methods(
    methods = methods,
    variables = variables,
    allowed_methods = c("binary", "continuous_parametric")
  )
  validate_choice(prob_ratio, c("1", "2"), "prob_ratio")
  validate_blocking_controls(controls_blocking)

  A <- data.table::copy(data.table::as.data.table(A))
  B <- data.table::copy(data.table::as.data.table(B))
  n_A_original <- nrow(A)
  n_B_original <- nrow(B)
  nonmatch_sample_size <- validate_nonmatch_sample_size(
    nonmatch_sample_size = nonmatch_sample_size,
    n_A = n_A_original,
    n_B = n_B_original
  )

  true_matches <- sanitize_true_matches(
    true_matches = true_matches,
    n_A = n_A_original,
    n_B = n_B_original,
    arg_name = "true_matches"
  )

  blocking_inputs <- prepare_blocking_inputs(
    A = A,
    B = B,
    blocking_x = blocking_x,
    blocking_y = blocking_y,
    blocking_variables = blocking_variables,
    blocking_sep = blocking_sep
  )

  if (verbose) {
    message("Running blocking::blocking().")
  }

  blocking_result <- run_blocking(
    blocking_x = blocking_inputs$x,
    blocking_y = blocking_inputs$y,
    controls_blocking = controls_blocking
  )
  blocking_table <- extract_blocking_result_table(blocking_result)
  reconstruction <- reconstruct_block_summary(
    blocking_table = blocking_table,
    n_A = n_A_original,
    n_B = n_B_original
  )
  block_summary <- reconstruction$block_summary
  excluded_records <- reconstruction$excluded_records

  if (verbose) {
    message(sprintf("Reconstructed %d disjoint blocks.", NROW(block_summary)))
  }

  data.table::set(A, j = "a", value = seq_len(n_A_original))
  data.table::set(B, j = "b", value = seq_len(n_B_original))
  stopifnot("There are no records with perfect agreement on the key variables.
            Please provide relevant datasets." =
              has_perfect_agreement(A, B, variables))

  preprocessing <- drop_constant_key_variables(
    A = A,
    B = B,
    variables = variables,
    comparators = comparators,
    methods = methods
  )
  A <- preprocessing$A
  B <- preprocessing$B
  variables <- preprocessing$variables
  comparators <- preprocessing$comparators
  methods <- preprocessing$methods
  constant_vars <- preprocessing$constant_vars

  for (var in constant_vars) {
    warning(paste("The variable", var, "has only one unique value and has been removed."))
  }

  exact_matches <- exact_match_pairs(A, B, variables)
  training_selection <- select_training_blocks(
    block_summary = block_summary,
    min_training_pairs = min_training_pairs,
    min_training_nonmatches = min_training_nonmatches,
    block_sampling_seed = block_sampling_seed
  )
  training_rule <- training_selection$training_rule
  training_blocks <- training_selection$training_blocks
  training_pairs <- make_block_pair_table(block_summary, training_blocks[["block"]])

  if (verbose) {
    message(sprintf(
      "Training MEC on %d pairs from %d block(s).",
      NROW(training_pairs),
      NROW(training_blocks)
    ))
  }

  training_vectors <- comparison_vectors(
    A = A,
    B = B,
    variables = variables,
    comparators = comparators,
    pairs = training_pairs
  )
  training_Omega <- training_vectors$Omega
  comparators <- training_vectors$comparators

  if (verbose) {
    message(sprintf(
      "Estimating nonmatch parameters from %d pair(s).",
      nonmatch_sample_size
    ))
  }

  nonmatch_params <- estimate_nonmatch_parameters(
    A = A,
    B = B,
    variables = variables,
    comparators = comparators,
    methods = methods,
    nonmatch_sample_size = nonmatch_sample_size,
    nonmatch_sampling_seed = nonmatch_sampling_seed,
    controls_nleqslv = controls_nleqslv
  )

  pooled_fit <- fit_mec_unsupervised_omega(
    A = A,
    B = B,
    variables = variables,
    comparators = comparators,
    methods = methods,
    Omega = training_Omega,
    exact_matches = exact_matches,
    start_params = start_params,
    nonmatch_params = nonmatch_params,
    nonpar_hurdle = nonpar_hurdle,
    delta = delta,
    eps = eps,
    max_iter_em = max_iter_em,
    tol_em = tol_em,
    controls_nleqslv = controls_nleqslv,
    controls_kliep = controls_kliep,
    context = "mec_blocking()"
  )
  pooled_model <- pooled_fit$model

  block_results <- lapply(seq_len(NROW(block_summary)), function(i) {
    current_summary <- block_summary[i]
    pair_table <- make_block_pair_table(current_summary)
    vectors <- comparison_vectors(
      A = A,
      B = B,
      variables = variables,
      comparators = comparators,
      pairs = pair_table
    )
    Omega_block <- score_mec_ratio(vectors$Omega, pooled_model)
    local_size <- estimate_local_mec_size(
      ratio = Omega_block[["ratio"]],
      n_pairs = NROW(Omega_block),
      n_A = current_summary[["n_A"]],
      n_B = current_summary[["n_B"]],
      fixed_method = fixed_method,
      prob_ratio = prob_ratio,
      prob_est = pooled_model$prob_est
    )
    data.table::set(Omega_block, j = "g_est", value = local_size$g_est)
    selection_summary <- summarize_mec_selection(
      a = Omega_block[["a"]],
      b = Omega_block[["b"]],
      ratio = Omega_block[["ratio"]],
      g_est = Omega_block[["g_est"]],
      n_M_est = local_size$n_M_est,
      duplicates_in_A = FALSE,
      set_construction = "size"
    )
    M_block <- Omega_block[selection_summary$selected_idx, c("a", "b", "block", "ratio", "g_est"), with = FALSE]
    block_estimate <- current_summary[, .(
      block,
      n_A,
      n_B,
      pair_count,
      nonmatches_min
    )]
    block_estimate[, n_M_est := local_size$n_M_est]
    block_estimate[, selected_pairs := NROW(M_block)]
    block_estimate[, flr_est := selection_summary$flr_est]
    block_estimate[, mmr_est := selection_summary$mmr_est]

    list(M_est = M_block, block_estimate = block_estimate)
  })

  M_est_full <- data.table::rbindlist(lapply(block_results, `[[`, "M_est"), use.names = TRUE, fill = TRUE)
  block_estimates <- data.table::rbindlist(lapply(block_results, `[[`, "block_estimate"), use.names = TRUE)
  n_M_est <- sum(block_estimates[["n_M_est"]])
  selected_g_est <- M_est_full[["g_est"]]
  flr_est <- if (NROW(M_est_full) == 0L) Inf else mean(1 - selected_g_est)
  mmr_est <- if (n_M_est == 0) {
    1
  } else {
    max(0, min(1, 1 - sum(selected_g_est) / n_M_est))
  }
  M_est <- M_est_full[, c("a", "b", "block", "ratio"), with = FALSE]

  block_pairs <- make_block_pair_table(block_summary)
  blocking_eval <- NULL
  eval_metrics <- NULL
  confusion <- NULL

  if (!is.null(true_matches)) {
    blocking_eval <- blocking_diagnostics(
      true_matches = true_matches,
      block_pairs = block_pairs,
      n_full_pairs = n_A_original * n_B_original
    )
    eval <- evaluation(M_est, true_matches, n_A_original * n_B_original)
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
      M_est = M_est,
      n_M_est = n_M_est,
      flr_est = flr_est,
      mmr_est = mmr_est,
      training_rule = training_rule,
      block_estimates = block_estimates,
      training_blocks = training_blocks,
      block_summary = block_summary,
      excluded_records = excluded_records,
      pooled_model = pooled_model,
      b_vars = pooled_model$b_vars,
      cpar_vars = pooled_model$cpar_vars,
      cnonpar_vars = pooled_model$cnonpar_vars,
      hm_vars = pooled_model$hm_vars,
      b_params = pooled_model$b_params,
      cpar_params = pooled_model$cpar_params,
      cnonpar_params = pooled_model$cnonpar_params,
      hm_params = pooled_model$hm_params,
      ratio_kliep = pooled_model$ratio_kliep,
      ratio_kliep_list = pooled_model$ratio_kliep_list,
      variables = variables,
      comparators = comparators,
      methods = methods,
      nonmatch_sample_size = nonmatch_params$sample_size,
      nonmatch_sampling_seed = nonmatch_sampling_seed,
      prob_ratio = prob_ratio,
      delta = delta,
      eps = eps,
      controls_nleqslv = controls_nleqslv,
      controls_blocking = controls_blocking,
      blocking_result = if (keep_blocking_result) blocking_result else NULL,
      training_Omega = if (keep_training_data) pooled_fit$Omega else NULL,
      blocking_eval = blocking_eval,
      eval_metrics = eval_metrics,
      confusion = confusion
    ),
    class = "mec_blocking"
  )
}
