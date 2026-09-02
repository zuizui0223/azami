#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ape)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
value_after <- function(flag, default = NULL) {
  index <- match(flag, args)
  if (is.na(index)) return(default)
  if (index == length(args)) stop(sprintf("Missing value after %s", flag))
  args[[index + 1]]
}

species_path <- value_after("--species")
registry_path <- value_after("--registry")
tree_dir <- value_after("--tree-dir")
out_dir <- value_after("--out-dir")
min_species <- as.integer(value_after("--min-species", "30"))
if (any(vapply(
  list(species_path, registry_path, tree_dir, out_dir),
  is.null,
  logical(1)
))) stop("Required: --species --registry --tree-dir --out-dir")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

species <- read.csv(species_path, stringsAsFactors = FALSE, check.names = FALSE)
registry <- read.csv(registry_path, stringsAsFactors = FALSE, check.names = FALSE)
required_registry <- c(
  "analysis_tier", "trait_group", "species_variable", "interpretation"
)
if (!all(required_registry %in% names(registry))) {
  stop("Endpoint registry is incomplete")
}
if (!"taxon_name" %in% names(species)) stop("Species table lacks taxon_name")
species$taxon_key <- gsub(" ", "_", trimws(species$taxon_name), fixed = TRUE)
if (anyDuplicated(species$taxon_key)) stop("taxon_key is not unique")

predictors <- c(
  "env_chelsa_bio01_species_median",
  "env_chelsa_bio04_species_median",
  "env_chelsa_bio12_species_median",
  "env_chelsa_bio15_species_median"
)
if (!all(predictors %in% names(species))) stop("Species table lacks CHELSA predictors")

scenario1 <- read.tree(file.path(tree_dir, "gbotb_lcvp_scenario1.tre"))
scenario3 <- read.tree(file.path(tree_dir, "gbotb_lcvp_scenario3.tre"))
scenario2 <- read.tree(file.path(tree_dir, "gbotb_lcvp_scenario2_randomized.trees"))
if (!inherits(scenario2, "multiPhylo")) class(scenario2) <- "multiPhylo"
n_random_expected <- length(scenario2)
if (n_random_expected < 1L) stop("No randomized historical trees were loaded")

standardize <- function(values) {
  values <- as.numeric(values)
  scale_value <- sd(values, na.rm = TRUE)
  if (!is.finite(scale_value) || scale_value <= 0) stop("Variable has zero variance")
  (values - mean(values, na.rm = TRUE)) / scale_value
}

prepare_tree_data <- function(tree, endpoint) {
  columns <- c("taxon_name", "taxon_key", endpoint, predictors)
  data <- species[, columns, drop = FALSE]
  for (column in c(endpoint, predictors)) {
    data[[column]] <- suppressWarnings(as.numeric(data[[column]]))
  }
  data <- data[complete.cases(data), , drop = FALSE]
  common <- intersect(tree$tip.label, data$taxon_key)
  if (length(common) < min_species) {
    stop(sprintf("Only %d complete species remain", length(common)))
  }
  drop <- setdiff(tree$tip.label, common)
  pruned <- if (length(drop)) drop.tip(tree, drop) else tree
  data <- data[match(pruned$tip.label, data$taxon_key), , drop = FALSE]
  if (!identical(data$taxon_key, pruned$tip.label)) stop("Tree/data order mismatch")
  if (is.null(pruned$edge.length)) stop("PGLS tree has no branch lengths")
  positive <- pruned$edge.length[is.finite(pruned$edge.length) & pruned$edge.length > 0]
  if (!length(positive)) stop("PGLS tree has no positive branch lengths")
  epsilon <- min(positive) * 1e-7
  pruned$edge.length[!is.finite(pruned$edge.length) | pruned$edge.length <= 0] <- epsilon
  data[[endpoint]] <- standardize(data[[endpoint]])
  for (predictor in predictors) data[[predictor]] <- standardize(data[[predictor]])
  rownames(data) <- data$taxon_key

  covariance <- ape::vcv.phylo(pruned, corr = FALSE)
  covariance <- covariance[data$taxon_key, data$taxon_key, drop = FALSE]
  covariance <- (covariance + t(covariance)) / 2
  if (any(!is.finite(covariance))) stop("Phylogenetic covariance contains non-finite values")
  diagonal_mean <- mean(diag(covariance))
  if (!is.finite(diagonal_mean) || diagonal_mean <= 0) {
    stop("Phylogenetic covariance has invalid diagonal")
  }
  covariance <- covariance / diagonal_mean
  list(tree = pruned, data = data, covariance = covariance)
}

safe_chol <- function(matrix) {
  scale <- mean(diag(matrix))
  if (!is.finite(scale) || scale <= 0) scale <- 1
  multipliers <- c(0, 1e-12, 1e-10, 1e-8, 1e-6, 1e-5, 1e-4)
  for (multiplier in multipliers) {
    candidate <- matrix
    if (multiplier > 0) {
      diag(candidate) <- diag(candidate) + multiplier * scale
    }
    factor <- tryCatch(chol(candidate), error = function(error) NULL)
    if (!is.null(factor)) {
      return(list(factor = factor, jitter = multiplier * scale))
    }
  }
  stop("Phylogenetic covariance is not positive definite after adaptive jitter")
}

solve_from_chol <- function(chol_factor, right_hand_side) {
  backsolve(
    chol_factor,
    forwardsolve(t(chol_factor), right_hand_side)
  )
}

fit_gls_at_lambda <- function(response, design, base_covariance, lambda) {
  if (!is.finite(lambda) || lambda < 0 || lambda > 1) {
    stop("Pagel lambda must be in [0,1]")
  }
  covariance <- lambda * base_covariance
  diag(covariance) <- diag(base_covariance)
  covariance <- (covariance + t(covariance)) / 2
  decomposition <- safe_chol(covariance)
  factor <- decomposition$factor

  inverse_design <- solve_from_chol(factor, design)
  inverse_response <- as.vector(solve_from_chol(factor, response))
  information <- crossprod(design, inverse_design)
  information_inverse <- tryCatch(
    solve(information),
    error = function(error) qr.solve(information, diag(nrow(information)))
  )
  beta <- as.vector(information_inverse %*% crossprod(design, inverse_response))
  names(beta) <- colnames(design)
  residual <- response - as.vector(design %*% beta)
  inverse_residual <- as.vector(solve_from_chol(factor, residual))
  residual_ss <- as.numeric(crossprod(residual, inverse_residual))
  n <- length(response)
  p <- ncol(design)
  if (!is.finite(residual_ss) || residual_ss <= 0 || n <= p) {
    stop("Invalid residual variance in matrix PGLS")
  }
  sigma2_ml <- residual_ss / n
  sigma2_unbiased <- residual_ss / (n - p)
  log_determinant <- 2 * sum(log(diag(factor)))
  log_likelihood <- -0.5 * (
    n * (log(2 * pi) + 1 + log(sigma2_ml)) + log_determinant
  )
  coefficient_covariance <- sigma2_unbiased * information_inverse
  standard_error <- sqrt(pmax(diag(coefficient_covariance), 0))
  names(standard_error) <- colnames(design)
  statistic <- beta / standard_error
  p_value <- 2 * pt(-abs(statistic), df = n - p)
  list(
    coefficients = beta,
    standard_error = standard_error,
    p_value = p_value,
    log_likelihood = log_likelihood,
    lambda = lambda,
    covariance_jitter = decomposition$jitter,
    n = n,
    df = n - p
  )
}

fit_bounded_pagel <- function(response, design, covariance) {
  evaluate <- function(lambda) {
    tryCatch(
      fit_gls_at_lambda(response, design, covariance, lambda),
      error = function(error) NULL
    )
  }
  grid <- seq(0, 1, length.out = 21)
  grid_fits <- lapply(grid, evaluate)
  grid_loglik <- vapply(
    grid_fits,
    function(fit) if (is.null(fit)) -Inf else fit$log_likelihood,
    numeric(1)
  )
  if (!any(is.finite(grid_loglik))) stop("No bounded Pagel-lambda grid fit succeeded")
  best_index <- which.max(grid_loglik)
  lower <- grid[max(1, best_index - 1)]
  upper <- grid[min(length(grid), best_index + 1)]
  optimized <- tryCatch(
    optimize(
      function(lambda) {
        fit <- evaluate(lambda)
        if (is.null(fit)) return(.Machine$double.xmax / 100)
        -fit$log_likelihood
      },
      interval = c(lower, upper),
      tol = 1e-6
    ),
    error = function(error) NULL
  )
  candidates <- grid_fits[!vapply(grid_fits, is.null, logical(1))]
  if (!is.null(optimized)) {
    optimized_fit <- evaluate(optimized$minimum)
    if (!is.null(optimized_fit)) candidates[[length(candidates) + 1]] <- optimized_fit
  }
  boundary_zero <- evaluate(0)
  boundary_one <- evaluate(1)
  if (!is.null(boundary_zero)) candidates[[length(candidates) + 1]] <- boundary_zero
  if (!is.null(boundary_one)) candidates[[length(candidates) + 1]] <- boundary_one
  log_likelihoods <- vapply(candidates, `[[`, numeric(1), "log_likelihood")
  candidates[[which.max(log_likelihoods)]]
}

model_matrix_for <- function(data, endpoint) {
  formula <- as.formula(sprintf(
    "%s ~ %s",
    endpoint,
    paste(predictors, collapse = " + ")
  ))
  list(
    formula = formula,
    response = as.numeric(data[[endpoint]]),
    design = model.matrix(formula, data = data)
  )
}

rows_from_matrix_fit <- function(
  fit,
  endpoint,
  endpoint_row,
  scenario,
  replicate,
  model
) {
  rows <- lapply(predictors, function(predictor) {
    estimate <- as.numeric(fit$coefficients[predictor])
    standard_error <- as.numeric(fit$standard_error[predictor])
    p_value <- as.numeric(fit$p_value[predictor])
    data.frame(
      analysis_tier = as.character(endpoint_row$analysis_tier),
      trait_group = as.character(endpoint_row$trait_group),
      endpoint = endpoint,
      endpoint_interpretation = as.character(endpoint_row$interpretation),
      scenario = scenario,
      replicate = replicate,
      covariance_model = model,
      n_species = fit$n,
      pagel_lambda = if (model == "PGLS_Pagel_lambda") fit$lambda else NA_real_,
      covariance_jitter = fit$covariance_jitter,
      predictor = predictor,
      estimate_standardized = estimate,
      standard_error = standard_error,
      ci95_low = estimate - 1.96 * standard_error,
      ci95_high = estimate + 1.96 * standard_error,
      p_value = p_value,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

fit_models <- function(tree, endpoint, endpoint_row, scenario, replicate, models) {
  prepared <- prepare_tree_data(tree, endpoint)
  data <- prepared$data
  matrices <- model_matrix_for(data, endpoint)
  rows <- list()
  statuses <- list()
  for (model in models) {
    result <- tryCatch({
      if (model == "OLS_standard_reference") {
        lm_fit <- lm(matrices$formula, data = data)
        coefficient_table <- summary(lm_fit)$coefficients
        model_rows <- lapply(predictors, function(predictor) {
          estimate <- as.numeric(coefficient_table[predictor, "Estimate"])
          standard_error <- as.numeric(coefficient_table[predictor, "Std. Error"])
          p_value <- as.numeric(coefficient_table[predictor, "Pr(>|t|)"])
          data.frame(
            analysis_tier = as.character(endpoint_row$analysis_tier),
            trait_group = as.character(endpoint_row$trait_group),
            endpoint = endpoint,
            endpoint_interpretation = as.character(endpoint_row$interpretation),
            scenario = scenario,
            replicate = replicate,
            covariance_model = model,
            n_species = nrow(data),
            pagel_lambda = NA_real_,
            covariance_jitter = 0,
            predictor = predictor,
            estimate_standardized = estimate,
            standard_error = standard_error,
            ci95_low = estimate - 1.96 * standard_error,
            ci95_high = estimate + 1.96 * standard_error,
            p_value = p_value,
            stringsAsFactors = FALSE
          )
        })
        list(
          rows = do.call(rbind, model_rows),
          lambda = NA_real_,
          jitter = 0
        )
      } else if (model == "PGLS_Brownian_fixed") {
        fit <- fit_gls_at_lambda(
          matrices$response,
          matrices$design,
          prepared$covariance,
          lambda = 1
        )
        list(
          rows = rows_from_matrix_fit(
            fit, endpoint, endpoint_row, scenario, replicate, model
          ),
          lambda = NA_real_,
          jitter = fit$covariance_jitter
        )
      } else if (model == "PGLS_Pagel_lambda") {
        fit <- fit_bounded_pagel(
          matrices$response,
          matrices$design,
          prepared$covariance
        )
        list(
          rows = rows_from_matrix_fit(
            fit, endpoint, endpoint_row, scenario, replicate, model
          ),
          lambda = fit$lambda,
          jitter = fit$covariance_jitter
        )
      } else {
        stop(sprintf("Unknown model: %s", model))
      }
    }, error = function(error) error)

    if (inherits(result, "error")) {
      statuses[[length(statuses) + 1]] <- data.frame(
        analysis_tier = as.character(endpoint_row$analysis_tier),
        trait_group = as.character(endpoint_row$trait_group),
        endpoint = endpoint,
        scenario = scenario,
        replicate = replicate,
        covariance_model = model,
        status = "failed",
        n_species = nrow(data),
        pagel_lambda = NA_real_,
        covariance_jitter = NA_real_,
        message = sprintf("%s: %s", class(result)[1], conditionMessage(result)),
        stringsAsFactors = FALSE
      )
    } else {
      rows[[length(rows) + 1]] <- result$rows
      statuses[[length(statuses) + 1]] <- data.frame(
        analysis_tier = as.character(endpoint_row$analysis_tier),
        trait_group = as.character(endpoint_row$trait_group),
        endpoint = endpoint,
        scenario = scenario,
        replicate = replicate,
        covariance_model = model,
        status = "success",
        n_species = nrow(data),
        pagel_lambda = result$lambda,
        covariance_jitter = result$jitter,
        message = "",
        stringsAsFactors = FALSE
      )
    }
  }
  list(
    rows = if (length(rows)) do.call(rbind, rows) else NULL,
    statuses = do.call(rbind, statuses)
  )
}

all_rows <- list()
all_status <- list()
for (index in seq_len(nrow(registry))) {
  endpoint_row <- registry[index, , drop = FALSE]
  endpoint <- as.character(endpoint_row$species_variable)
  if (!endpoint %in% names(species)) next

  deterministic_models <- c(
    "OLS_standard_reference",
    "PGLS_Brownian_fixed",
    "PGLS_Pagel_lambda"
  )
  for (scenario_item in list(
    list(tree = scenario1, label = "GBOTB_LCVP_S1", replicate = NA_integer_),
    list(tree = scenario3, label = "GBOTB_LCVP_S3", replicate = NA_integer_)
  )) {
    fitted <- fit_models(
      scenario_item$tree,
      endpoint,
      endpoint_row,
      scenario_item$label,
      scenario_item$replicate,
      deterministic_models
    )
    if (!is.null(fitted$rows)) all_rows[[length(all_rows) + 1]] <- fitted$rows
    all_status[[length(all_status) + 1]] <- fitted$statuses
  }

  # Randomized genus-grafting uncertainty is propagated for main endpoints only.
  if (as.character(endpoint_row$analysis_tier) == "main") {
    for (replicate in seq_along(scenario2)) {
      fitted <- fit_models(
        scenario2[[replicate]],
        endpoint,
        endpoint_row,
        "GBOTB_LCVP_S2_randomized",
        replicate,
        c("PGLS_Pagel_lambda")
      )
      if (!is.null(fitted$rows)) all_rows[[length(all_rows) + 1]] <- fitted$rows
      all_status[[length(all_status) + 1]] <- fitted$statuses
    }
  }
}

coefficients <- if (length(all_rows)) do.call(rbind, all_rows) else data.frame()
statuses <- if (length(all_status)) do.call(rbind, all_status) else data.frame()
if (!nrow(coefficients)) stop("No historical-constraint model succeeded")

replicate_key <- ifelse(
  is.na(coefficients$replicate),
  "deterministic",
  as.character(coefficients$replicate)
)
family_key <- interaction(
  coefficients$analysis_tier,
  coefficients$scenario,
  replicate_key,
  coefficients$covariance_model,
  drop = TRUE,
  lex.order = TRUE
)
coefficients$p_fdr_within_tree_model_tier <- ave(
  coefficients$p_value,
  family_key,
  FUN = function(values) p.adjust(values, method = "BH")
)
coefficients$historical_sensitivity_signal_fdr_0_05 <- (
  coefficients$p_fdr_within_tree_model_tier < 0.05
)
write.csv(
  coefficients,
  file.path(out_dir, "historical_constraint_model_coefficients.csv"),
  row.names = FALSE,
  na = ""
)
write.csv(
  statuses,
  file.path(out_dir, "historical_constraint_model_status.csv"),
  row.names = FALSE,
  na = ""
)

random <- coefficients[
  coefficients$scenario == "GBOTB_LCVP_S2_randomized" &
    coefficients$covariance_model == "PGLS_Pagel_lambda",
  ,
  drop = FALSE
]
random_summary <- data.frame()
if (nrow(random)) {
  split_key <- interaction(
    random$analysis_tier,
    random$trait_group,
    random$endpoint,
    random$predictor,
    drop = TRUE,
    lex.order = TRUE
  )
  random_summary <- do.call(rbind, lapply(split(random, split_key), function(group) {
    estimates <- group$estimate_standardized
    lambdas <- group$pagel_lambda
    estimate_q025 <- quantile(estimates, 0.025, na.rm = TRUE, names = FALSE)
    estimate_q975 <- quantile(estimates, 0.975, na.rm = TRUE, names = FALSE)
    n_successful <- nrow(group)
    positive_fraction <- mean(estimates > 0, na.rm = TRUE)
    negative_fraction <- mean(estimates < 0, na.rm = TRUE)
    data.frame(
      analysis_tier = group$analysis_tier[1],
      trait_group = group$trait_group[1],
      endpoint = group$endpoint[1],
      predictor = group$predictor[1],
      n_expected_random_trees = n_random_expected,
      n_successful_random_trees = n_successful,
      complete_random_tree_support = n_successful == n_random_expected,
      estimate_median = median(estimates, na.rm = TRUE),
      estimate_q025 = estimate_q025,
      estimate_q975 = estimate_q975,
      fraction_positive = positive_fraction,
      fraction_negative = negative_fraction,
      sign_stability = max(positive_fraction, negative_fraction),
      fraction_nominal_p_lt_0_05 = mean(group$p_value < 0.05, na.rm = TRUE),
      fraction_fdr_signal = mean(
        group$historical_sensitivity_signal_fdr_0_05,
        na.rm = TRUE
      ),
      pagel_lambda_median = median(lambdas, na.rm = TRUE),
      pagel_lambda_q025 = quantile(lambdas, 0.025, na.rm = TRUE, names = FALSE),
      pagel_lambda_q975 = quantile(lambdas, 0.975, na.rm = TRUE, names = FALSE),
      random_tree_interval_excludes_zero = (
        estimate_q025 > 0 || estimate_q975 < 0
      ),
      stringsAsFactors = FALSE
    )
  }))
  random_summary$stable_direction_across_random_trees <- (
    random_summary$sign_stability >= 0.90
  )
  random_summary$historical_sensitivity_robust_candidate <- (
    random_summary$complete_random_tree_support &
      random_summary$stable_direction_across_random_trees &
      random_summary$random_tree_interval_excludes_zero &
      random_summary$fraction_nominal_p_lt_0_05 >= 0.80
  )
  write.csv(
    random_summary,
    file.path(out_dir, "historical_constraint_random_tree_summary.csv"),
    row.names = FALSE,
    na = ""
  )
}

lambda_rows <- statuses[
  statuses$status == "success" &
    statuses$covariance_model == "PGLS_Pagel_lambda" &
    is.finite(statuses$pagel_lambda),
  ,
  drop = FALSE
]
lambda_summary <- data.frame()
if (nrow(lambda_rows)) {
  lambda_key <- interaction(
    lambda_rows$analysis_tier,
    lambda_rows$trait_group,
    lambda_rows$endpoint,
    lambda_rows$scenario,
    drop = TRUE,
    lex.order = TRUE
  )
  lambda_summary <- do.call(rbind, lapply(split(lambda_rows, lambda_key), function(group) {
    values <- group$pagel_lambda
    data.frame(
      analysis_tier = group$analysis_tier[1],
      trait_group = group$trait_group[1],
      endpoint = group$endpoint[1],
      scenario = group$scenario[1],
      n = length(values),
      lambda_median = median(values),
      lambda_q025 = quantile(values, 0.025, names = FALSE),
      lambda_q975 = quantile(values, 0.975, names = FALSE),
      stringsAsFactors = FALSE
    )
  }))
  write.csv(
    lambda_summary,
    file.path(out_dir, "historical_constraint_pagel_lambda_summary.csv"),
    row.names = FALSE,
    na = ""
  )
}

expected_fits <- nrow(registry) * 2 * 3 + sum(registry$analysis_tier == "main") * n_random_expected
report <- list(
  created_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  n_registered_endpoints = nrow(registry),
  n_main_endpoints = sum(registry$analysis_tier == "main"),
  n_auxiliary_endpoints = sum(registry$analysis_tier == "auxiliary"),
  n_random_trees = n_random_expected,
  n_expected_model_fits = expected_fits,
  n_successful_model_fits = sum(statuses$status == "success"),
  n_failed_model_fits = sum(statuses$status == "failed"),
  n_coefficient_rows = nrow(coefficients),
  n_random_tree_summary_rows = nrow(random_summary),
  n_historical_sensitivity_robust_candidates = if (nrow(random_summary)) {
    sum(random_summary$historical_sensitivity_robust_candidate)
  } else 0L,
  lambda_estimation = paste(
    "Matrix PGLS with maximum-likelihood Pagel lambda explicitly bounded to [0,1];",
    "adaptive diagonal jitter is reported and used only for numerical positive-definiteness."
  ),
  interpretation = list(
    within_species_primary = paste(
      "The earlier within-species models remain the least tree-dependent",
      "historical control because they remove between-species differences."
    ),
    pgls_role = paste(
      "These PGLS models are sensitivity analyses over dated mega-tree",
      "grafting scenarios, not claims that GBOTB is the true Cirsium species tree."
    ),
    random_tree_role = paste(
      "Scenario-2 random trees propagate uncertainty in the placement of taxa",
      "missing from the backbone. Robust candidates require complete support",
      "across all random trees, stable direction and intervals excluding zero."
    ),
    unresolved_biology = paste(
      "Hybridization, chloroplast capture and allopolyploid parentage are not",
      "represented by a bifurcating mega-tree and require exclusion/evidence",
      "sensitivities when those data become available."
    )
  )
)
write_json(
  report,
  file.path(out_dir, "historical_constraint_model_report.json"),
  pretty = TRUE,
  auto_unbox = TRUE,
  null = "null"
)
cat(toJSON(report, pretty = TRUE, auto_unbox = TRUE), "\n")
