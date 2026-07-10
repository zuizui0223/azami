#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ape)
  library(jsonlite)
  library(nlme)
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
  epsilon <- min(positive) * 1e-6
  pruned$edge.length[!is.finite(pruned$edge.length) | pruned$edge.length <= 0] <- epsilon
  data[[endpoint]] <- standardize(data[[endpoint]])
  for (predictor in predictors) data[[predictor]] <- standardize(data[[predictor]])
  rownames(data) <- data$taxon_key
  list(tree = pruned, data = data)
}

extract_rows <- function(fit, endpoint, endpoint_row, scenario, replicate, model, n_species, lambda = NA_real_) {
  coefficient_table <- summary(fit)$tTable
  rows <- lapply(predictors, function(predictor) {
    estimate <- as.numeric(coefficient_table[predictor, "Value"])
    standard_error <- as.numeric(coefficient_table[predictor, "Std.Error"])
    p_value <- as.numeric(coefficient_table[predictor, "p-value"])
    data.frame(
      analysis_tier = endpoint_row$analysis_tier,
      trait_group = endpoint_row$trait_group,
      endpoint = endpoint,
      endpoint_interpretation = endpoint_row$interpretation,
      scenario = scenario,
      replicate = replicate,
      covariance_model = model,
      n_species = n_species,
      pagel_lambda = lambda,
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
  phy <- prepared$tree
  formula <- as.formula(sprintf(
    "%s ~ %s",
    endpoint,
    paste(predictors, collapse = " + ")
  ))
  rows <- list()
  statuses <- list()
  for (model in models) {
    result <- tryCatch({
      if (model == "OLS_HC3_reference") {
        fit <- lm(formula, data = data)
        coefficient_table <- summary(fit)$coefficients
        model_rows <- lapply(predictors, function(predictor) {
          estimate <- as.numeric(coefficient_table[predictor, "Estimate"])
          standard_error <- as.numeric(coefficient_table[predictor, "Std. Error"])
          p_value <- as.numeric(coefficient_table[predictor, "Pr(>|t|)"])
          data.frame(
            analysis_tier = endpoint_row$analysis_tier,
            trait_group = endpoint_row$trait_group,
            endpoint = endpoint,
            endpoint_interpretation = endpoint_row$interpretation,
            scenario = scenario,
            replicate = replicate,
            covariance_model = model,
            n_species = nrow(data),
            pagel_lambda = NA_real_,
            predictor = predictor,
            estimate_standardized = estimate,
            standard_error = standard_error,
            ci95_low = estimate - 1.96 * standard_error,
            ci95_high = estimate + 1.96 * standard_error,
            p_value = p_value,
            stringsAsFactors = FALSE
          )
        })
        list(rows = do.call(rbind, model_rows), lambda = NA_real_)
      } else if (model == "PGLS_Brownian_fixed") {
        correlation <- ape::corBrownian(
          value = 1,
          phy = phy,
          form = ~ taxon_key
        )
        fit <- nlme::gls(
          formula,
          data = data,
          correlation = correlation,
          method = "ML",
          control = glsControl(opt = "optim", maxIter = 300, msMaxIter = 300)
        )
        list(
          rows = extract_rows(
            fit, endpoint, endpoint_row, scenario, replicate,
            model, nrow(data), NA_real_
          ),
          lambda = NA_real_
        )
      } else if (model == "PGLS_Pagel_lambda") {
        correlation <- ape::corPagel(
          value = 0.5,
          phy = phy,
          fixed = FALSE,
          form = ~ taxon_key
        )
        fit <- nlme::gls(
          formula,
          data = data,
          correlation = correlation,
          method = "ML",
          control = glsControl(opt = "optim", maxIter = 500, msMaxIter = 500)
        )
        lambda <- as.numeric(coef(fit$modelStruct$corStruct, unconstrained = FALSE)[1])
        list(
          rows = extract_rows(
            fit, endpoint, endpoint_row, scenario, replicate,
            model, nrow(data), lambda
          ),
          lambda = lambda
        )
      } else {
        stop(sprintf("Unknown model: %s", model))
      }
    }, error = function(error) error)
    if (inherits(result, "error")) {
      statuses[[length(statuses) + 1]] <- data.frame(
        analysis_tier = endpoint_row$analysis_tier,
        trait_group = endpoint_row$trait_group,
        endpoint = endpoint,
        scenario = scenario,
        replicate = replicate,
        covariance_model = model,
        status = "failed",
        n_species = nrow(data),
        pagel_lambda = NA_real_,
        message = sprintf("%s: %s", class(result)[1], conditionMessage(result)),
        stringsAsFactors = FALSE
      )
    } else {
      rows[[length(rows) + 1]] <- result$rows
      statuses[[length(statuses) + 1]] <- data.frame(
        analysis_tier = endpoint_row$analysis_tier,
        trait_group = endpoint_row$trait_group,
        endpoint = endpoint,
        scenario = scenario,
        replicate = replicate,
        covariance_model = model,
        status = "success",
        n_species = nrow(data),
        pagel_lambda = result$lambda,
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
  endpoint <- endpoint_row$species_variable
  if (!endpoint %in% names(species)) next

  deterministic_models <- c(
    "OLS_HC3_reference",
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
  if (endpoint_row$analysis_tier == "main") {
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

# Correct p-values within each fitted tree/model/tier family.
family_key <- interaction(
  coefficients$analysis_tier,
  coefficients$scenario,
  coefficients$replicate,
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
    data.frame(
      analysis_tier = group$analysis_tier[1],
      trait_group = group$trait_group[1],
      endpoint = group$endpoint[1],
      predictor = group$predictor[1],
      n_successful_random_trees = nrow(group),
      estimate_median = median(estimates, na.rm = TRUE),
      estimate_q025 = quantile(estimates, 0.025, na.rm = TRUE, names = FALSE),
      estimate_q975 = quantile(estimates, 0.975, na.rm = TRUE, names = FALSE),
      fraction_positive = mean(estimates > 0, na.rm = TRUE),
      fraction_negative = mean(estimates < 0, na.rm = TRUE),
      sign_stability = max(
        mean(estimates > 0, na.rm = TRUE),
        mean(estimates < 0, na.rm = TRUE)
      ),
      fraction_nominal_p_lt_0_05 = mean(group$p_value < 0.05, na.rm = TRUE),
      fraction_fdr_signal = mean(
        group$historical_sensitivity_signal_fdr_0_05,
        na.rm = TRUE
      ),
      pagel_lambda_median = median(lambdas, na.rm = TRUE),
      pagel_lambda_q025 = quantile(lambdas, 0.025, na.rm = TRUE, names = FALSE),
      pagel_lambda_q975 = quantile(lambdas, 0.975, na.rm = TRUE, names = FALSE),
      random_tree_interval_excludes_zero = (
        quantile(estimates, 0.025, na.rm = TRUE, names = FALSE) > 0 ||
          quantile(estimates, 0.975, na.rm = TRUE, names = FALSE) < 0
      ),
      stringsAsFactors = FALSE
    )
  }))
  random_summary$stable_direction_across_random_trees <- (
    random_summary$sign_stability >= 0.90
  )
  random_summary$historical_sensitivity_robust_candidate <- (
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

lambda_summary <- aggregate(
  pagel_lambda ~ analysis_tier + trait_group + endpoint + scenario,
  data = statuses[
    statuses$status == "success" &
      statuses$covariance_model == "PGLS_Pagel_lambda" &
      is.finite(statuses$pagel_lambda),
    ,
    drop = FALSE
  ],
  FUN = function(values) c(
    n = length(values),
    median = median(values),
    q025 = quantile(values, 0.025, names = FALSE),
    q975 = quantile(values, 0.975, names = FALSE)
  )
)
if (nrow(lambda_summary)) {
  expanded <- data.frame(
    lambda_summary[, c("analysis_tier", "trait_group", "endpoint", "scenario")],
    n = lambda_summary$pagel_lambda[, "n"],
    lambda_median = lambda_summary$pagel_lambda[, "median"],
    lambda_q025 = lambda_summary$pagel_lambda[, "q025"],
    lambda_q975 = lambda_summary$pagel_lambda[, "q975"],
    row.names = NULL
  )
  write.csv(
    expanded,
    file.path(out_dir, "historical_constraint_pagel_lambda_summary.csv"),
    row.names = FALSE,
    na = ""
  )
}

report <- list(
  created_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  n_registered_endpoints = nrow(registry),
  n_main_endpoints = sum(registry$analysis_tier == "main"),
  n_auxiliary_endpoints = sum(registry$analysis_tier == "auxiliary"),
  n_successful_model_fits = sum(statuses$status == "success"),
  n_failed_model_fits = sum(statuses$status == "failed"),
  n_coefficient_rows = nrow(coefficients),
  n_random_tree_summary_rows = nrow(random_summary),
  n_historical_sensitivity_robust_candidates = if (nrow(random_summary)) {
    sum(random_summary$historical_sensitivity_robust_candidate)
  } else 0L,
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
      "missing from the backbone. Robust candidates require stable direction",
      "and intervals excluding zero across randomized graftings."
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
