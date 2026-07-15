#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ape)
  library(jsonlite)
  library(phytools)
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
min_species <- as.integer(value_after("--min-species", "20"))
k_nsim <- as.integer(value_after("--k-nsim", "999"))
if (any(vapply(list(species_path, registry_path, tree_dir, out_dir), is.null, logical(1)))) {
  stop("Required: --species --registry --tree-dir --out-dir")
}
if (!is.finite(min_species) || min_species < 10L) stop("--min-species must be >= 10")
if (!is.finite(k_nsim) || k_nsim < 0L) stop("--k-nsim must be non-negative")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

species <- read.csv(species_path, stringsAsFactors = FALSE, check.names = FALSE)
registry <- read.csv(registry_path, stringsAsFactors = FALSE, check.names = FALSE)
audit <- read.csv(file.path(tree_dir, "gbotb_lcvp_backbone_audit.csv"), stringsAsFactors = FALSE)
required_registry <- c("analysis_tier", "trait_group", "species_variable", "interpretation")
if (!all(required_registry %in% names(registry))) stop("Endpoint registry is incomplete")
if (!"taxon_name" %in% names(species)) stop("Species table lacks taxon_name")
if (!all(c("taxon_name", "directly_present_in_GBOTB_extended_LCVP") %in% names(audit))) {
  stop("Backbone audit is incomplete")
}

species$taxon_key <- gsub(" ", "_", trimws(species$taxon_name), fixed = TRUE)
audit$taxon_key <- gsub(" ", "_", trimws(audit$taxon_name), fixed = TRUE)
direct_keys <- audit$taxon_key[as.logical(audit$directly_present_in_GBOTB_extended_LCVP)]
if (length(direct_keys) != 54L) stop(sprintf("Expected 54 direct tips, found %d", length(direct_keys)))

endpoints <- unique(registry$species_variable[nzchar(registry$species_variable)])
endpoints <- endpoints[endpoints %in% names(species)]
if (!length(endpoints)) stop("No registered species endpoints are present")

scenario1 <- read.tree(file.path(tree_dir, "gbotb_lcvp_scenario1.tre"))
scenario3 <- read.tree(file.path(tree_dir, "gbotb_lcvp_scenario3.tre"))
scenario2 <- read.tree(file.path(tree_dir, "gbotb_lcvp_scenario2_randomized.trees"))
if (!inherits(scenario2, "multiPhylo")) class(scenario2) <- "multiPhylo"
if (length(scenario2) < 1L) stop("No randomized trees found")

repair_tree <- function(tree) {
  if (is.null(tree$edge.length)) stop("Tree lacks branch lengths")
  positive <- tree$edge.length[is.finite(tree$edge.length) & tree$edge.length > 0]
  if (!length(positive)) stop("Tree lacks positive branch lengths")
  epsilon <- min(positive) * 1e-7
  tree$edge.length[!is.finite(tree$edge.length) | tree$edge.length <= 0] <- epsilon
  tree
}

prepare_trait <- function(tree, endpoint, allowed_keys = NULL) {
  values <- suppressWarnings(as.numeric(species[[endpoint]]))
  names(values) <- species$taxon_key
  values <- values[is.finite(values)]
  common <- intersect(tree$tip.label, names(values))
  if (!is.null(allowed_keys)) common <- intersect(common, allowed_keys)
  if (length(common) < min_species) return(NULL)
  pruned <- drop.tip(tree, setdiff(tree$tip.label, common))
  pruned <- repair_tree(pruned)
  x <- values[pruned$tip.label]
  if (length(unique(x)) < 3L || !is.finite(sd(x)) || sd(x) <= 0) return(NULL)
  list(tree = pruned, x = x)
}

extract_phylosig <- function(fit, estimate_name, run_test) {
  if (is.list(fit)) {
    if (is.null(fit[[estimate_name]])) {
      stop(sprintf("phylosig list result lacks %s", estimate_name))
    }
    estimate <- as.numeric(fit[[estimate_name]])[[1]]
    p_value <- if (run_test && !is.null(fit[["P"]])) as.numeric(fit[["P"]])[[1]] else NA_real_
  } else if (is.atomic(fit) && length(fit) >= 1L) {
    estimate <- as.numeric(fit)[[1]]
    p_value <- NA_real_
  } else {
    stop(sprintf("Unexpected phylosig return type for %s: %s", estimate_name, paste(class(fit), collapse = "/")))
  }
  if (!is.finite(estimate)) stop(sprintf("Non-finite %s estimate", estimate_name))
  list(estimate = estimate, p_value = p_value)
}

# Guard both documented return shapes before the expensive tree loop.
stopifnot(identical(extract_phylosig(0.5, "K", FALSE)$estimate, 0.5))
stopifnot(identical(extract_phylosig(list(K = 0.5, P = 0.1), "K", TRUE)$p_value, 0.1))
stopifnot(identical(extract_phylosig(list(lambda = 0.7, P = 0.2), "lambda", TRUE)$estimate, 0.7))

safe_signal <- function(tree, x, run_test) {
  k_fit <- tryCatch(
    phytools::phylosig(tree, x, method = "K", test = run_test, nsim = if (run_test) k_nsim else 0),
    error = function(error) error
  )
  lambda_fit <- tryCatch(
    phytools::phylosig(tree, x, method = "lambda", test = run_test),
    error = function(error) error
  )
  if (inherits(k_fit, "error") || inherits(lambda_fit, "error")) {
    message <- paste(
      if (inherits(k_fit, "error")) conditionMessage(k_fit) else "",
      if (inherits(lambda_fit, "error")) conditionMessage(lambda_fit) else ""
    )
    return(list(status = "failed", error = trimws(message)))
  }
  parsed <- tryCatch({
    k <- extract_phylosig(k_fit, "K", run_test)
    lambda <- extract_phylosig(lambda_fit, "lambda", run_test)
    list(K = k$estimate, K_p_value = k$p_value, lambda = lambda$estimate, lambda_p_value = lambda$p_value)
  }, error = function(error) error)
  if (inherits(parsed, "error")) {
    return(list(status = "failed", error = conditionMessage(parsed)))
  }
  c(list(status = "success"), parsed, list(error = ""))
}

registry_lookup <- registry[match(endpoints, registry$species_variable), , drop = FALSE]
rows <- list()
add_fit <- function(endpoint, endpoint_row, tree, scenario, replicate, information_tier, allowed_keys, run_test) {
  prepared <- prepare_trait(tree, endpoint, allowed_keys)
  if (is.null(prepared)) {
    rows[[length(rows) + 1L]] <<- data.frame(
      analysis_tier = endpoint_row$analysis_tier, trait_group = endpoint_row$trait_group,
      endpoint = endpoint, interpretation = endpoint_row$interpretation,
      scenario = scenario, replicate = replicate, information_tier = information_tier,
      n_species = 0L, K = NA_real_, K_p_value = NA_real_, lambda = NA_real_,
      lambda_p_value = NA_real_, status = "insufficient_data", error = "",
      stringsAsFactors = FALSE
    )
    return(invisible(NULL))
  }
  fit <- safe_signal(prepared$tree, prepared$x, run_test)
  rows[[length(rows) + 1L]] <<- data.frame(
    analysis_tier = endpoint_row$analysis_tier, trait_group = endpoint_row$trait_group,
    endpoint = endpoint, interpretation = endpoint_row$interpretation,
    scenario = scenario, replicate = replicate, information_tier = information_tier,
    n_species = length(prepared$x),
    K = if (fit$status == "success") fit$K else NA_real_,
    K_p_value = if (fit$status == "success") fit$K_p_value else NA_real_,
    lambda = if (fit$status == "success") fit$lambda else NA_real_,
    lambda_p_value = if (fit$status == "success") fit$lambda_p_value else NA_real_,
    status = fit$status, error = fit$error, stringsAsFactors = FALSE
  )
}

for (i in seq_along(endpoints)) {
  endpoint <- endpoints[[i]]
  endpoint_row <- registry_lookup[i, , drop = FALSE]
  add_fit(endpoint, endpoint_row, scenario1, "S1", NA_integer_, "all_tips", NULL, TRUE)
  add_fit(endpoint, endpoint_row, scenario3, "S3", NA_integer_, "all_tips", NULL, TRUE)
  add_fit(endpoint, endpoint_row, scenario1, "S1", NA_integer_, "direct_backbone_only", direct_keys, TRUE)
  for (j in seq_along(scenario2)) {
    add_fit(endpoint, endpoint_row, scenario2[[j]], "S2_random", j, "all_tips", NULL, FALSE)
  }
}

results <- do.call(rbind, rows)
write.csv(results, file.path(out_dir, "phylogenetic_signal_all_tree_fits.csv"), row.names = FALSE, na = "")
successful <- results[results$status == "success", , drop = FALSE]
random <- successful[successful$scenario == "S2_random" & successful$information_tier == "all_tips", , drop = FALSE]
quantile_safe <- function(x, probability) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  as.numeric(quantile(x, probability, names = FALSE, type = 8))
}
summary_rows <- lapply(split(random, random$endpoint), function(group) {
  data.frame(
    analysis_tier = group$analysis_tier[[1]], trait_group = group$trait_group[[1]],
    endpoint = group$endpoint[[1]], interpretation = group$interpretation[[1]],
    n_successful_random_trees = nrow(group), n_species_min = min(group$n_species),
    n_species_max = max(group$n_species), K_median = median(group$K, na.rm = TRUE),
    K_q025 = quantile_safe(group$K, 0.025), K_q975 = quantile_safe(group$K, 0.975),
    lambda_median = median(group$lambda, na.rm = TRUE),
    lambda_q025 = quantile_safe(group$lambda, 0.025),
    lambda_q975 = quantile_safe(group$lambda, 0.975), stringsAsFactors = FALSE
  )
})
random_summary <- do.call(rbind, summary_rows)
write.csv(random_summary, file.path(out_dir, "phylogenetic_signal_random_tree_summary.csv"), row.names = FALSE, na = "")
deterministic <- successful[successful$scenario %in% c("S1", "S3"), , drop = FALSE]
write.csv(deterministic, file.path(out_dir, "phylogenetic_signal_deterministic_and_direct.csv"), row.names = FALSE, na = "")

report <- list(
  generated_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  n_input_taxa = nrow(species), n_direct_backbone_taxa = length(direct_keys),
  n_endpoints = length(endpoints), n_random_trees = length(scenario2),
  min_species = min_species, K_permutation_replicates_for_deterministic_fits = k_nsim,
  n_total_fits = nrow(results), n_successful_fits = sum(results$status == "success"),
  n_failed_fits = sum(results$status == "failed"),
  n_insufficient_data_fits = sum(results$status == "insufficient_data"),
  interpretation = paste(
    "K and lambda are historical-signal sensitivities across a dated vascular-plant",
    "backbone and alternative grafting hypotheses. They do not resolve reticulate",
    "evolution, chloroplast capture or allopolyploid parentage."
  )
)
write_json(report, file.path(out_dir, "phylogenetic_signal_report.json"), pretty = TRUE, auto_unbox = TRUE)
cat(toJSON(report, pretty = TRUE, auto_unbox = TRUE), "\n")
