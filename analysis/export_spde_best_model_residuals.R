#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(INLA)
  library(Matrix)
  library(sf)
})

# Diagnostic-only refit of the frozen best-WAIC per-trait SPDE specifications.
# The accepted coefficient tables are not replaced. This script exists because
# the original compact workflow discarded fitted model objects and did not
# archive observation-level fitted values or residuals.

traits <- c(
  "orientation_angle_degrees_median",
  "corolla_lab_lightness_median",
  "corolla_lab_chroma_median",
  "corolla_hue_sin_median",
  "corolla_hue_cos_median",
  "shape_aspect_ratio_median",
  "shape_circularity_median",
  "shape_solidity_median",
  "shape_width_cv_median"
)

argv <- commandArgs(trailingOnly = TRUE)
getarg <- function(key, default = NULL) {
  i <- match(paste0("--", key), argv)
  if (is.na(i)) default else argv[[i + 1L]]
}

environment_file <- getarg("environment")
model_summary_file <- getarg("model-summary")
predictor_selection_file <- getarg("predictor-selection")
run_metadata_file <- getarg("run-metadata")
out_dir <- getarg("out-dir")
min_species_n <- as.integer(getarg("min-species-n", "5"))

if (any(vapply(
  list(environment_file, model_summary_file, predictor_selection_file, run_metadata_file, out_dir),
  is.null,
  logical(1)
))) {
  stop(
    "Required: --environment --model-summary --predictor-selection ",
    "--run-metadata --out-dir"
  )
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

df <- read.csv(environment_file, check.names = FALSE)
models <- read.csv(model_summary_file, check.names = FALSE)
selection <- read.csv(predictor_selection_file, check.names = FALSE)
run_metadata <- read.csv(run_metadata_file, check.names = FALSE)

required_observation <- c("obs_id", "taxon_name", "latitude", "longitude", traits)
missing_observation <- setdiff(required_observation, names(df))
if (length(missing_observation)) {
  stop("Environment table lacks required columns: ", paste(missing_observation, collapse = ", "))
}
required_models <- c("trait", "model_group", "status", "n_observations", "n_species", "waic")
missing_models <- setdiff(required_models, names(models))
if (length(missing_models)) {
  stop("Model summary lacks required columns: ", paste(missing_models, collapse = ", "))
}
required_selection <- c("model_group", "predictor", "status")
missing_selection <- setdiff(required_selection, names(selection))
if (length(missing_selection)) {
  stop("Predictor selection lacks required columns: ", paste(missing_selection, collapse = ", "))
}
required_metadata <- c("max_edge_km", "cutoff_km", "offset_km")
missing_metadata <- setdiff(required_metadata, names(run_metadata))
if (length(missing_metadata)) {
  stop("SPDE run metadata lacks required columns: ", paste(missing_metadata, collapse = ", "))
}
if (anyDuplicated(df$obs_id)) stop("Environment table must be unique by obs_id")

numeric_candidates <- unique(c(
  "latitude", "longitude", traits,
  selection$predictor[selection$status == "used"]
))
for (v in intersect(numeric_candidates, names(df))) {
  df[[v]] <- suppressWarnings(as.numeric(df[[v]]))
}

df <- df[is.finite(df$latitude) & is.finite(df$longitude), , drop = FALSE]
keep_species <- names(which(table(df$taxon_name) >= min_species_n))
df <- df[df$taxon_name %in% keep_species, , drop = FALSE]

# Reconstruct the exact within-species centered/scaled predictor matrix used by
# the frozen grouped-SPDE workflow. Predictor membership is read from the
# archived selection table rather than recomputed.
used_predictors <- unique(selection$predictor[selection$status == "used"])
missing_predictors <- setdiff(used_predictors, names(df))
if (length(missing_predictors)) {
  stop("Selected predictors absent from environment table: ", paste(missing_predictors, collapse = ", "))
}
Xall <- matrix(
  NA_real_,
  nrow(df),
  length(used_predictors),
  dimnames = list(NULL, used_predictors)
)
for (j in seq_along(used_predictors)) {
  p <- used_predictors[[j]]
  v <- df[[p]]
  centered <- v - ave(v, df$taxon_name, FUN = function(z) mean(z, na.rm = TRUE))
  s <- sd(centered, na.rm = TRUE)
  if (is.finite(s) && s > 0) Xall[, j] <- centered / s
}
if (any(vapply(seq_len(ncol(Xall)), function(j) all(!is.finite(Xall[, j])), logical(1)))) {
  stop("At least one archived selected predictor cannot be reconstructed")
}

# The frozen workflow used one global Equal Earth mesh shared by all traits.
max_edge <- as.numeric(run_metadata$max_edge_km[[1]])
cutoff <- as.numeric(run_metadata$cutoff_km[[1]])
offset <- as.numeric(run_metadata$offset_km[[1]])
pts <- st_as_sf(df, coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)
loc <- st_coordinates(st_transform(
  pts,
  "+proj=eqearth +datum=WGS84 +units=m +no_defs"
)) / 1000
mesh <- inla.mesh.2d(
  loc = loc,
  max.edge = c(max_edge, max_edge * 2),
  cutoff = cutoff,
  offset = c(offset, offset * 2)
)
spde <- inla.spde2.pcmatern(
  mesh,
  alpha = 2,
  prior.range = c(500, 0.5),
  prior.sigma = c(1, 0.01)
)
A <- inla.spde.make.A(mesh, loc = loc)
spatial_index <- inla.spde.make.index("spatial", mesh$n)

ok_models <- models[models$status == "ok" & is.finite(models$waic), , drop = FALSE]
if (!all(traits %in% ok_models$trait)) stop("Not every primary trait has a frozen successful model")
best_rows <- do.call(rbind, lapply(traits, function(trait) {
  z <- ok_models[ok_models$trait == trait, , drop = FALSE]
  z[which.min(z$waic), , drop = FALSE]
}))
rownames(best_rows) <- NULL

# Original workflow forced the same complete-case cohort across all four model
# groups for a trait by requiring every predictor used by any group.
union_preds <- unique(selection$predictor[selection$status == "used"])

prediction_tables <- list()
refit_rows <- list()

for (i in seq_len(nrow(best_rows))) {
  frozen <- best_rows[i, , drop = FALSE]
  trait <- as.character(frozen$trait[[1]])
  group <- as.character(frozen$model_group[[1]])
  preds <- selection$predictor[
    selection$model_group == group & selection$status == "used"
  ]
  preds <- unique(as.character(preds))
  if (!length(preds)) stop("No archived predictors for best model group: ", group)

  base_work <- is.finite(df[[trait]]) & complete.cases(Xall[, union_preds, drop = FALSE])
  counts <- table(df$taxon_name[base_work])
  ok_species <- names(counts[counts >= min_species_n])
  base_work <- base_work & df$taxon_name %in% ok_species

  n_obs <- sum(base_work)
  n_species <- length(unique(df$taxon_name[base_work]))
  if (n_obs != as.integer(frozen$n_observations[[1]])) {
    stop(
      sprintf(
        "%s observation count differs from frozen model: %d versus %d",
        trait, n_obs, as.integer(frozen$n_observations[[1]])
      )
    )
  }
  if (n_species != as.integer(frozen$n_species[[1]])) {
    stop(
      sprintf(
        "%s species count differs from frozen model: %d versus %d",
        trait, n_species, as.integer(frozen$n_species[[1]])
      )
    )
  }

  y_raw <- df[[trait]][base_work]
  y_mean <- mean(y_raw)
  y_sd <- sd(y_raw)
  y <- (y_raw - y_mean) / y_sd
  species <- as.integer(factor(df$taxon_name[base_work]))
  Awork <- A[base_work, , drop = FALSE]
  X <- as.data.frame(Xall[base_work, preds, drop = FALSE])
  names(X) <- make.names(names(X), unique = TRUE)

  stk <- inla.stack(
    data = list(y = y),
    A = list(1, Awork),
    effects = list(
      c(list(Intercept = rep(1, length(y)), species = species), X),
      spatial_index
    ),
    tag = "est"
  )
  form <- as.formula(paste(
    "y ~ 0 + Intercept +",
    paste(names(X), collapse = " + "),
    "+ f(species, model='iid') + f(spatial, model=spde)"
  ))

  fit <- inla(
    form,
    data = inla.stack.data(stk),
    family = "gaussian",
    control.predictor = list(A = inla.stack.A(stk), compute = TRUE),
    control.compute = list(waic = TRUE, dic = TRUE, cpo = TRUE, config = FALSE),
    verbose = FALSE
  )
  stack_index <- inla.stack.index(stk, "est")$data
  fitted <- fit$summary.fitted.values[stack_index, "mean"]
  if (length(fitted) != length(y)) stop("Unexpected fitted-value length for ", trait)

  prediction_tables[[trait]] <- data.frame(
    obs_id = as.character(df$obs_id[base_work]),
    taxon_name = as.character(df$taxon_name[base_work]),
    endpoint = trait,
    model_group = group,
    observed = as.numeric(y),
    observed_raw = as.numeric(y_raw),
    fitted = as.numeric(fitted),
    residual = as.numeric(y - fitted),
    latitude = as.numeric(df$latitude[base_work]),
    longitude = as.numeric(df$longitude[base_work]),
    stringsAsFactors = FALSE
  )

  refit_rows[[trait]] <- data.frame(
    trait = trait,
    model_group = group,
    n_observations = n_obs,
    n_species = n_species,
    n_predictors = length(preds),
    n_mesh_vertices = mesh$n,
    frozen_waic = as.numeric(frozen$waic[[1]]),
    diagnostic_refit_waic = as.numeric(fit$waic$waic),
    waic_difference = as.numeric(fit$waic$waic - frozen$waic[[1]]),
    response_mean_raw = y_mean,
    response_sd_raw = y_sd,
    inla_version = as.character(packageVersion("INLA")),
    stringsAsFactors = FALSE
  )

  rm(fit, stk)
  gc()
}

predictions <- do.call(rbind, prediction_tables)
rownames(predictions) <- NULL
if (anyDuplicated(predictions[c("obs_id", "endpoint")])) {
  stop("Diagnostic predictions are not unique by obs_id x endpoint")
}
write.csv(
  predictions,
  file.path(out_dir, "spde_best_model_observation_predictions.csv"),
  row.names = FALSE
)
write.csv(
  do.call(rbind, refit_rows),
  file.path(out_dir, "spde_best_model_refit_summary.csv"),
  row.names = FALSE
)
write.csv(
  best_rows,
  file.path(out_dir, "frozen_best_waic_model_registry.csv"),
  row.names = FALSE
)

cat(sprintf(
  "Exported %d endpoint-observation residual rows across %d traits.\n",
  nrow(predictions),
  length(unique(predictions$endpoint))
))
