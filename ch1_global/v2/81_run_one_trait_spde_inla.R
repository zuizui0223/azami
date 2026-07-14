#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(INLA)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) return(default)
  if (i == length(args)) stop(sprintf("Missing value for %s", flag))
  args[[i + 1L]]
}

input <- get_arg("--environment")
trait <- get_arg("--trait")
out_dir <- get_arg("--out-dir")
mesh_max_edge_km <- as.numeric(get_arg("--mesh-max-edge-km", "750"))
mesh_cutoff_km <- as.numeric(get_arg("--mesh-cutoff-km", "50"))
mesh_offset_km <- as.numeric(get_arg("--mesh-offset-km", "1000"))
min_species_n <- as.integer(get_arg("--min-species-n", "5"))

if (is.null(input) || is.null(trait) || is.null(out_dir)) {
  stop("Required: --environment FILE --trait NAME --out-dir DIR")
}
if (!file.exists(input)) stop(sprintf("Input not found: %s", input))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

climate <- c("chelsa_bio01", "chelsa_bio04", "chelsa_bio12", "chelsa_bio15")
topography <- c("topo_elevation", "topo_slope", "topo_roughness")
soil <- paste0(
  "soil_",
  c("bdod", "cec", "cfvo", "clay", "sand", "silt", "nitrogen", "phh2o", "soc", "ocd"),
  "_0_30cm"
)
groups <- list(
  climate = climate,
  climate_topography = c(climate, topography),
  climate_soil = c(climate, soil),
  climate_topography_soil = c(climate, topography, soil)
)

required_base <- c("obs_id", "taxon_name", "latitude", "longitude", trait)
df <- read.csv(input, stringsAsFactors = FALSE, check.names = FALSE)
missing_base <- setdiff(required_base, names(df))
if (length(missing_base)) stop(sprintf("Missing columns: %s", paste(missing_base, collapse = ", ")))

for (nm in unique(c(trait, climate, topography, soil, "latitude", "longitude"))) {
  if (nm %in% names(df)) df[[nm]] <- suppressWarnings(as.numeric(df[[nm]]))
}

# Keep only taxa with enough observations and genuine within-species response variation.
base_ok <- is.finite(df[[trait]]) & is.finite(df$latitude) & is.finite(df$longitude) & !is.na(df$taxon_name)
df <- df[base_ok, , drop = FALSE]
counts <- table(df$taxon_name)
keep_taxa <- names(counts[counts >= min_species_n])
df <- df[df$taxon_name %in% keep_taxa, , drop = FALSE]
var_by_species <- tapply(df[[trait]], df$taxon_name, function(x) length(unique(x[is.finite(x)])))
keep_taxa <- names(var_by_species[var_by_species >= 2])
df <- df[df$taxon_name %in% keep_taxa, , drop = FALSE]
if (nrow(df) < 100 || length(unique(df$taxon_name)) < 3) {
  stop(sprintf("Insufficient data after filtering: n=%d species=%d", nrow(df), length(unique(df$taxon_name))))
}

# Equal Earth projection (EPSG:8857) implemented directly to avoid a heavy sf dependency.
# Coordinates are returned in kilometres. This gives one globally continuous planar mesh;
# it is an approximation to geodesic distance and is recorded explicitly in the output.
equal_earth_xy_km <- function(lon_deg, lat_deg) {
  lon <- lon_deg * pi / 180
  lat <- lat_deg * pi / 180
  A1 <- 1.340264
  A2 <- -0.081106
  A3 <- 0.000893
  A4 <- 0.003796
  M <- sqrt(3) / 2
  theta <- asin(M * sin(lat))
  theta2 <- theta * theta
  theta6 <- theta2 * theta2 * theta2
  denom <- 3 * (9 * A4 * theta6 + 7 * A3 * theta2 * theta2 + 3 * A2 * theta2 + A1)
  x <- 2 * sqrt(3) * lon * cos(theta) / denom
  y <- theta * (A4 * theta6 + A3 * theta2 * theta2 + A2 * theta2 + A1)
  earth_radius_km <- 6371.0088
  cbind(x = x * earth_radius_km, y = y * earth_radius_km)
}

loc_all <- equal_earth_xy_km(df$longitude, df$latitude)

# Build one mesh per trait and reuse it across predictor-group fits. The mesh is based only
# on observation locations, so model-group comparisons share exactly the same spatial basis.
mesh <- inla.mesh.2d(
  loc = loc_all,
  max.edge = c(mesh_max_edge_km, mesh_max_edge_km * 2),
  cutoff = mesh_cutoff_km,
  offset = c(mesh_offset_km, mesh_offset_km * 2)
)

# Penalised-complexity priors: P(range < 500 km)=0.05 and P(sigma > 1)=0.05 after the
# response has been standardised within each model fit.
spde <- inla.spde2.pcmatern(
  mesh = mesh,
  alpha = 2,
  prior.range = c(500, 0.05),
  prior.sigma = c(1, 0.05)
)

fit_one_group <- function(group_name, predictors) {
  cols <- unique(c(required_base, predictors))
  work <- df[, cols, drop = FALSE]
  complete <- is.finite(work[[trait]]) & is.finite(work$latitude) & is.finite(work$longitude)
  for (p in predictors) complete <- complete & is.finite(work[[p]])
  work <- work[complete, , drop = FALSE]
  if (nrow(work) < 100 || length(unique(work$taxon_name)) < 3) {
    return(list(status = "insufficient", summary = data.frame(
      trait = trait, model_group = group_name, status = "insufficient",
      n_observations = nrow(work), n_species = length(unique(work$taxon_name))
    )))
  }

  # Mundlak-style within-species decomposition for the environmental predictors.
  # Only within-species deviations enter the fixed effects; the species random intercept
  # absorbs persistent among-species differences in trait means.
  for (p in predictors) {
    species_mean <- ave(work[[p]], work$taxon_name, FUN = mean)
    xw <- work[[p]] - species_mean
    s <- sd(xw)
    work[[paste0("xw__", p)]] <- if (is.finite(s) && s > 0) xw / s else 0
  }

  y_mean <- mean(work[[trait]])
  y_sd <- sd(work[[trait]])
  if (!is.finite(y_sd) || y_sd <= 0) {
    return(list(status = "insufficient", summary = data.frame(
      trait = trait, model_group = group_name, status = "zero_response_variance",
      n_observations = nrow(work), n_species = length(unique(work$taxon_name))
    )))
  }
  work$y_std <- (work[[trait]] - y_mean) / y_sd
  work$species_id <- as.integer(factor(work$taxon_name))
  species_n <- ave(work$obs_id, work$taxon_name, FUN = length)
  weights <- 1 / species_n
  weights <- weights / mean(weights)

  loc <- equal_earth_xy_km(work$longitude, work$latitude)
  A <- inla.spde.make.A(mesh = mesh, loc = loc)
  spatial_index <- inla.spde.make.index("spatial", spde$n.spde)
  xcols <- paste0("xw__", predictors)

  effects_fixed <- c(
    list(intercept = rep(1, nrow(work)), species_id = work$species_id),
    as.list(work[xcols])
  )
  stk <- inla.stack(
    data = list(y = work$y_std),
    A = list(1, A),
    effects = list(effects_fixed, spatial_index),
    tag = "est"
  )

  rhs <- c("-1", "intercept", xcols, "f(species_id, model='iid')", "f(spatial, model=spde)")
  formula <- as.formula(paste("y ~", paste(rhs, collapse = " + ")))

  fit <- inla(
    formula,
    data = inla.stack.data(stk, spde = spde),
    family = "gaussian",
    weights = weights,
    control.predictor = list(A = inla.stack.A(stk), compute = TRUE),
    control.compute = list(waic = TRUE, dic = TRUE, cpo = TRUE, config = FALSE),
    control.inla = list(strategy = "adaptive"),
    verbose = FALSE
  )

  fixed <- fit$summary.fixed
  fixed$term <- rownames(fixed)
  rownames(fixed) <- NULL
  fixed$trait <- trait
  fixed$model_group <- group_name
  fixed$estimate_original_units_per_1sd_within_predictor <- fixed$mean * y_sd
  fixed$lower_0.025_original_units <- fixed$`0.025quant` * y_sd
  fixed$upper_0.975_original_units <- fixed$`0.975quant` * y_sd

  hyper <- fit$summary.hyperpar
  hyper$parameter <- rownames(hyper)
  rownames(hyper) <- NULL
  hyper$trait <- trait
  hyper$model_group <- group_name

  cpo_fail <- if (!is.null(fit$cpo$failure)) sum(fit$cpo$failure > 0, na.rm = TRUE) else NA_integer_
  summary <- data.frame(
    trait = trait,
    model_group = group_name,
    status = "ok",
    n_observations = nrow(work),
    n_species = length(unique(work$taxon_name)),
    n_mesh_vertices = mesh$n,
    waic = fit$waic$waic,
    dic = fit$dic$dic,
    mean_log_cpo = if (!is.null(fit$cpo$cpo)) mean(log(fit$cpo$cpo), na.rm = TRUE) else NA_real_,
    cpo_failures = cpo_fail,
    response_mean = y_mean,
    response_sd = y_sd,
    stringsAsFactors = FALSE
  )

  list(status = "ok", summary = summary, fixed = fixed, hyper = hyper)
}

summaries <- list()
fixed_all <- list()
hyper_all <- list()
for (group_name in names(groups)) {
  message(sprintf("Fitting %s / %s", trait, group_name))
  result <- fit_one_group(group_name, groups[[group_name]])
  summaries[[group_name]] <- result$summary
  if (identical(result$status, "ok")) {
    fixed_all[[group_name]] <- result$fixed
    hyper_all[[group_name]] <- result$hyper
  }
}

summary_df <- do.call(rbind, summaries)
write.csv(summary_df, file.path(out_dir, "spde_inla_model_summary.csv"), row.names = FALSE)
if (length(fixed_all)) {
  write.csv(do.call(rbind, fixed_all), file.path(out_dir, "spde_inla_fixed_effects.csv"), row.names = FALSE)
}
if (length(hyper_all)) {
  write.csv(do.call(rbind, hyper_all), file.path(out_dir, "spde_inla_hyperparameters.csv"), row.names = FALSE)
}

mesh_info <- data.frame(
  trait = trait,
  projection = "Equal Earth EPSG:8857 analytic approximation",
  coordinate_units = "km",
  n_mesh_vertices = mesh$n,
  mesh_max_edge_km = mesh_max_edge_km,
  mesh_cutoff_km = mesh_cutoff_km,
  mesh_offset_km = mesh_offset_km,
  spde_prior_range_km = 500,
  spde_prior_range_probability = 0.05,
  spde_prior_sigma = 1,
  spde_prior_sigma_probability = 0.05,
  stringsAsFactors = FALSE
)
write.csv(mesh_info, file.path(out_dir, "spde_inla_mesh_metadata.csv"), row.names = FALSE)

cat(sprintf("Completed SPDE-INLA analysis for %s\n", trait))
print(summary_df)
