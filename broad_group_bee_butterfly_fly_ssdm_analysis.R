############################################################
## Broad pollinator group SSDM + improved grid analysis
##
## Groups:
##   bee        = all bee genera in current SDM set
##   butterfly  = all butterfly genera in current SDM set
##   fly        = hoverfly / fly genera in current SDM set
##
## Optional:
##   hawkmoth is saved separately but excluded from main 3-group model.
##
## Aim:
##   Run the same improved analysis as before, but using simple broad groups:
##     bee, butterfly, fly
##
## Main response:
##   cbind(n_nodding_species, n_upward_species)
##
## Model families:
##   1. binomial GLM
##   2. beta-binomial GLM via glmmTMB
##   3. spatial GAM via mgcv
##   4. spatial block CV
##
## Inputs:
##   chelsa_pollinator_enmeval_rebuild_no_swe/
##     species_enmeval_run_summary.csv
##     cirsium_candidateEnvPCA_pollinator_metric_reduced/
##       01_cirsium_occurrences_with_candidateEnv_and_guildSDM.csv
##       03_candidateEnv_PCA_model.rds
##
## Output:
##   chelsa_pollinator_enmeval_rebuild_no_swe/
##     broad_group_bee_butterfly_fly_ssdm_analysis/
############################################################

## =========================================================
## 0. Packages
## =========================================================

pkgs <- c(
  "terra", "dplyr", "readr", "stringr", "tibble", "tidyr",
  "ggplot2", "pROC", "patchwork", "forcats",
  "mgcv", "glmmTMB"
)

for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, dependencies = TRUE)
  }
}

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
  library(tidyr)
  library(ggplot2)
  library(pROC)
  library(patchwork)
  library(forcats)
  library(mgcv)
  library(glmmTMB)
})

set.seed(42)

## =========================================================
## 1. Paths and settings
## =========================================================

ROOT_DIR <- "chelsa_pollinator_enmeval_rebuild_no_swe"

SPECIES_SUMMARY_FILE <- file.path(
  ROOT_DIR,
  "species_enmeval_run_summary.csv"
)

BASE_ANALYSIS_DIR <- file.path(
  ROOT_DIR,
  "cirsium_candidateEnvPCA_pollinator_metric_reduced"
)

CIRSIUM_OCC_FILE <- file.path(
  BASE_ANALYSIS_DIR,
  "01_cirsium_occurrences_with_candidateEnv_and_guildSDM.csv"
)

PCA_MODEL_FILE <- file.path(
  BASE_ANALYSIS_DIR,
  "03_candidateEnv_PCA_model.rds"
)

OUT_DIR <- file.path(
  ROOT_DIR,
  "broad_group_bee_butterfly_fly_ssdm_analysis"
)

STACK_DIR <- file.path(OUT_DIR, "broad_group_stacks")
FIG_DIR <- file.path(OUT_DIR, "figures")
TMP_DIR <- file.path(OUT_DIR, "terra_tmp")

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(STACK_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TMP_DIR, showWarnings = FALSE, recursive = TRUE)

terra::terraOptions(memfrac = 0.45, tempdir = TMP_DIR)

GRID_SIZE_DEG <- 0.5
MIN_TOTAL_SPECIES_PER_GRID <- 2

FORCE_REBUILD_STACKS <- FALSE

## Spatial block CV
BLOCK_SIZE_DEG <- 2.0
N_CV_REPEAT <- 20

## Include hawkmoth in main model?
## Main requested analysis is bee / butterfly / fly, so hawkmoth = FALSE.
INCLUDE_HAWKMOTH_IN_MAIN <- FALSE

## =========================================================
## 2. Broad pollinator classification
## =========================================================

BEE_GENERA <- c(
  "Bombus", "Xylocopa", "Eucera",
  "Apis", "Megachile", "Ceratina", "Lasioglossum",
  "Andrena", "Colletes", "Nomada", "Osmia",
  "Halictus", "Hylaeus", "Sphecodes", "Seladonia",
  "Anthophora"
)

BUTTERFLY_GENERA <- c(
  "Pieris", "Aglais", "Aporia", "Araschnia", "Argynnis",
  "Fabriciana", "Lycaena", "Nymphalis", "Ochlodes",
  "Papilio", "Polygonia", "Satyrium", "Speyeria"
)

FLY_GENERA <- c(
  "Syrphus", "Eristalis", "Melanostoma", "Volucella",
  "Dasysyrphus", "Didea", "Episyrphus", "Eupeodes",
  "Ferdinandea", "Leucozona", "Meliscaeva", "Merodon",
  "Platycheirus", "Rhingia", "Sphaerophoria"
)

HAWKMOTH_GENERA <- c("Macroglossum")

BROAD_GROUPS <- list(
  bee = BEE_GENERA,
  butterfly = BUTTERFLY_GENERA,
  fly = FLY_GENERA,
  hawkmoth = HAWKMOTH_GENERA
)

classify_broad_group <- function(genus) {
  dplyr::case_when(
    genus %in% BEE_GENERA ~ "bee",
    genus %in% BUTTERFLY_GENERA ~ "butterfly",
    genus %in% FLY_GENERA ~ "fly",
    genus %in% HAWKMOTH_GENERA ~ "hawkmoth",
    TRUE ~ "unknown"
  )
}

## =========================================================
## 3. Helper functions
## =========================================================

stop_if_missing <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)
}

safe_name <- function(x) {
  x |>
    as.character() |>
    stringr::str_replace_all("[^A-Za-z0-9_\\-]+", "_") |>
    stringr::str_replace_all("_+", "_") |>
    stringr::str_replace_all("^_|_$", "")
}

extract_genus <- function(species_name) {
  genus <- stringr::str_extract(as.character(species_name), "^[A-Za-z]+")
  ifelse(is.na(genus) | genus == "", "Unknown", genus)
}

to_numeric_coord <- function(x) {
  if (is.list(x)) x <- unlist(x)
  x <- as.character(x)
  x <- stringr::str_replace_all(x, ",", ".")
  suppressWarnings(as.numeric(x))
}

safe_extract_df <- function(r, xy_or_vect) {
  out <- terra::extract(r, xy_or_vect)
  out <- as.data.frame(out)

  if (ncol(out) > terra::nlyr(r)) {
    id_like <- names(out)[1]
    if (tolower(id_like) %in% c("id", "id.1") || id_like == "ID") {
      out <- out[, -1, drop = FALSE]
    }
  }

  out
}

safe_scale <- function(x) {
  x <- as.numeric(x)
  s <- stats::sd(x, na.rm = TRUE)
  m <- mean(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(0, length(x)))
  (x - m) / s
}

safe_auc <- function(y, p) {
  tryCatch(as.numeric(pROC::auc(y, p, quiet = TRUE)), error = function(e) NA_real_)
}

logloss_binomial <- function(y, n, p, eps = 1e-8) {
  p <- pmin(pmax(p, eps), 1 - eps)
  -mean(y * log(p) + (n - y) * log(1 - p), na.rm = TRUE)
}

brier_prop <- function(y_prop, p) {
  mean((y_prop - p)^2, na.rm = TRUE)
}

make_formula <- function(rhs) {
  stats::as.formula(paste("cbind(n_nodding_species, n_upward_species) ~", rhs))
}

get_rhs <- function(vars) {
  vars <- vars[!is.na(vars) & vars != ""]
  if (length(vars) == 0) return("1")
  paste(vars, collapse = " + ")
}

empty_summary <- function(model, method, error = NA_character_) {
  tibble::tibble(
    model = model,
    method = method,
    AIC = NA_real_,
    BIC = NA_real_,
    n_grid = NA_integer_,
    logLik = NA_real_,
    df = NA_real_,
    AUC_presence = NA_real_,
    logloss = NA_real_,
    brier = NA_real_,
    deviance_explained = NA_real_,
    error = error
  )
}

empty_coef <- function(model, method, error = NA_character_) {
  tibble::tibble(
    model = model,
    method = method,
    term = NA_character_,
    estimate = NA_real_,
    std.error = NA_real_,
    statistic = NA_real_,
    p.value = NA_real_,
    error = error
  )
}

empty_smooth <- function(model, method, error = NA_character_) {
  tibble::tibble(
    model = model,
    method = method,
    term = NA_character_,
    edf = NA_real_,
    statistic = NA_real_,
    p.value = NA_real_,
    error = error
  )
}

safe_arrange_aic <- function(df) {
  if (!"AIC" %in% names(df)) df$AIC <- NA_real_

  min_aic <- suppressWarnings(min(df$AIC, na.rm = TRUE))
  if (!is.finite(min_aic)) {
    df$delta_AIC <- NA_real_
  } else {
    df$delta_AIC <- df$AIC - min_aic
  }

  df %>% dplyr::arrange(.data$AIC)
}

coef_from_glm <- function(m, model_name, method_name) {
  sm <- summary(m)$coefficients
  as.data.frame(sm) %>%
    tibble::rownames_to_column("term") %>%
    tibble::as_tibble() %>%
    dplyr::rename(
      estimate = Estimate,
      std.error = `Std. Error`,
      statistic = `z value`,
      p.value = `Pr(>|z|)`
    ) %>%
    dplyr::mutate(model = model_name, method = method_name, .before = 1)
}

coef_from_glmmTMB <- function(m, model_name, method_name) {
  sm <- summary(m)$coefficients$cond
  as.data.frame(sm) %>%
    tibble::rownames_to_column("term") %>%
    tibble::as_tibble() %>%
    dplyr::rename(
      estimate = Estimate,
      std.error = `Std. Error`,
      statistic = `z value`,
      p.value = `Pr(>|z|)`
    ) %>%
    dplyr::mutate(model = model_name, method = method_name, .before = 1)
}

coef_from_gam_param <- function(m, model_name, method_name) {
  sm <- summary(m)$p.table
  as.data.frame(sm) %>%
    tibble::rownames_to_column("term") %>%
    tibble::as_tibble() %>%
    dplyr::rename(
      estimate = Estimate,
      std.error = `Std. Error`,
      statistic = `z value`,
      p.value = `Pr(>|z|)`
    ) %>%
    dplyr::mutate(model = model_name, method = method_name, .before = 1)
}

smooth_from_gam <- function(m, model_name, method_name) {
  sm <- summary(m)$s.table
  if (is.null(sm) || nrow(sm) == 0) {
    return(empty_smooth(model_name, method_name))
  }

  as.data.frame(sm) %>%
    tibble::rownames_to_column("term") %>%
    tibble::as_tibble() %>%
    dplyr::rename(
      edf = edf,
      statistic = `Chi.sq`,
      p.value = `p-value`
    ) %>%
    dplyr::mutate(model = model_name, method = method_name, .before = 1)
}

fit_glm_binom <- function(formula, data, name) {
  method <- "binomial_GLM"

  tryCatch({
    m <- stats::glm(formula, data = data, family = stats::binomial)
    pred <- stats::predict(m, type = "response")
    y_presence <- ifelse(data$n_nodding_species > 0, 1, 0)

    list(
      model = m,
      summary = tibble::tibble(
        model = name,
        method = method,
        AIC = stats::AIC(m),
        BIC = stats::BIC(m),
        n_grid = stats::nobs(m),
        logLik = as.numeric(stats::logLik(m)),
        df = attr(stats::logLik(m), "df"),
        AUC_presence = safe_auc(y_presence, pred),
        logloss = logloss_binomial(data$n_nodding_species, data$n_total_species, pred),
        brier = brier_prop(data$prop_nodding_species, pred),
        deviance_explained = NA_real_,
        error = NA_character_
      ),
      coef = coef_from_glm(m, name, method)
    )
  }, error = function(e) {
    list(
      model = NULL,
      summary = empty_summary(name, method, conditionMessage(e)),
      coef = empty_coef(name, method, conditionMessage(e))
    )
  })
}

fit_glmmTMB_betabinom <- function(formula, data, name) {
  method <- "beta_binomial_glmmTMB"

  tryCatch({
    m <- glmmTMB::glmmTMB(formula, data = data, family = glmmTMB::betabinomial(link = "logit"))
    pred <- stats::predict(m, type = "response")
    y_presence <- ifelse(data$n_nodding_species > 0, 1, 0)

    list(
      model = m,
      summary = tibble::tibble(
        model = name,
        method = method,
        AIC = stats::AIC(m),
        BIC = stats::BIC(m),
        n_grid = nrow(data),
        logLik = as.numeric(stats::logLik(m)),
        df = attr(stats::logLik(m), "df"),
        AUC_presence = safe_auc(y_presence, pred),
        logloss = logloss_binomial(data$n_nodding_species, data$n_total_species, pred),
        brier = brier_prop(data$prop_nodding_species, pred),
        deviance_explained = NA_real_,
        error = NA_character_
      ),
      coef = coef_from_glmmTMB(m, name, method)
    )
  }, error = function(e) {
    list(
      model = NULL,
      summary = empty_summary(name, method, conditionMessage(e)),
      coef = empty_coef(name, method, conditionMessage(e))
    )
  })
}

fit_gam_binom <- function(formula, data, name) {
  method <- "spatial_GAM_binomial"

  tryCatch({
    m <- mgcv::gam(formula, data = data, family = stats::binomial, method = "REML")
    pred <- stats::predict(m, type = "response")
    y_presence <- ifelse(data$n_nodding_species > 0, 1, 0)

    list(
      model = m,
      summary = tibble::tibble(
        model = name,
        method = method,
        AIC = stats::AIC(m),
        BIC = stats::BIC(m),
        n_grid = nrow(data),
        logLik = as.numeric(stats::logLik(m)),
        df = sum(m$edf),
        AUC_presence = safe_auc(y_presence, pred),
        logloss = logloss_binomial(data$n_nodding_species, data$n_total_species, pred),
        brier = brier_prop(data$prop_nodding_species, pred),
        deviance_explained = summary(m)$dev.expl,
        error = NA_character_
      ),
      coef = coef_from_gam_param(m, name, "spatial_GAM_binomial_parametric"),
      smooth = smooth_from_gam(m, name, "spatial_GAM_binomial_smooth")
    )
  }, error = function(e) {
    list(
      model = NULL,
      summary = empty_summary(name, method, conditionMessage(e)),
      coef = empty_coef(name, "spatial_GAM_binomial_parametric", conditionMessage(e)),
      smooth = empty_smooth(name, "spatial_GAM_binomial_smooth", conditionMessage(e))
    )
  })
}

make_broad_stack <- function(group_name, genera, sp_tbl) {
  group_safe <- safe_name(group_name)
  group_dir <- file.path(STACK_DIR, group_safe)
  dir.create(group_dir, showWarnings = FALSE, recursive = TRUE)

  sum_file <- file.path(group_dir, paste0("broad_", group_safe, "_sum_cloglog.tif"))
  mean_file <- file.path(group_dir, paste0("broad_", group_safe, "_mean_cloglog.tif"))
  rich_file <- file.path(group_dir, paste0("broad_", group_safe, "_binary_richness_10p.tif"))

  sub <- sp_tbl |> dplyr::filter(genus %in% genera)

  if (nrow(sub) == 0) {
    warning("No species SDM found for broad group: ", group_name)
    return(NULL)
  }

  if (!FORCE_REBUILD_STACKS && file.exists(sum_file) && file.exists(mean_file)) {
    return(tibble::tibble(
      broad_group = group_name,
      group_safe = group_safe,
      n_species_sdm = dplyr::n_distinct(sub$species),
      n_genera = dplyr::n_distinct(sub$genus),
      genera_available = paste(sort(unique(sub$genus)), collapse = "; "),
      species_list = paste(sort(unique(sub$species)), collapse = "; "),
      sum_cloglog_file = sum_file,
      mean_cloglog_file = mean_file,
      binary_richness_file = ifelse(file.exists(rich_file), rich_file, NA_character_)
    ))
  }

  message("Building broad SSDM: ", group_name)

  pred_files <- sub$pred_file[file.exists(sub$pred_file)]
  preds <- terra::rast(pred_files)
  names(preds) <- safe_name(tools::file_path_sans_ext(basename(pred_files)))

  sum_r <- terra::app(preds, fun = sum, na.rm = TRUE)
  mean_r <- terra::app(preds, fun = mean, na.rm = TRUE)

  names(sum_r) <- paste0("broad_", group_safe, "_sum_cloglog")
  names(mean_r) <- paste0("broad_", group_safe, "_mean_cloglog")

  terra::writeRaster(sum_r, sum_file, overwrite = TRUE)
  terra::writeRaster(mean_r, mean_file, overwrite = TRUE)

  rich_out <- NA_character_
  bin_files <- sub$bin_file[sub$bin_file_exists]

  if (length(bin_files) > 0) {
    bins <- terra::rast(bin_files)
    rich_r <- terra::app(bins, fun = sum, na.rm = TRUE)
    names(rich_r) <- paste0("broad_", group_safe, "_binary_richness_10p")
    terra::writeRaster(rich_r, rich_file, overwrite = TRUE)
    rich_out <- rich_file
  }

  gc(verbose = FALSE)

  tibble::tibble(
    broad_group = group_name,
    group_safe = group_safe,
    n_species_sdm = dplyr::n_distinct(sub$species),
    n_genera = dplyr::n_distinct(sub$genus),
    genera_available = paste(sort(unique(sub$genus)), collapse = "; "),
    species_list = paste(sort(unique(sub$species)), collapse = "; "),
    sum_cloglog_file = sum_file,
    mean_cloglog_file = mean_file,
    binary_richness_file = rich_out
  )
}

## =========================================================
## 4. Load pollinator species SDMs and build broad SSDMs
## =========================================================

stop_if_missing(SPECIES_SUMMARY_FILE)
stop_if_missing(CIRSIUM_OCC_FILE)
stop_if_missing(PCA_MODEL_FILE)

sp_sum <- readr::read_csv(SPECIES_SUMMARY_FILE, show_col_types = FALSE) |>
  dplyr::filter(status %in% c("success", "skipped_existing")) |>
  dplyr::filter(!is.na(pred_file), file.exists(pred_file)) |>
  dplyr::mutate(
    genus = extract_genus(species),
    broad_group = classify_broad_group(genus),
    bin_file_exists = !is.na(bin_file) & file.exists(bin_file)
  )

classification_summary <- sp_sum |>
  dplyr::group_by(broad_group, genus) |>
  dplyr::summarise(
    n_species_sdm = dplyr::n_distinct(species),
    species_list = paste(sort(unique(species)), collapse = "; "),
    .groups = "drop"
  ) |>
  dplyr::arrange(broad_group, genus)

readr::write_csv(classification_summary, file.path(OUT_DIR, "00_pollinator_broad_group_classification.csv"))

message("Pollinator broad-group classification:")
print(classification_summary, n = Inf, width = Inf)

unknowns <- classification_summary |> dplyr::filter(broad_group == "unknown")
if (nrow(unknowns) > 0) {
  warning("Some genera are unknown. Check 00_pollinator_broad_group_classification.csv")
}

stack_results <- list()

for (grp in names(BROAD_GROUPS)) {
  stack_results[[grp]] <- make_broad_stack(
    group_name = grp,
    genera = BROAD_GROUPS[[grp]],
    sp_tbl = sp_sum
  )
}

stack_summary <- dplyr::bind_rows(stack_results) |>
  dplyr::arrange(broad_group)

readr::write_csv(stack_summary, file.path(OUT_DIR, "01_broad_group_ssdm_stack_summary.csv"))

message("\nBroad-group SSDM stack summary:")
print(stack_summary, n = Inf, width = Inf)

## =========================================================
## 5. Extract broad SSDMs to Cirsium occurrence points
## =========================================================

occ0 <- readr::read_csv(CIRSIUM_OCC_FILE, show_col_types = FALSE) |>
  dplyr::mutate(
    decimalLongitude = to_numeric_coord(decimalLongitude),
    decimalLatitude = to_numeric_coord(decimalLatitude)
  ) |>
  dplyr::filter(
    is.finite(decimalLongitude),
    is.finite(decimalLatitude),
    head_orientation_binary %in% c("upward", "nodding"),
    !is.na(species_final),
    species_final != ""
  ) |>
  dplyr::mutate(
    y_nodding = ifelse(head_orientation_binary == "nodding", 1, 0)
  )

old_guild_cols <- names(occ0)[stringr::str_detect(names(occ0), "^guild_")]
old_hyp_cols <- names(occ0)[stringr::str_detect(names(occ0), "^hyp_")]
occ_base <- occ0 |> dplyr::select(-dplyr::all_of(c(old_guild_cols, old_hyp_cols)))

all_stack_files <- c(
  stack_summary$sum_cloglog_file,
  stack_summary$mean_cloglog_file,
  stack_summary$binary_richness_file
)
all_stack_files <- all_stack_files[!is.na(all_stack_files) & file.exists(all_stack_files)]

if (length(all_stack_files) == 0) stop("No broad SSDM rasters found.")

broad_r <- terra::rast(all_stack_files)
names(broad_r) <- safe_name(tools::file_path_sans_ext(basename(all_stack_files)))

xy <- cbind(occ_base$decimalLongitude, occ_base$decimalLatitude)
colnames(xy) <- c("lon", "lat")

broad_vals <- safe_extract_df(broad_r, xy)
names(broad_vals) <- names(broad_r)

occ_broad <- dplyr::bind_cols(occ_base, broad_vals, .name_repair = "unique")

readr::write_csv(occ_broad, file.path(OUT_DIR, "02_cirsium_occurrences_with_broad_group_SSDM.csv"))

## =========================================================
## 6. Grid table
## =========================================================

occ_grid <- occ_broad |>
  dplyr::mutate(
    grid_lon = floor(decimalLongitude / GRID_SIZE_DEG) * GRID_SIZE_DEG + GRID_SIZE_DEG / 2,
    grid_lat = floor(decimalLatitude / GRID_SIZE_DEG) * GRID_SIZE_DEG + GRID_SIZE_DEG / 2,
    grid_id = paste0(round(grid_lon, 3), "_", round(grid_lat, 3))
  )

grid_species <- occ_grid |>
  dplyr::distinct(
    grid_id, grid_lon, grid_lat,
    species_final, japanese_name, head_orientation_binary,
    .keep_all = TRUE
  )

grid_counts <- grid_species |>
  dplyr::group_by(grid_id, grid_lon, grid_lat) |>
  dplyr::summarise(
    n_total_species = dplyr::n_distinct(species_final),
    n_nodding_species = dplyr::n_distinct(species_final[head_orientation_binary == "nodding"]),
    n_upward_species = dplyr::n_distinct(species_final[head_orientation_binary == "upward"]),
    prop_nodding_species = n_nodding_species / n_total_species,
    logit_nodding_species = log((n_nodding_species + 0.5) / (n_upward_species + 0.5)),
    .groups = "drop"
  )

env_cols <- names(occ_grid)[stringr::str_detect(names(occ_grid), "^pcaenvCand_")]
broad_cols <- names(occ_grid)[stringr::str_detect(names(occ_grid), "^broad_")]

grid_predictors <- occ_grid |>
  dplyr::group_by(grid_id, grid_lon, grid_lat) |>
  dplyr::summarise(
    n_occ = dplyr::n(),
    dplyr::across(
      dplyr::all_of(c(env_cols, broad_cols)),
      ~ mean(.x, na.rm = TRUE),
      .names = "{.col}"
    ),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    dplyr::across(where(is.numeric), ~ ifelse(is.nan(.x), NA_real_, .x))
  )

grid_dat <- grid_counts |>
  dplyr::left_join(grid_predictors, by = c("grid_id", "grid_lon", "grid_lat"))

readr::write_csv(grid_dat, file.path(OUT_DIR, "03_grid_nodding_species_with_broad_group_SSDM_raw.csv"))

## =========================================================
## 7. Environmental PCA projection
## =========================================================

pca <- readRDS(PCA_MODEL_FILE)
pca_vars <- rownames(pca$rotation)

missing_pca_vars <- setdiff(pca_vars, names(grid_dat))
if (length(missing_pca_vars) > 0) {
  stop("Grid table missing PCA variables:\n", paste(missing_pca_vars, collapse = "\n"))
}

env_mat <- grid_dat |>
  dplyr::select(dplyr::all_of(pca_vars)) |>
  as.data.frame()

for (nm in names(env_mat)) {
  x <- env_mat[[nm]]
  x[!is.finite(x)] <- mean(x, na.rm = TRUE)
  env_mat[[nm]] <- x
}

pc_scores <- predict(pca, newdata = env_mat) |>
  as.data.frame() |>
  tibble::as_tibble()

pc_keep <- intersect(c("PC1", "PC2", "PC3"), names(pc_scores))
pc_scores <- pc_scores[, pc_keep, drop = FALSE]
names(pc_scores) <- paste0("grid_", names(pc_scores))

grid_dat <- dplyr::bind_cols(grid_dat, pc_scores)

for (nm in names(pc_scores)) {
  grid_dat[[paste0("z_", nm)]] <- safe_scale(grid_dat[[nm]])
}

for (nm in broad_cols) {
  if (stringr::str_detect(nm, "sum_cloglog|binary_richness")) {
    grid_dat[[nm]][is.na(grid_dat[[nm]])] <- 0
  }
  grid_dat[[paste0("z_", nm)]] <- safe_scale(grid_dat[[nm]])
}

analysis_grid <- grid_dat |>
  dplyr::filter(n_total_species >= MIN_TOTAL_SPECIES_PER_GRID) |>
  tidyr::drop_na(
    n_nodding_species,
    n_upward_species,
    n_total_species,
    z_grid_PC1,
    z_grid_PC2,
    z_grid_PC3
  ) |>
  dplyr::mutate(
    prop_nodding_species = n_nodding_species / n_total_species,
    y_presence = ifelse(n_nodding_species > 0, 1, 0),
    z_lon = safe_scale(grid_lon),
    z_lat = safe_scale(grid_lat)
  )

readr::write_csv(analysis_grid, file.path(OUT_DIR, "04_analysis_grid_data_broad_group_SSDM.csv"))

message("\nAnalysis grid summary:")
print(analysis_grid |> dplyr::summarise(
  n_grid = dplyr::n(),
  mean_total_species = mean(n_total_species),
  total_nodding_grid_counts = sum(n_nodding_species),
  total_upward_grid_counts = sum(n_upward_species),
  mean_prop = mean(prop_nodding_species)
))

## =========================================================
## 8. Define model variables
## =========================================================

ENV_VARS <- c("z_grid_PC1", "z_grid_PC2", "z_grid_PC3")
ENV_VARS <- ENV_VARS[ENV_VARS %in% names(analysis_grid)]

POLL_VARS_MAIN <- c(
  "z_broad_bee_sum_cloglog",
  "z_broad_butterfly_sum_cloglog",
  "z_broad_fly_sum_cloglog"
)
POLL_VARS_MAIN <- POLL_VARS_MAIN[POLL_VARS_MAIN %in% names(analysis_grid)]

POLL_VARS_WITH_HAWKMOTH <- c(POLL_VARS_MAIN, "z_broad_hawkmoth_sum_cloglog")
POLL_VARS_WITH_HAWKMOTH <- POLL_VARS_WITH_HAWKMOTH[POLL_VARS_WITH_HAWKMOTH %in% names(analysis_grid)]

if (!INCLUDE_HAWKMOTH_IN_MAIN) {
  POLL_VARS <- POLL_VARS_MAIN
} else {
  POLL_VARS <- POLL_VARS_WITH_HAWKMOTH
}

for (v in c(ENV_VARS, POLL_VARS, POLL_VARS_WITH_HAWKMOTH)) {
  analysis_grid[[v]][is.na(analysis_grid[[v]])] <- 0
}

predictor_summary <- tibble::tibble(
  class = c(rep("env", length(ENV_VARS)), rep("poll_main", length(POLL_VARS)), rep("poll_with_hawkmoth", length(POLL_VARS_WITH_HAWKMOTH))),
  variable = c(ENV_VARS, POLL_VARS, POLL_VARS_WITH_HAWKMOTH)
) |> dplyr::distinct()

readr::write_csv(predictor_summary, file.path(OUT_DIR, "05_predictors_used_broad_group.csv"))

message("Predictors used:")
print(predictor_summary, n = Inf)

## =========================================================
## 9. Model specifications
## =========================================================

rhs_env <- get_rhs(ENV_VARS)
rhs_poll <- get_rhs(POLL_VARS)
rhs_poll_hawk <- get_rhs(POLL_VARS_WITH_HAWKMOTH)

rhs_env_poll <- get_rhs(c(ENV_VARS, POLL_VARS))
rhs_env_poll_hawk <- get_rhs(c(ENV_VARS, POLL_VARS_WITH_HAWKMOTH))

rhs_space <- "z_lon + z_lat + I(z_lon^2) + I(z_lat^2) + z_lon:z_lat"

rhs_env_space <- paste(rhs_env, rhs_space, sep = " + ")
rhs_env_poll_space <- paste(rhs_env_poll, rhs_space, sep = " + ")
rhs_env_poll_hawk_space <- paste(rhs_env_poll_hawk, rhs_space, sep = " + ")

model_specs <- list(
  null = "1",
  space_poly = rhs_space,
  env = rhs_env,
  poll_broad = rhs_poll,
  poll_broad_with_hawkmoth = rhs_poll_hawk,
  env_plus_poll_broad = rhs_env_poll,
  env_plus_poll_broad_with_hawkmoth = rhs_env_poll_hawk,
  env_plus_space = rhs_env_space,
  env_plus_poll_broad_plus_space = rhs_env_poll_space,
  env_plus_poll_broad_with_hawkmoth_plus_space = rhs_env_poll_hawk_space
)

model_specs <- model_specs[!vapply(model_specs, function(x) is.na(x) || x == "", logical(1))]

gam_specs <- list(
  gam_space = "s(grid_lon, grid_lat, k = 40)",
  gam_env_plus_space = paste(rhs_env, "s(grid_lon, grid_lat, k = 40)", sep = " + "),
  gam_env_poll_broad_plus_space = paste(rhs_env_poll, "s(grid_lon, grid_lat, k = 40)", sep = " + "),
  gam_env_poll_broad_with_hawkmoth_plus_space = paste(rhs_env_poll_hawk, "s(grid_lon, grid_lat, k = 40)", sep = " + ")
)

## =========================================================
## 10. Fit GLM and beta-binomial
## =========================================================

glm_results <- list()
bb_results <- list()

for (nm in names(model_specs)) {
  f <- make_formula(model_specs[[nm]])
  message("Fitting GLM and beta-binomial: ", nm)
  glm_results[[nm]] <- fit_glm_binom(f, analysis_grid, nm)
  bb_results[[nm]] <- fit_glmmTMB_betabinom(f, analysis_grid, nm)
}

glm_summary <- dplyr::bind_rows(lapply(glm_results, `[[`, "summary")) %>% safe_arrange_aic()
bb_summary <- dplyr::bind_rows(lapply(bb_results, `[[`, "summary")) %>% safe_arrange_aic()

glm_coef <- dplyr::bind_rows(lapply(glm_results, `[[`, "coef"))
bb_coef <- dplyr::bind_rows(lapply(bb_results, `[[`, "coef"))

readr::write_csv(glm_summary, file.path(OUT_DIR, "06_binomial_GLM_model_comparison.csv"))
readr::write_csv(glm_coef, file.path(OUT_DIR, "06_binomial_GLM_coefficients.csv"))
readr::write_csv(bb_summary, file.path(OUT_DIR, "07_beta_binomial_model_comparison.csv"))
readr::write_csv(bb_coef, file.path(OUT_DIR, "07_beta_binomial_coefficients.csv"))

message("\nBinomial GLM model comparison:")
print(glm_summary, n = Inf, width = Inf)

message("\nBeta-binomial model comparison:")
print(bb_summary, n = Inf, width = Inf)

## =========================================================
## 11. Fit spatial GAM
## =========================================================

gam_results <- list()

for (nm in names(gam_specs)) {
  f <- make_formula(gam_specs[[nm]])
  message("Fitting spatial GAM: ", nm)
  gam_results[[nm]] <- fit_gam_binom(f, analysis_grid, nm)
}

gam_summary <- dplyr::bind_rows(lapply(gam_results, `[[`, "summary")) %>% safe_arrange_aic()
gam_coef <- dplyr::bind_rows(lapply(gam_results, `[[`, "coef"))
gam_smooth <- dplyr::bind_rows(lapply(gam_results, `[[`, "smooth"))

readr::write_csv(gam_summary, file.path(OUT_DIR, "08_spatial_GAM_model_comparison.csv"))
readr::write_csv(gam_coef, file.path(OUT_DIR, "08_spatial_GAM_parametric_coefficients.csv"))
readr::write_csv(gam_smooth, file.path(OUT_DIR, "08_spatial_GAM_smooth_terms.csv"))

message("\nSpatial GAM model comparison:")
print(gam_summary, n = Inf, width = Inf)

message("\nSpatial GAM parametric coefficients:")
print(gam_coef %>% dplyr::filter(term != "(Intercept)") %>% dplyr::arrange(p.value), n = Inf, width = Inf)

message("\nSpatial GAM smooth terms:")
print(gam_smooth, n = Inf, width = Inf)

## =========================================================
## 12. Incremental tests
## =========================================================

lrt_rows <- list()

if (!is.null(glm_results$env$model) && !is.null(glm_results$env_plus_poll_broad$model)) {
  lrt_rows$glm_poll_after_env <- as.data.frame(
    stats::anova(glm_results$env$model, glm_results$env_plus_poll_broad$model, test = "Chisq")
  ) %>%
    tibble::rownames_to_column("step") %>%
    tibble::as_tibble() %>%
    dplyr::mutate(test = "binomial_GLM_broad_pollinator_after_env")
}

if (!is.null(glm_results$env_plus_space$model) && !is.null(glm_results$env_plus_poll_broad_plus_space$model)) {
  lrt_rows$glm_poll_after_env_space <- as.data.frame(
    stats::anova(glm_results$env_plus_space$model, glm_results$env_plus_poll_broad_plus_space$model, test = "Chisq")
  ) %>%
    tibble::rownames_to_column("step") %>%
    tibble::as_tibble() %>%
    dplyr::mutate(test = "binomial_GLM_broad_pollinator_after_env_space")
}

if (!is.null(bb_results$env$model) && !is.null(bb_results$env_plus_poll_broad$model)) {
  lrt_rows$bb_poll_after_env <- as.data.frame(
    stats::anova(bb_results$env$model, bb_results$env_plus_poll_broad$model)
  ) %>%
    tibble::rownames_to_column("step") %>%
    tibble::as_tibble() %>%
    dplyr::mutate(test = "beta_binomial_broad_pollinator_after_env")
}

if (!is.null(bb_results$env_plus_space$model) && !is.null(bb_results$env_plus_poll_broad_plus_space$model)) {
  lrt_rows$bb_poll_after_env_space <- as.data.frame(
    stats::anova(bb_results$env_plus_space$model, bb_results$env_plus_poll_broad_plus_space$model)
  ) %>%
    tibble::rownames_to_column("step") %>%
    tibble::as_tibble() %>%
    dplyr::mutate(test = "beta_binomial_broad_pollinator_after_env_space")
}

lrt_tbl <- dplyr::bind_rows(lrt_rows)
readr::write_csv(lrt_tbl, file.path(OUT_DIR, "09_incremental_tests_broad_pollinator_after_env_space.csv"))

message("\nIncremental tests:")
print(lrt_tbl, n = Inf, width = Inf)

## =========================================================
## 13. Spatial block CV
## =========================================================

make_blocks <- function(df, block_size) {
  df %>%
    dplyr::mutate(
      block_lon = floor(grid_lon / block_size),
      block_lat = floor(grid_lat / block_size),
      block_id = paste(block_lon, block_lat, sep = "_")
    )
}

dat_cv <- make_blocks(analysis_grid, BLOCK_SIZE_DEG)
blocks <- unique(dat_cv$block_id)

cv_model_specs <- list(
  env = rhs_env,
  env_plus_poll_broad = rhs_env_poll,
  env_plus_space = rhs_env_space,
  env_plus_poll_broad_plus_space = rhs_env_poll_space
)

cv_records <- list()
counter <- 1

for (rep_i in seq_len(N_CV_REPEAT)) {
  shuffled_blocks <- sample(blocks)
  fold_id <- rep(seq_len(5), length.out = length(shuffled_blocks))
  fold_map <- tibble::tibble(block_id = shuffled_blocks, fold = fold_id)

  dd <- dat_cv %>% dplyr::left_join(fold_map, by = "block_id")

  for (fold in sort(unique(dd$fold))) {
    train <- dd %>% dplyr::filter(fold != !!fold)
    test <- dd %>% dplyr::filter(fold == !!fold)

    if (nrow(train) < 20 || nrow(test) < 5) next

    for (mn in names(cv_model_specs)) {
      form <- make_formula(cv_model_specs[[mn]])

      fit <- tryCatch(stats::glm(form, data = train, family = stats::binomial), error = function(e) NULL)
      if (is.null(fit)) next

      pred <- tryCatch(stats::predict(fit, newdata = test, type = "response"), error = function(e) rep(NA_real_, nrow(test)))
      if (all(!is.finite(pred))) next

      cv_records[[counter]] <- tibble::tibble(
        repeat_id = rep_i,
        fold = fold,
        model = mn,
        n_test = nrow(test),
        logloss = logloss_binomial(test$n_nodding_species, test$n_total_species, pred),
        brier = brier_prop(test$prop_nodding_species, pred),
        AUC_presence = safe_auc(ifelse(test$n_nodding_species > 0, 1, 0), pred)
      )

      counter <- counter + 1
    }
  }
}

cv_tbl <- dplyr::bind_rows(cv_records)

if (nrow(cv_tbl) > 0) {
  cv_summary <- cv_tbl %>%
    dplyr::group_by(model) %>%
    dplyr::summarise(
      n_folds = dplyr::n(),
      mean_logloss = mean(logloss, na.rm = TRUE),
      se_logloss = sd(logloss, na.rm = TRUE) / sqrt(dplyr::n()),
      mean_brier = mean(brier, na.rm = TRUE),
      se_brier = sd(brier, na.rm = TRUE) / sqrt(dplyr::n()),
      mean_AUC_presence = mean(AUC_presence, na.rm = TRUE),
      se_AUC_presence = sd(AUC_presence, na.rm = TRUE) / sqrt(sum(is.finite(AUC_presence))),
      .groups = "drop"
    ) %>%
    dplyr::arrange(.data$mean_logloss)
} else {
  cv_summary <- tibble::tibble(
    model = character(),
    n_folds = integer(),
    mean_logloss = numeric(),
    se_logloss = numeric(),
    mean_brier = numeric(),
    se_brier = numeric(),
    mean_AUC_presence = numeric(),
    se_AUC_presence = numeric()
  )
}

readr::write_csv(cv_tbl, file.path(OUT_DIR, "10_spatial_block_CV_raw.csv"))
readr::write_csv(cv_summary, file.path(OUT_DIR, "10_spatial_block_CV_summary.csv"))

message("\nSpatial block CV summary:")
print(cv_summary, n = Inf, width = Inf)

## =========================================================
## 14. Best model coefficients
## =========================================================

best_method <- "beta_binomial_glmmTMB"
best_table <- bb_summary

if (all(is.na(best_table$AIC))) {
  best_method <- "binomial_GLM"
  best_table <- glm_summary
}

best_table_nonNA <- best_table %>% dplyr::filter(is.finite(.data$AIC))

if (nrow(best_table_nonNA) > 0) {
  best_name <- best_table_nonNA$model[1]
} else {
  best_method <- "binomial_GLM"
  best_table_nonNA <- glm_summary %>% dplyr::filter(is.finite(.data$AIC))
  best_name <- best_table_nonNA$model[1]
}

best_coef <- if (best_method == "beta_binomial_glmmTMB") {
  bb_coef %>% dplyr::filter(model == best_name)
} else {
  glm_coef %>% dplyr::filter(model == best_name)
}

readr::write_csv(best_coef, file.path(OUT_DIR, "11_best_model_coefficients.csv"))

message("\nBest model by AIC:")
print(best_table_nonNA[1, ], width = Inf)

message("\nBest model coefficients:")
print(best_coef, n = Inf, width = Inf)

## =========================================================
## 15. Figures
## =========================================================

plot_model_compare <- function(df, title, out_file) {
  if (!"AIC" %in% names(df) || all(is.na(df$AIC))) return(NULL)

  plot_df <- df %>%
    dplyr::filter(is.finite(.data$AIC)) %>%
    dplyr::arrange(.data$AIC)

  if (nrow(plot_df) == 0) return(NULL)

  p <- plot_df %>%
    dplyr::mutate(model = factor(model, levels = rev(model))) %>%
    ggplot2::ggplot(ggplot2::aes(x = delta_AIC, y = model)) +
    ggplot2::geom_col() +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::labs(x = "ΔAIC", y = NULL, title = title)

  ggplot2::ggsave(out_file, p, width = 8, height = 5, dpi = 300)
}

plot_model_compare(glm_summary, "Binomial GLM model comparison", file.path(FIG_DIR, "model_compare_binomial_GLM.png"))
plot_model_compare(bb_summary, "Beta-binomial model comparison", file.path(FIG_DIR, "model_compare_beta_binomial.png"))
plot_model_compare(gam_summary, "Spatial GAM model comparison", file.path(FIG_DIR, "model_compare_spatial_GAM.png"))

if (nrow(cv_summary) > 0) {
  p_cv <- cv_summary %>%
    dplyr::mutate(model = factor(model, levels = rev(model))) %>%
    ggplot2::ggplot(ggplot2::aes(x = mean_logloss, y = model)) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_errorbarh(
      ggplot2::aes(xmin = mean_logloss - se_logloss, xmax = mean_logloss + se_logloss),
      height = 0.15
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::labs(x = "Mean spatial-block CV log loss", y = NULL, title = "Spatial block cross-validation")

  ggplot2::ggsave(file.path(FIG_DIR, "spatial_block_CV_logloss.png"), p_cv, width = 8, height = 5, dpi = 300)
}

if (nrow(best_coef) > 0 && "estimate" %in% names(best_coef)) {
  coef_plot <- best_coef %>%
    dplyr::filter(term != "(Intercept)") %>%
    dplyr::mutate(
      conf.low = estimate - 1.96 * std.error,
      conf.high = estimate + 1.96 * std.error,
      term = forcats::fct_reorder(term, estimate)
    )

  if (nrow(coef_plot) > 0) {
    p_coef <- ggplot2::ggplot(coef_plot, ggplot2::aes(x = estimate, y = term)) +
      ggplot2::geom_vline(xintercept = 0, linetype = 2) +
      ggplot2::geom_errorbarh(ggplot2::aes(xmin = conf.low, xmax = conf.high), height = 0.15) +
      ggplot2::geom_point(size = 2) +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::labs(
        x = "Coefficient",
        y = NULL,
        title = paste0("Best model: ", best_name, " / ", best_method)
      )

    ggplot2::ggsave(file.path(FIG_DIR, "best_model_coefficients_broad_group.png"), p_coef, width = 9, height = 6, dpi = 300)
  }
}

p_map <- ggplot2::ggplot(
  analysis_grid,
  ggplot2::aes(grid_lon, grid_lat, color = prop_nodding_species, size = n_total_species)
) +
  ggplot2::geom_point(alpha = 0.85) +
  ggplot2::coord_equal() +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::scale_color_viridis_c(option = "C", limits = c(0, 1)) +
  ggplot2::labs(
    x = "Longitude",
    y = "Latitude",
    color = "Proportion\nnodding",
    size = "Total\nspecies",
    title = "Observed regional proportion of nodding species"
  )

ggplot2::ggsave(file.path(FIG_DIR, "map_observed_prop_nodding_broad_group.png"), p_map, width = 7, height = 7, dpi = 300)

## =========================================================
## 16. Done
## =========================================================

message("\n🎉 DONE — broad bee / butterfly / fly SSDM analysis")
message("Output dir: ", OUT_DIR)
message("Key outputs:")
message("  00_pollinator_broad_group_classification.csv")
message("  01_broad_group_ssdm_stack_summary.csv")
message("  04_analysis_grid_data_broad_group_SSDM.csv")
message("  06_binomial_GLM_model_comparison.csv")
message("  07_beta_binomial_model_comparison.csv")
message("  08_spatial_GAM_model_comparison.csv")
message("  08_spatial_GAM_parametric_coefficients.csv")
message("  08_spatial_GAM_smooth_terms.csv")
message("  09_incremental_tests_broad_pollinator_after_env_space.csv")
message("  10_spatial_block_CV_summary.csv")
message("  11_best_model_coefficients.csv")
message("Figures:")
message("  figures/model_compare_*.png")
message("  figures/spatial_block_CV_logloss.png")
message("  figures/best_model_coefficients_broad_group.png")
message("  figures/map_observed_prop_nodding_broad_group.png")
