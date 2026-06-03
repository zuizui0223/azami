############################################################
## FULL PIPELINE FROM SCRATCH
## CHELSA variables excluding SWE -> crop Japan -> VIFSTEP10
## -> ENMeval/maxnet species-specific pollinator SDMs
## -> guild-level stacked pollinator surfaces
## -> extract guild SDMs to Cirsium occurrence/species table
##
## Main reason for rebuilding:
##   SWE had very high NA rate in occurrence points and caused ENMeval
##   to remove most occurrences. Therefore SWE is excluded from candidates
##   BEFORE VIFSTEP.
##
## Reference workflow:
##   Based on user's sdm2.R:
##     - prepare aligned environmental stack
##     - VIFSTEP on continuous variables
##     - ENMevaluate(algorithm = "maxnet")
##     - select best model
##     - predict cloglog map
##     - save response curves and metadata
##
## Input:
##   pollinator_gbif_outputs/gbif_pollinator_occurrences_20plus_unique_coords.csv
##
## Optional:
##   cirsium_head_orientation_worldclim_ssdm_outputs/01_cirsium_occurrences_with_worldclim_ssdm.csv
##
## Output:
##   chelsa_pollinator_enmeval_rebuild_no_swe/
############################################################

## =========================================================
## 0. Packages
## =========================================================

pkgs <- c(
  "terra", "raster", "geodata", "sf",
  "dplyr", "readr", "stringr", "tibble", "purrr", "tidyr",
  "usdm", "ENMeval", "maxnet", "tools"
)

for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, dependencies = TRUE)
  }
}

suppressPackageStartupMessages({
  library(terra)
  library(raster)
  library(geodata)
  library(sf)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
  library(purrr)
  library(tidyr)
  library(usdm)
  library(ENMeval)
  library(maxnet)
  library(tools)
})

set.seed(42)

## =========================================================
## 1. Global settings
## =========================================================

ROOT_DIR <- "chelsa_pollinator_enmeval_rebuild_no_swe"

## IMPORTANT:
## Do not download CHELSA again.
## Reuse the files already downloaded/cropped by the previous pipeline.
PREV_ROOT_DIR <- "chelsa_pollinator_head_orientation_pipeline"
PREV_ENV_RAW_DIR <- file.path(PREV_ROOT_DIR, "env_raw_downloads")
PREV_ENV_CROP_DIR <- file.path(PREV_ROOT_DIR, "env_candidates", "cropped_layers_core")

## New outputs go here, but raw/cropped inputs are reused from PREV_*.
ENV_RAW_DIR <- PREV_ENV_RAW_DIR
ENV_CROP_DIR <- file.path(ROOT_DIR, "env_cropped_layers_no_swe")
ENV_CAND_DIR <- file.path(ROOT_DIR, "env_candidates")
ENV_VIF_DIR <- file.path(ROOT_DIR, "env_vif")
OCC_OUT_DIR <- file.path(ROOT_DIR, "occurrences")
SPECIES_OUT_DIR <- file.path(ROOT_DIR, "species_models")
GROUP_OUT_DIR <- file.path(ROOT_DIR, "guild_stacks")
M_AREA_DIR <- file.path(ROOT_DIR, "species_accessible_area")
LOG_DIR <- file.path(ROOT_DIR, "logs")
TMP_DIR <- file.path(ROOT_DIR, "terra_tmp")

dirs <- c(
  ROOT_DIR, ENV_RAW_DIR, ENV_CROP_DIR, ENV_CAND_DIR, ENV_VIF_DIR,
  OCC_OUT_DIR, SPECIES_OUT_DIR, GROUP_OUT_DIR, M_AREA_DIR, LOG_DIR, TMP_DIR
)
for (d in dirs) dir.create(d, showWarnings = FALSE, recursive = TRUE)

terra::terraOptions(memfrac = 0.45, tempdir = TMP_DIR)

## Input files
VISITOR_FILE <- file.path(
  "pollinator_gbif_outputs",
  "gbif_pollinator_occurrences_20plus_unique_coords.csv"
)

CIRSIUM_FILE <- file.path(
  "cirsium_head_orientation_worldclim_ssdm_outputs",
  "01_cirsium_occurrences_with_worldclim_ssdm.csv"
)

## CHELSA settings
CHELSA_BASE_URL <- "https://os.unil.cloud.switch.ch/chelsa02/chelsa/global/bioclim"
CHELSA_PERIOD <- "1981-2010"
CHELSA_VERSION <- "V.2.1"

## Japan extent
JAPAN_EXT <- terra::ext(122, 146.5, 24, 46)

## Rebuild switches
## This version never downloads CHELSA again.
FORCE_REDOWNLOAD_CHELSA <- FALSE
FORCE_REBUILD_CROPPED_LAYERS <- FALSE
FORCE_REBUILD_VIF_STACK <- TRUE

## Use already cropped layers from previous pipeline whenever possible.
USE_PREVIOUS_CROPPED_LAYERS <- TRUE
ALLOW_DOWNLOAD_IF_MISSING <- FALSE
FORCE_RERUN_SPECIES <- FALSE
FORCE_REBUILD_GROUP_STACKS <- TRUE

## VIF settings
VIF_THRESHOLD <- 10
VIF_SAMPLE_SIZE <- 10000

## ENMeval settings
BG_SIZE_PER_SPECIES <- 10000
ENMEVAL_CORES <- 4
ENMEVAL_PARALLEL <- TRUE

TUNE_FC <- c("L", "LQ", "LQH")
TUNE_RM <- 1:5

MIN_UNIQUE_CELLS_PER_SPECIES <- 20
MAX_OCC_PER_SPECIES <- Inf

## Species-specific accessible area / M
BUFFER_KM_CANDIDATES <- c(50, 100, 150, 200)
MIN_BG_POINTS <- 500
BUFFER_CRS <- "EPSG:3857"

USE_BLOCK_FIRST <- TRUE
RANDOM_KFOLDS <- 5

## Known problematic variables to exclude BEFORE VIFSTEP.
## SWE caused massive NA rates, so it is not downloaded/used.
EXCLUDE_ENV_VARS <- c("swe")

## If any remaining layer has too much NA for a species, drop it only for that species model.
MAX_LAYER_NA_PROP_OCC <- 0.20
MIN_OCC_AFTER_NA_FILTER <- 5

## Outputs
ENV_CANDIDATE_FILELIST <- file.path(ENV_CAND_DIR, "chelsa_candidates_no_swe_filelist.csv")
ENV_CANDIDATE_STACK_TIF <- file.path(ENV_CAND_DIR, "chelsa_candidates_no_swe_japan.tif")

ENV_VIF_TIF <- file.path(ENV_VIF_DIR, "chelsa_vif10_no_swe_japan.tif")
ENV_VIF_GRD <- file.path(ENV_VIF_DIR, "chelsa_vif10_no_swe_japan_rasterstack.grd")
ENV_VIF_SELECTED_CSV <- file.path(ENV_VIF_DIR, "vifstep10_selected_variables_no_swe.csv")

## =========================================================
## 2. Helpers
## =========================================================

safe_name <- function(x) {
  x |>
    as.character() |>
    stringr::str_replace_all("[^A-Za-z0-9_\\-]+", "_") |>
    stringr::str_replace_all("_+", "_") |>
    stringr::str_replace_all("^_|_$", "")
}

stop_if_missing <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)
}

to_numeric_coord <- function(x) {
  if (is.list(x)) x <- unlist(x)
  x <- as.character(x)
  x <- stringr::str_replace_all(x, ",", ".")
  suppressWarnings(as.numeric(x))
}

make_xy_df <- function(lon, lat) {
  data.frame(
    lon = to_numeric_coord(lon),
    lat = to_numeric_coord(lat)
  ) |>
    dplyr::filter(is.finite(lon), is.finite(lat))
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

safe_extract_first <- function(r, xy) {
  xy <- as.matrix(xy[, c("lon", "lat"), drop = FALSE])
  out <- safe_extract_df(r, xy)
  as.numeric(out[[1]])
}

chelsa_url <- function(var) {
  paste0(
    CHELSA_BASE_URL, "/", var, "/", CHELSA_PERIOD,
    "/CHELSA_", var, "_", CHELSA_PERIOD, "_", CHELSA_VERSION, ".tif"
  )
}

download_one <- function(url, dest, force = FALSE) {
  if (!force && file.exists(dest) && file.info(dest)$size > 1000) {
    message("Cached: ", basename(dest))
    return(TRUE)
  }

  if (force && file.exists(dest)) unlink(dest)

  message("Downloading: ", basename(dest))

  ok <- tryCatch({
    utils::download.file(url, dest, mode = "wb", quiet = FALSE)
    file.exists(dest) && file.info(dest)$size > 1000
  }, error = function(e) {
    message("FAILED: ", basename(dest), " | ", conditionMessage(e))
    FALSE
  })

  if (!ok && file.exists(dest)) unlink(dest)
  ok
}

crop_mask_write_layer <- function(src_file, clean_name, jpn, out_dir, template_file = NULL, force = FALSE) {
  out_file <- file.path(out_dir, paste0(clean_name, "_japan.tif"))

  if (!force && file.exists(out_file) && file.info(out_file)$size > 1000) {
    message("Cropped cached: ", basename(out_file))
    return(out_file)
  }

  if (force && file.exists(out_file)) unlink(out_file)

  message("Cropping/masking: ", clean_name)

  r <- terra::rast(src_file)
  jpn2 <- terra::project(jpn, terra::crs(r))

  r <- terra::crop(r, JAPAN_EXT)
  r <- terra::crop(r, jpn2)
  r <- terra::mask(r, jpn2)

  if (!is.null(template_file) && file.exists(template_file)) {
    template <- terra::rast(template_file)
    if (!terra::compareGeom(template, r, stopOnError = FALSE)) {
      if (terra::crs(r) == terra::crs(template)) {
        r <- terra::resample(r, template, method = "bilinear")
      } else {
        r <- terra::project(r, template, method = "bilinear")
      }
    }
  }

  names(r) <- clean_name

  terra::writeRaster(
    r,
    out_file,
    overwrite = TRUE,
    wopt = list(gdal = c("COMPRESS=LZW", "TILED=YES"), datatype = "FLT4S")
  )

  rm(r)
  gc(verbose = FALSE)

  out_file
}

standardize_occ_columns <- function(df) {
  if (!"species" %in% names(df)) {
    if ("query_name" %in% names(df)) {
      df <- dplyr::mutate(df, species = query_name)
    } else if ("scientificName" %in% names(df)) {
      df <- dplyr::mutate(df, species = scientificName)
    } else if ("visitor_name" %in% names(df)) {
      df <- dplyr::mutate(df, species = visitor_name)
    } else {
      stop("No species/query_name/scientificName/visitor_name column found.")
    }
  }

  if (!"longitude" %in% names(df)) {
    if ("decimalLongitude" %in% names(df)) {
      df <- dplyr::mutate(df, longitude = decimalLongitude)
    } else {
      stop("No longitude/decimalLongitude column found.")
    }
  }

  if (!"latitude" %in% names(df)) {
    if ("decimalLatitude" %in% names(df)) {
      df <- dplyr::mutate(df, latitude = decimalLatitude)
    } else {
      stop("No latitude/decimalLatitude column found.")
    }
  }

  df
}

classify_visitor_fine <- function(species_name, previous_group = NA_character_) {
  genus <- stringr::str_extract(as.character(species_name), "^[A-Za-z]+")

  bee_hanging_effective <- c("Bombus", "Xylocopa", "Anthophora", "Eucera")

  bee_non_hanging <- c(
    "Apis", "Megachile", "Osmia", "Andrena", "Colletes", "Nomada", "Ceratina",
    "Lasioglossum", "Halictus", "Hylaeus", "Sphecodes", "Seladonia"
  )

  butterfly <- c(
    "Papilio", "Pieris", "Colias", "Gonepteryx",
    "Vanessa", "Aglais", "Nymphalis", "Polygonia",
    "Argynnis", "Fabriciana", "Speyeria", "Cynthia",
    "Aporia", "Erebia", "Araschnia", "Parnassius",
    "Lycaena", "Celastrina", "Zizeeria", "Lampides",
    "Ypthima", "Lethe", "Neptis", "Melitaea",
    "Erynnis", "Hesperia", "Ochlodes"
  )

  hoverfly <- c(
    "Episyrphus", "Eristalis", "Melanostoma", "Syrphus",
    "Sphaerophoria", "Eupeodes", "Helophilus", "Volucella",
    "Metasyrphus", "Paragus", "Cheilosia", "Platycheirus"
  )

  hawkmoth <- c("Macroglossum", "Theretra", "Deilephila", "Sphinx", "Agrius")

  wasp <- c(
    "Vespa", "Vespula", "Polistes", "Eumenes",
    "Ammophila", "Sphex", "Cerceris", "Ichneumon", "Pimpla"
  )

  beetle <- c(
    "Cetonia", "Protaetia", "Oxythyrea", "Mordella",
    "Cantharis", "Cryptocephalus", "Chrysolina", "Coccinella"
  )

  other_fly <- c(
    "Lucilia", "Calliphora", "Sarcophaga", "Musca",
    "Anthomyia", "Tachina", "Phormia", "Scathophaga",
    "Empis", "Bibio"
  )

  dplyr::case_when(
    genus %in% bee_hanging_effective ~ "bee_hanging_effective",
    genus %in% bee_non_hanging ~ "bee_non_hanging",
    genus %in% butterfly ~ "butterfly",
    genus %in% hoverfly ~ "hoverfly",
    genus %in% hawkmoth ~ "hawkmoth",
    genus %in% wasp ~ "wasp",
    genus %in% beetle ~ "beetle",
    genus %in% other_fly ~ "other_fly",
    !is.na(previous_group) & previous_group %in% c("lepidoptera", "butterfly") ~ "butterfly",
    !is.na(previous_group) & previous_group %in% c("hoverfly") ~ "hoverfly",
    !is.na(previous_group) & previous_group %in% c("bee_or_wasp", "bee") ~ "bee_unknown",
    TRUE ~ "unknown"
  )
}

select_best <- function(res) {
  if ("cbi.val.avg" %in% names(res) && any(is.finite(res$cbi.val.avg))) {
    max_cbi <- max(res$cbi.val.avg, na.rm = TRUE)
    cand <- which(res$cbi.val.avg == max_cbi)

    if ("auc.val.avg" %in% names(res)) {
      auc_sub <- res$auc.val.avg[cand]
      if (any(is.finite(auc_sub))) cand <- cand[which.max(auc_sub)]
    }

    if ("or.10p.val.avg" %in% names(res)) {
      or_sub <- res$or.10p.val.avg[cand]
      if (any(is.finite(or_sub))) cand <- cand[which.min(or_sub)]
    }

    return(cand[1])
  }

  if ("AICc" %in% names(res) && any(is.finite(res$AICc))) return(which.min(res$AICc))
  if ("auc.val.avg" %in% names(res) && any(is.finite(res$auc.val.avg))) return(which.max(res$auc.val.avg))

  1L
}

predict_maxnet_rast <- function(env_rast, model, categoricals = NULL) {
  terra::predict(
    env_rast,
    model,
    fun = function(m, newdata, ...) {
      if (!is.null(categoricals)) {
        for (cn in categoricals) {
          if (cn %in% colnames(newdata)) {
            newdata[[cn]] <- factor(as.character(newdata[[cn]]))
          }
        }
      }
      as.numeric(predict(m, newdata = as.data.frame(newdata), type = "cloglog"))
    },
    na.rm = TRUE
  )
}

make_species_m_area <- function(occ_xy, env_full, buffer_km, species_safe) {
  pts <- terra::vect(
    occ_xy,
    geom = c("lon", "lat"),
    crs = "EPSG:4326"
  )

  pts_m <- terra::project(pts, BUFFER_CRS)
  buff_m <- terra::buffer(pts_m, width = buffer_km * 1000)
  buff_m <- terra::aggregate(buff_m)
  buff_ll <- terra::project(buff_m, terra::crs(env_full))

  m_file <- file.path(M_AREA_DIR, paste0(species_safe, "_M_buffer_", buffer_km, "km.gpkg"))
  try(terra::writeVector(buff_ll, m_file, overwrite = TRUE), silent = TRUE)

  env_m <- terra::crop(env_full, buff_ll)
  env_m <- terra::mask(env_m, buff_ll)

  list(poly = buff_ll, env = env_m, file = m_file)
}

extend_to_template <- function(r, template) {
  out <- terra::extend(r, template)
  if (!terra::compareGeom(out, template, stopOnError = FALSE)) {
    out <- terra::resample(out, template, method = "bilinear")
  }
  out
}

screen_env_layers_for_species <- function(env_m, occ_xy, sp_dir, sp_safe, log_msg) {
  vals <- safe_extract_df(env_m, as.matrix(occ_xy[, c("lon", "lat"), drop = FALSE]))

  na_by_layer <- tibble::tibble(
    layer = names(vals),
    n_NA = colSums(is.na(vals)),
    prop_NA = colMeans(is.na(vals))
  ) |>
    dplyr::arrange(dplyr::desc(prop_NA))

  readr::write_csv(
    na_by_layer,
    file.path(sp_dir, paste0(sp_safe, "_NA_by_layer_before_ENMeval.csv"))
  )

  bad_layers <- na_by_layer |>
    dplyr::filter(prop_NA > MAX_LAYER_NA_PROP_OCC) |>
    dplyr::pull(layer)

  bad_layers <- intersect(bad_layers, names(env_m))

  if (length(bad_layers) > 0) {
    log_msg("Dropping high-NA predictors for this species: ", paste(bad_layers, collapse = ", "))
    env_m <- env_m[[setdiff(names(env_m), bad_layers)]]
    vals <- vals[, setdiff(names(vals), bad_layers), drop = FALSE]
  }

  ok <- stats::complete.cases(vals)

  log_msg("Occurrences before predictor-NA filtering: ", nrow(occ_xy))
  log_msg("Occurrences after predictor-NA filtering: ", sum(ok))
  log_msg("Removed by remaining NA predictors: ", sum(!ok))

  readr::write_csv(
    tibble::tibble(
      species = sp_safe,
      n_occ_before = nrow(occ_xy),
      n_occ_after = sum(ok),
      n_removed = sum(!ok),
      dropped_layers = paste(bad_layers, collapse = ";")
    ),
    file.path(sp_dir, paste0(sp_safe, "_NA_filter_summary.csv"))
  )

  list(
    env = env_m,
    occ = occ_xy[ok, , drop = FALSE],
    dropped_layers = bad_layers,
    na_by_layer = na_by_layer
  )
}

## =========================================================
## 3. CHELSA variable table, excluding SWE before VIF
## =========================================================

vars_core <- tibble::tribble(
  ~var, ~hypothesis, ~meaning, ~expected_nodding_direction,
  "bio01", "cool_montane", "Annual mean temperature", "lower",
  "bio04", "seasonality", "Temperature seasonality", "unknown",
  "bio06", "cold_snow", "Minimum temperature of coldest month", "lower",
  "bio10", "summer_temperature", "Mean temperature of warmest quarter", "lower",
  "bio11", "cold_snow", "Mean temperature of coldest quarter", "lower",
  "bio13", "rain_wet_protection", "Precipitation of wettest month", "higher",
  "bio15", "precipitation_seasonality", "Precipitation seasonality", "unknown",
  "bio18", "flowering_rain_proxy", "Precipitation of warmest quarter", "higher",
  "bio19", "snow_cold_season_precip", "Precipitation of coldest quarter", "higher",

  "cltmean", "cloud_wetness", "Mean cloud area fraction", "higher",
  "cmimean", "water_balance", "Mean climate moisture index", "higher",
  "hursmean", "humidity_wetness", "Mean relative humidity", "higher",
  "vpdmean", "drying_power", "Mean vapor pressure deficit", "lower",
  "petmean", "evapotranspiration", "Mean potential evapotranspiration", "lower_or_unknown",
  "rsdsmean", "radiation_flower_temperature", "Mean shortwave radiation", "unknown",
  "sfcWindmean", "wind_exposure", "Mean near-surface wind speed", "higher",

  "scd", "snow", "Snow cover days", "higher",
  ## SWE intentionally excluded
  ## "swe", "snow", "Snow water equivalent", "higher",
  "fcf", "frost", "Frost change frequency", "higher",
  "gdd5", "growing_degree_days", "Growing degree days above 5C", "lower",
  "gsl", "growing_season", "Growing season length", "shorter",
  "gsp", "growing_season_precip", "Growing season precipitation", "higher",
  "gst", "growing_season_temperature", "Growing season temperature", "lower",
  "swb", "site_water_balance", "Site water balance", "higher"
) |>
  dplyr::filter(!var %in% EXCLUDE_ENV_VARS)

CHELSA_VARS_TBL <- vars_core |>
  dplyr::mutate(
    url = chelsa_url(var),
    local_file = file.path(
      ENV_RAW_DIR,
      paste0("CHELSA_", var, "_", CHELSA_PERIOD, "_", CHELSA_VERSION, ".tif")
    ),
    clean_name = stringr::str_replace(var, "^bio", "BIO")
  )

readr::write_csv(CHELSA_VARS_TBL, file.path(ROOT_DIR, "selected_chelsa_variables_no_swe.csv"))

## =========================================================
## 4. Build candidate env stack from scratch
## =========================================================

message("\n=== Preparing environmental layers from existing files ===")
jpn <- NULL

message("\n=== Reusing existing CHELSA cropped/raw files, excluding SWE ===")

cropped_files <- character(0)
cropped_names <- character(0)
download_log <- list()
template_file <- NULL

for (i in seq_len(nrow(CHELSA_VARS_TBL))) {
  v <- CHELSA_VARS_TBL$var[i]
  clean_nm <- CHELSA_VARS_TBL$clean_name[i]
  url <- CHELSA_VARS_TBL$url[i]
  dest <- CHELSA_VARS_TBL$local_file[i]

  cropped_file <- NA_character_
  ok <- FALSE
  source_type <- NA_character_

  ## 1) Prefer already-cropped Japan layer from previous pipeline.
  prev_crop <- file.path(PREV_ENV_CROP_DIR, paste0(clean_nm, "_japan.tif"))

  if (USE_PREVIOUS_CROPPED_LAYERS && file.exists(prev_crop) && file.info(prev_crop)$size > 1000) {
    message("Using previous cropped layer: ", basename(prev_crop))
    cropped_file <- prev_crop
    ok <- TRUE
    source_type <- "previous_cropped"
  }

  ## 2) If cropped layer is not available, use already-downloaded raw CHELSA file and crop/mask.
  if (!ok && file.exists(dest) && file.info(dest)$size > 1000) {
    message("Using previous raw download and cropping: ", basename(dest))
    if (is.null(jpn)) {
      jpn <- geodata::gadm(country = "JPN", level = 0, path = file.path(ENV_RAW_DIR, "gadm"))
    }
    cropped_file <- crop_mask_write_layer(
      src_file = dest,
      clean_name = clean_nm,
      jpn = jpn,
      out_dir = ENV_CROP_DIR,
      template_file = template_file,
      force = FORCE_REBUILD_CROPPED_LAYERS
    )
    ok <- TRUE
    source_type <- "previous_raw_cropped_now"
  }

  ## 3) Optional download fallback. Default is FALSE.
  if (!ok && ALLOW_DOWNLOAD_IF_MISSING) {
    ok_download <- download_one(url, dest, force = FORCE_REDOWNLOAD_CHELSA)
    if (ok_download) {
      if (is.null(jpn)) {
        jpn <- geodata::gadm(country = "JPN", level = 0, path = file.path(ENV_RAW_DIR, "gadm"))
      }
      cropped_file <- crop_mask_write_layer(
        src_file = dest,
        clean_name = clean_nm,
        jpn = jpn,
        out_dir = ENV_CROP_DIR,
        template_file = template_file,
        force = FORCE_REBUILD_CROPPED_LAYERS
      )
      ok <- TRUE
      source_type <- "downloaded_now"
    }
  }

  if (!ok) {
    stop(
      "Missing existing CHELSA file for ", clean_nm, "\n",
      "Expected cropped: ", prev_crop, "\n",
      "Expected raw: ", dest, "\n",
      "This no-download script will not fetch it automatically. ",
      "Set ALLOW_DOWNLOAD_IF_MISSING <- TRUE only if you really want to download."
    )
  }

  if (is.null(template_file)) template_file <- cropped_file

  cropped_files <- c(cropped_files, cropped_file)
  cropped_names <- c(cropped_names, clean_nm)

  download_log[[i]] <- tibble::tibble(
    var = v,
    clean_name = clean_nm,
    url = url,
    local_file = dest,
    cropped_file = cropped_file,
    source_type = source_type,
    success = ok
  )
}

readr::write_csv(dplyr::bind_rows(download_log), file.path(ROOT_DIR, "download_log_chelsa_no_swe.csv"))

## Topography
message("\n=== Topography: elevation / slope / roughness ===")
elev_file <- file.path(ENV_CROP_DIR, "elevation_japan.tif")
slope_file <- file.path(ENV_CROP_DIR, "slope_japan.tif")
rough_file <- file.path(ENV_CROP_DIR, "roughness_japan.tif")

topo_exists <- file.exists(elev_file) && file.exists(slope_file) && file.exists(rough_file)

if (!topo_exists || FORCE_REBUILD_CROPPED_LAYERS) {
  prev_elev <- file.path(PREV_ENV_CROP_DIR, "elevation_japan.tif")
  prev_slope <- file.path(PREV_ENV_CROP_DIR, "slope_japan.tif")
  prev_rough <- file.path(PREV_ENV_CROP_DIR, "roughness_japan.tif")

  if (file.exists(prev_elev) && file.exists(prev_slope) && file.exists(prev_rough)) {
    message("Using previous cropped topography layers.")
    elev_file <- prev_elev
    slope_file <- prev_slope
    rough_file <- prev_rough
  } else {
    message("Previous cropped topography not found; using cached geodata elevation if available.")
    if (is.null(jpn)) {
      jpn <- geodata::gadm(country = "JPN", level = 0, path = file.path(ENV_RAW_DIR, "gadm"))
    }

    elev <- geodata::elevation_30s(country = "JPN", path = file.path(ENV_RAW_DIR, "elev"))
    jpn_e <- terra::project(jpn, terra::crs(elev))

    elev <- terra::crop(elev, JAPAN_EXT)
    elev <- terra::crop(elev, jpn_e)
    elev <- terra::mask(elev, jpn_e)

    if (!is.null(template_file) && file.exists(template_file)) {
      template <- terra::rast(template_file)
      if (!terra::compareGeom(template, elev, stopOnError = FALSE)) {
        if (terra::crs(elev) == terra::crs(template)) {
          elev <- terra::resample(elev, template, method = "bilinear")
        } else {
          elev <- terra::project(elev, template, method = "bilinear")
        }
      }
    }

    names(elev) <- "elevation"

    slope <- terra::terrain(elev, v = "slope", unit = "degrees")
    names(slope) <- "slope"

    roughness <- terra::terrain(elev, v = "roughness")
    names(roughness) <- "roughness"

    terra::writeRaster(elev, elev_file, overwrite = TRUE, wopt = list(gdal = c("COMPRESS=LZW", "TILED=YES"), datatype = "FLT4S"))
    terra::writeRaster(slope, slope_file, overwrite = TRUE, wopt = list(gdal = c("COMPRESS=LZW", "TILED=YES"), datatype = "FLT4S"))
    terra::writeRaster(roughness, rough_file, overwrite = TRUE, wopt = list(gdal = c("COMPRESS=LZW", "TILED=YES"), datatype = "FLT4S"))

    rm(elev, slope, roughness)
    gc(verbose = FALSE)
  }
}

cropped_files <- c(cropped_files, elev_file, slope_file, rough_file)
cropped_names <- c(cropped_names, "elevation", "slope", "roughness")

env_candidates <- terra::rast(cropped_files)
names(env_candidates) <- cropped_names

readr::write_csv(
  tibble::tibble(variable = names(env_candidates), cropped_file = cropped_files),
  ENV_CANDIDATE_FILELIST
)

message("Candidate env layers:")
print(names(env_candidates))

## Optional: save full candidate stack
terra::writeRaster(
  env_candidates,
  ENV_CANDIDATE_STACK_TIF,
  overwrite = TRUE,
  wopt = list(gdal = c("COMPRESS=LZW", "TILED=YES"), datatype = "FLT4S")
)

## =========================================================
## 5. VIFSTEP10 after excluding SWE
## =========================================================

message("\n=== Running VIFSTEP10 without SWE ===")

bg_vif <- terra::spatSample(
  env_candidates,
  size = VIF_SAMPLE_SIZE,
  method = "random",
  na.rm = TRUE,
  as.df = TRUE
)

bg_vif <- bg_vif |>
  dplyr::select(where(is.numeric)) |>
  tidyr::drop_na()

var_sd <- vapply(bg_vif, stats::sd, numeric(1), na.rm = TRUE)
keep <- names(var_sd)[is.finite(var_sd) & var_sd > 0]
bg_vif2 <- bg_vif[, keep, drop = FALSE]

vif_res <- usdm::vifstep(bg_vif2, th = VIF_THRESHOLD)
selected_vars <- vif_res@results$Variables

readr::write_csv(as.data.frame(vif_res@results), file.path(ENV_VIF_DIR, "vifstep10_results_no_swe.csv"))
readr::write_csv(tibble::tibble(selected_variable = selected_vars), ENV_VIF_SELECTED_CSV)

selected_summary <- CHELSA_VARS_TBL |>
  dplyr::mutate(selected_by_vif = clean_name %in% selected_vars) |>
  dplyr::select(var, clean_name, hypothesis, meaning, expected_nodding_direction, selected_by_vif)

topo_summary <- tibble::tibble(
  var = c("elevation", "slope", "roughness"),
  clean_name = c("elevation", "slope", "roughness"),
  hypothesis = c("topography", "wind_slope_exposure", "terrain_complexity"),
  meaning = c("Elevation", "Slope", "Terrain roughness"),
  expected_nodding_direction = c("higher_or_unknown", "higher", "higher"),
  selected_by_vif = c("elevation", "slope", "roughness") %in% selected_vars
)

readr::write_csv(
  dplyr::bind_rows(selected_summary, topo_summary),
  file.path(ENV_VIF_DIR, "selected_variable_hypothesis_summary_no_swe.csv")
)

message("Selected variables by VIF:")
print(selected_vars)

env <- env_candidates[[selected_vars]]

terra::writeRaster(
  env,
  ENV_VIF_TIF,
  overwrite = TRUE,
  wopt = list(gdal = c("COMPRESS=LZW", "TILED=YES"), datatype = "FLT4S")
)

Env_r <- raster::stack(ENV_VIF_TIF)
raster::writeRaster(Env_r, filename = ENV_VIF_GRD, overwrite = TRUE)

message("Saved VIF stack:")
message("  ", ENV_VIF_TIF)
message("  ", ENV_VIF_GRD)

## Reload with terra
env <- terra::rast(ENV_VIF_TIF)

## Ensure no SWE
if ("swe" %in% names(env)) stop("swe is still in env stack. Check exclusion.")

message("Final env layers:")
print(names(env))

## Categorical detection
cat_idx <- grep("worldcover|landcover|categor", names(env), ignore.case = TRUE)
if (length(cat_idx) == 0) {
  categoricals <- NULL
} else {
  categoricals <- names(env)[cat_idx]
  for (cn in categoricals) env[[cn]] <- as.factor(env[[cn]])
}

message("Categoricals: ", ifelse(is.null(categoricals), "(none)", paste(categoricals, collapse = ", ")))

## =========================================================
## 6. Prepare visitor occurrences
## =========================================================

message("\n=== Preparing visitor occurrences ===")
stop_if_missing(VISITOR_FILE)

visitor_df <- readr::read_csv(VISITOR_FILE, show_col_types = FALSE) |>
  standardize_occ_columns() |>
  dplyr::mutate(
    longitude = to_numeric_coord(longitude),
    latitude = to_numeric_coord(latitude)
  )

previous_group_col <- if ("pollinator_group" %in% names(visitor_df)) "pollinator_group" else NA_character_

visitor_df <- visitor_df |>
  dplyr::filter(!is.na(species), species != "") |>
  dplyr::filter(is.finite(longitude), is.finite(latitude)) |>
  dplyr::filter(longitude >= 122, longitude <= 146.5, latitude >= 24, latitude <= 46) |>
  dplyr::mutate(
    previous_group = if (!is.na(previous_group_col)) as.character(.data[[previous_group_col]]) else NA_character_,
    guild = mapply(classify_visitor_fine, species, previous_group) |> as.character()
  ) |>
  dplyr::filter(guild != "unknown") |>
  dplyr::distinct(guild, species, longitude, latitude, .keep_all = TRUE)

readr::write_csv(visitor_df, file.path(OCC_OUT_DIR, "visitor_occurrences_classified.csv"))

xy_mat <- cbind(visitor_df$longitude, visitor_df$latitude)
colnames(xy_mat) <- c("lon", "lat")
visitor_df$cell_id <- terra::cellFromXY(env, xy_mat)

visitor_clean <- visitor_df |>
  dplyr::filter(!is.na(cell_id)) |>
  dplyr::distinct(guild, species, cell_id, .keep_all = TRUE)

if (is.finite(MAX_OCC_PER_SPECIES)) {
  visitor_clean <- visitor_clean |>
    dplyr::group_by(guild, species) |>
    dplyr::group_modify(function(.x, .y) {
      n_take <- min(nrow(.x), MAX_OCC_PER_SPECIES)
      dplyr::slice_sample(.x, n = n_take)
    }) |>
    dplyr::ungroup()
}

species_counts <- visitor_clean |>
  dplyr::count(guild, species, name = "n_unique_cells") |>
  dplyr::arrange(guild, species)

readr::write_csv(species_counts, file.path(OCC_OUT_DIR, "visitor_species_unique_cell_counts.csv"))

species_to_run <- species_counts |>
  dplyr::filter(n_unique_cells >= MIN_UNIQUE_CELLS_PER_SPECIES)

readr::write_csv(species_to_run, file.path(OCC_OUT_DIR, "visitor_species_to_run.csv"))

message("Species to run by guild:")
print(species_to_run |> dplyr::count(guild, name = "n_species"))

## =========================================================
## 7. ENMeval species function
## =========================================================

run_enmeval_one_species <- function(sp, guild, occ_xy, env_full, categoricals) {
  sp_safe <- safe_name(sp)
  guild_safe <- safe_name(guild)
  sp_dir <- file.path(SPECIES_OUT_DIR, guild_safe, sp_safe)
  dir.create(sp_dir, showWarnings = FALSE, recursive = TRUE)

  pred_file <- file.path(sp_dir, paste0(sp_safe, "_best_cloglog_full_template.tif"))
  pred_local_file <- file.path(sp_dir, paste0(sp_safe, "_best_cloglog_local_M.tif"))
  bin_file  <- file.path(sp_dir, paste0(sp_safe, "_binary_10p_full_template.tif"))
  rds_file  <- file.path(sp_dir, paste0(sp_safe, "_enmeval.rds"))
  result_file <- file.path(sp_dir, paste0(sp_safe, "_results.csv"))
  meta_file <- file.path(sp_dir, paste0(sp_safe, "_best_meta.csv"))
  log_file <- file.path(LOG_DIR, paste0("log_", guild_safe, "__", sp_safe, ".txt"))

  log_msg <- function(...) {
    msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
    message(msg)
    cat(msg, "\n", file = log_file, append = TRUE)
  }

  occ_xy <- make_xy_df(occ_xy$lon, occ_xy$lat)

  if (file.exists(pred_file) && file.exists(meta_file) && !FORCE_RERUN_SPECIES) {
    log_msg("SKIP existing prediction: ", sp)
    return(tibble::tibble(
      guild = guild,
      species = sp,
      status = "skipped_existing",
      n_occ = nrow(occ_xy),
      pred_file = pred_file,
      pred_local_file = pred_local_file,
      bin_file = bin_file,
      meta_file = meta_file
    ))
  }

  log_msg("==== Running ENMevaluate: ", sp, " | guild=", guild, " ====")
  log_msg("Effective occurrences before M-area cleaning: ", nrow(occ_xy))

  if (nrow(occ_xy) < 5) {
    log_msg("SKIP too few occurrences < 5")
    return(tibble::tibble(guild = guild, species = sp, status = "skipped_too_few_occ", n_occ = nrow(occ_xy)))
  }

  ## Species-specific M area with adaptive buffer
  m_obj <- NULL
  env_m <- NULL
  bg_xy <- NULL
  occ_xy_m <- NULL
  used_buffer_km <- NA_real_

  for (buf_km in BUFFER_KM_CANDIDATES) {
    log_msg("Trying species-specific M buffer: ", buf_km, " km")

    m_tmp <- make_species_m_area(
      occ_xy = occ_xy,
      env_full = env_full,
      buffer_km = buf_km,
      species_safe = paste0(sp_safe, "_", buf_km, "km")
    )

    env_tmp <- m_tmp$env

    cid <- terra::cellFromXY(env_tmp, as.matrix(occ_xy[, c("lon", "lat")]))
    occ_tmp <- occ_xy[!duplicated(cid) & !is.na(cid), , drop = FALSE]

    if (nrow(occ_tmp) < 5) {
      log_msg("  -> too few occurrences inside M: ", nrow(occ_tmp))
      next
    }

    bg_tmp <- terra::spatSample(
      env_tmp,
      size = BG_SIZE_PER_SPECIES,
      method = "random",
      values = FALSE,
      xy = TRUE,
      na.rm = TRUE
    ) |>
      as.data.frame()

    colnames(bg_tmp) <- c("lon", "lat")

    log_msg("  -> n_occ inside M: ", nrow(occ_tmp), "; n_bg: ", nrow(bg_tmp))

    if (nrow(bg_tmp) >= MIN_BG_POINTS) {
      m_obj <- m_tmp
      env_m <- env_tmp
      bg_xy <- bg_tmp
      occ_xy_m <- occ_tmp
      used_buffer_km <- buf_km
      break
    }
  }

  if (is.null(m_obj)) {
    log_msg("SKIP: no buffer size produced enough background/occurrence points")
    return(tibble::tibble(
      guild = guild,
      species = sp,
      status = "skipped_no_valid_M",
      n_occ = nrow(occ_xy),
      buffer_km = NA_real_
    ))
  }

  occ_xy <- occ_xy_m

  log_msg("Selected M buffer: ", used_buffer_km, " km")
  log_msg("Effective occurrences inside M before predictor-NA screening: ", nrow(occ_xy))
  log_msg("Background points inside M before predictor-NA screening: ", nrow(bg_xy))

  screen_obj <- screen_env_layers_for_species(
    env_m = env_m,
    occ_xy = occ_xy,
    sp_dir = sp_dir,
    sp_safe = sp_safe,
    log_msg = log_msg
  )

  env_m <- screen_obj$env
  occ_xy <- screen_obj$occ

  if (nrow(occ_xy) < MIN_OCC_AFTER_NA_FILTER) {
    log_msg("SKIP too few occurrences after predictor-NA filtering: ", nrow(occ_xy))
    return(tibble::tibble(
      guild = guild,
      species = sp,
      status = "skipped_too_few_occ_after_NA_filter",
      n_occ = nrow(occ_xy),
      buffer_km = used_buffer_km,
      dropped_layers = paste(screen_obj$dropped_layers, collapse = ";")
    ))
  }

  bg_xy <- terra::spatSample(
    env_m,
    size = BG_SIZE_PER_SPECIES,
    method = "random",
    values = FALSE,
    xy = TRUE,
    na.rm = TRUE
  ) |>
    as.data.frame()

  colnames(bg_xy) <- c("lon", "lat")

  if (nrow(bg_xy) < MIN_BG_POINTS) {
    log_msg("SKIP too few background points after predictor screening: ", nrow(bg_xy))
    return(tibble::tibble(
      guild = guild,
      species = sp,
      status = "skipped_too_few_bg_after_NA_filter",
      n_occ = nrow(occ_xy),
      n_bg = nrow(bg_xy),
      buffer_km = used_buffer_km
    ))
  }

  log_msg("Effective occurrences after predictor-NA screening: ", nrow(occ_xy))
  log_msg("Background points after predictor-NA screening: ", nrow(bg_xy))

  categoricals_m <- if (is.null(categoricals)) NULL else intersect(categoricals, names(env_m))

  ## ENMevaluate
  e <- NULL

  if (USE_BLOCK_FIRST && nrow(occ_xy) >= 20) {
    log_msg("Trying partitions = block")
    e <- tryCatch(
      ENMeval::ENMevaluate(
        occs = occ_xy,
        envs = env_m,
        bg = bg_xy,
        algorithm = "maxnet",
        partitions = "block",
        tune.args = list(fc = TUNE_FC, rm = TUNE_RM),
        categoricals = categoricals_m,
        parallel = ENMEVAL_PARALLEL,
        numCores = ENMEVAL_CORES
      ),
      error = function(err) {
        log_msg("block FAILED: ", conditionMessage(err))
        NULL
      }
    )
  }

  if (is.null(e)) {
    k <- min(RANDOM_KFOLDS, nrow(occ_xy))
    if (k >= 2) {
      log_msg("Trying partitions = randomkfold, k=", k)
      e <- tryCatch(
        ENMeval::ENMevaluate(
          occs = occ_xy,
          envs = env_m,
          bg = bg_xy,
          algorithm = "maxnet",
          partitions = "randomkfold",
          partition.settings = list(kfolds = k),
          tune.args = list(fc = TUNE_FC, rm = TUNE_RM),
          categoricals = categoricals_m,
          parallel = ENMEVAL_PARALLEL,
          numCores = ENMEVAL_CORES
        ),
        error = function(err) {
          log_msg("randomkfold FAILED: ", conditionMessage(err))
          NULL
        }
      )
    }
  }

  if (is.null(e)) {
    log_msg("ENMevaluate FAILED all attempts")
    return(tibble::tibble(
      guild = guild,
      species = sp,
      status = "enmeval_failed",
      n_occ = nrow(occ_xy),
      n_bg = nrow(bg_xy),
      buffer_km = used_buffer_km,
      dropped_layers = paste(screen_obj$dropped_layers, collapse = ";"),
      M_file = m_obj$file
    ))
  }

  res <- ENMeval::eval.results(e)
  readr::write_csv(res, result_file)
  saveRDS(e, rds_file)

  best_i <- select_best(res)
  model <- ENMeval::eval.models(e)[[best_i]]

  log_msg("Best row: ", best_i)
  if ("fc" %in% names(res) && "rm" %in% names(res)) {
    log_msg("Best fc=", res$fc[best_i], " rm=", res$rm[best_i])
  }

  pred_local <- predict_maxnet_rast(env_m, model, categoricals_m)
  terra::writeRaster(pred_local, pred_local_file, overwrite = TRUE)

  pred_full <- extend_to_template(pred_local, env_full)
  terra::writeRaster(pred_full, pred_file, overwrite = TRUE)

  occ_vals <- safe_extract_first(pred_local, occ_xy)
  occ_vals <- occ_vals[is.finite(occ_vals)]
  th10 <- if (length(occ_vals) > 0) {
    as.numeric(stats::quantile(occ_vals, probs = 0.10, na.rm = TRUE))
  } else {
    NA_real_
  }

  if (is.finite(th10)) {
    bin_local <- pred_local >= th10
    bin_full <- extend_to_template(bin_local, env_full)
    terra::writeRaster(bin_full, bin_file, overwrite = TRUE)
  } else {
    bin_file <- NA_character_
  }

  png(
    filename = file.path(sp_dir, paste0(sp_safe, "_response_curves.png")),
    width = 2000, height = 1600, res = 200
  )
  try(plot(model, type = "cloglog"), silent = TRUE)
  dev.off()

  best_meta <- dplyr::bind_cols(
    tibble::tibble(
      guild = guild,
      species = sp,
      best_row = best_i,
      threshold_10p = th10,
      n_occ = nrow(occ_xy),
      n_bg = nrow(bg_xy),
      buffer_km = used_buffer_km,
      dropped_layers = paste(screen_obj$dropped_layers, collapse = ";"),
      M_file = m_obj$file,
      pred_local_file = pred_local_file
    ),
    res[best_i, , drop = FALSE]
  )

  readr::write_csv(best_meta, meta_file)

  log_msg("SAVED full-template prediction: ", pred_file)

  tibble::tibble(
    guild = guild,
    species = sp,
    status = "success",
    n_occ = nrow(occ_xy),
    n_bg = nrow(bg_xy),
    best_row = best_i,
    threshold_10p = th10,
    buffer_km = used_buffer_km,
    dropped_layers = paste(screen_obj$dropped_layers, collapse = ";"),
    M_file = m_obj$file,
    pred_file = pred_file,
    pred_local_file = pred_local_file,
    bin_file = bin_file,
    rds_file = rds_file,
    meta_file = meta_file
  )
}

## =========================================================
## 8. ENMeval species loop
## =========================================================

results <- list()

for (i in seq_len(nrow(species_to_run))) {
  guild_i <- species_to_run$guild[i]
  sp_i <- species_to_run$species[i]

  tmp_i <- visitor_clean |>
    dplyr::filter(guild == guild_i, species == sp_i)

  occ_i <- make_xy_df(tmp_i$longitude, tmp_i$latitude)

  res_i <- tryCatch(
    run_enmeval_one_species(
      sp = sp_i,
      guild = guild_i,
      occ_xy = occ_i,
      env_full = env,
      categoricals = categoricals
    ),
    error = function(e) {
      message("FAILED species ", sp_i, ": ", conditionMessage(e))
      tibble::tibble(
        guild = guild_i,
        species = sp_i,
        status = "script_error",
        error = conditionMessage(e),
        n_occ = nrow(occ_i)
      )
    }
  )

  results[[paste(guild_i, sp_i, sep = "__")]] <- res_i

  readr::write_csv(
    dplyr::bind_rows(results),
    file.path(ROOT_DIR, "species_enmeval_run_summary_incremental.csv")
  )

  gc(verbose = FALSE)
}

species_summary <- dplyr::bind_rows(results)
readr::write_csv(species_summary, file.path(ROOT_DIR, "species_enmeval_run_summary.csv"))

## =========================================================
## 9. Build guild stacks
## =========================================================

make_group_stack <- function(guild, pred_files, bin_files) {
  guild_safe <- safe_name(guild)
  guild_dir <- file.path(GROUP_OUT_DIR, guild_safe)
  dir.create(guild_dir, showWarnings = FALSE, recursive = TRUE)

  pred_files <- pred_files[file.exists(pred_files)]
  bin_files <- bin_files[file.exists(bin_files)]

  if (length(pred_files) == 0) {
    warning("No prediction files for guild: ", guild)
    return(NULL)
  }

  message("\n==== Building group stack: ", guild, " ====")
  message("n species predictions: ", length(pred_files))

  preds <- terra::rast(pred_files)
  names(preds) <- tools::file_path_sans_ext(basename(pred_files)) |> safe_name()

  sum_r <- terra::app(preds, fun = sum, na.rm = TRUE)
  mean_r <- terra::app(preds, fun = mean, na.rm = TRUE)

  sum_file <- file.path(guild_dir, paste0("guild_", guild_safe, "_sum_cloglog_speciesM_no_swe.tif"))
  mean_file <- file.path(guild_dir, paste0("guild_", guild_safe, "_mean_cloglog_speciesM_no_swe.tif"))

  terra::writeRaster(sum_r, sum_file, overwrite = TRUE)
  terra::writeRaster(mean_r, mean_file, overwrite = TRUE)

  richness_file <- NA_character_

  if (length(bin_files) > 0) {
    bins <- terra::rast(bin_files)
    richness <- terra::app(bins, fun = sum, na.rm = TRUE)
    richness_file <- file.path(guild_dir, paste0("guild_", guild_safe, "_binary_richness_10p_speciesM_no_swe.tif"))
    terra::writeRaster(richness, richness_file, overwrite = TRUE)
  }

  tibble::tibble(
    guild = guild,
    n_species_predictions = length(pred_files),
    sum_cloglog_file = sum_file,
    mean_cloglog_file = mean_file,
    binary_richness_10p_file = richness_file
  )
}

success_tbl <- species_summary |>
  dplyr::filter(status %in% c("success", "skipped_existing")) |>
  dplyr::filter(!is.na(pred_file), file.exists(pred_file))

guilds <- sort(unique(success_tbl$guild))
group_results <- list()

for (g in guilds) {
  sub <- success_tbl |> dplyr::filter(guild == g)
  group_results[[g]] <- make_group_stack(
    guild = g,
    pred_files = sub$pred_file,
    bin_files = sub$bin_file
  )
}

group_summary <- dplyr::bind_rows(group_results)
readr::write_csv(group_summary, file.path(GROUP_OUT_DIR, "guild_stack_summary_no_swe.csv"))

## =========================================================
## 10. Optional: extract guild stacks to Cirsium points/species
## =========================================================

if (file.exists(CIRSIUM_FILE) && nrow(group_summary) > 0) {
  message("\nExtracting guild stacks to Cirsium occurrence points...")

  cir <- readr::read_csv(CIRSIUM_FILE, show_col_types = FALSE)

  if (all(c("decimalLongitude", "decimalLatitude") %in% names(cir))) {
    guild_rasters <- c(
      group_summary$sum_cloglog_file,
      group_summary$mean_cloglog_file,
      group_summary$binary_richness_10p_file
    )
    guild_rasters <- guild_rasters[!is.na(guild_rasters) & file.exists(guild_rasters)]

    if (length(guild_rasters) > 0) {
      gr <- terra::rast(guild_rasters)
      names(gr) <- tools::file_path_sans_ext(basename(guild_rasters)) |> safe_name()

      cir_xy <- cbind(
        to_numeric_coord(cir$decimalLongitude),
        to_numeric_coord(cir$decimalLatitude)
      )
      colnames(cir_xy) <- c("lon", "lat")

      vals <- safe_extract_df(gr, cir_xy)

      cir_out <- dplyr::bind_cols(cir, vals)
      readr::write_csv(cir_out, file.path(ROOT_DIR, "cirsium_occurrences_with_guild_sdm_speciesM_no_swe.csv"))

      if (all(c("species_final", "head_orientation_binary") %in% names(cir_out))) {
        sdm_cols <- names(vals)

        cir_species <- cir_out |>
          dplyr::group_by(species_final, japanese_name, head_orientation_binary) |>
          dplyr::summarise(
            n_occ = dplyr::n(),
            dplyr::across(dplyr::all_of(sdm_cols), ~ mean(.x, na.rm = TRUE)),
            .groups = "drop"
          )

        readr::write_csv(cir_species, file.path(ROOT_DIR, "cirsium_species_mean_guild_sdm_speciesM_no_swe.csv"))
      }
    }
  }
}

message("\n🎉 ALL DONE — reused existing CHELSA files, rebuilt no-SWE VIF stack + species-buffer ENMeval SDMs + guild stacks!")
message("Root output: ", ROOT_DIR)
message("VIF stack: ", ENV_VIF_TIF)
message("Species summary: ", file.path(ROOT_DIR, "species_enmeval_run_summary.csv"))
message("Guild stack summary: ", file.path(GROUP_OUT_DIR, "guild_stack_summary_no_swe.csv"))
