############################################################
## Pollinator SSDM pipeline for Cirsium head-orientation study
## FIXED VERSION
##
## Fixes:
##   - geodata::gadm() already returns a SpatVector, so do NOT use vect(jpn)
##   - added tidyr to packages because tidyr::drop_na() is used
##   - safer VIF sampling
##
## Input:
##   pollinator_gbif_outputs/gbif_pollinator_occurrences_20plus_unique_coords.csv
##
## Outputs:
##   env_for_ssdm/
##     env_stack_bio1_19_elev_slope_roughness_japan.tif
##     env_stack_vif10_japan.tif
##     env_vifstep10_selected_variables.csv
##
##   pollinator_ssdm_outputs/
##     species_sdm/*.tif
##     ssdm_by_guild/ssdm_mean_<guild>.tif
##     ssdm_by_guild/ssdm_sum_<guild>.tif
##     ssdm_by_guild/ssdm_binary_richness_<guild>.tif
##     sdm_model_summary.csv
############################################################

## =========================
## 0. Packages
## =========================

pkgs <- c(
  "terra",
  "geodata",
  "dplyr",
  "readr",
  "stringr",
  "purrr",
  "tibble",
  "tidyr",
  "usdm",
  "maxnet"
)

for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
}

library(terra)
library(geodata)
library(dplyr)
library(readr)
library(stringr)
library(purrr)
library(tibble)
library(tidyr)
library(usdm)
library(maxnet)

set.seed(42)

## =========================
## 1. Settings
## =========================

OCC_FILE <- "pollinator_gbif_outputs/gbif_pollinator_occurrences_20plus_unique_coords.csv"

ENV_DIR <- "env_for_ssdm"
OUT_DIR <- "pollinator_ssdm_outputs"
SPECIES_SDM_DIR <- file.path(OUT_DIR, "species_sdm")
GUILD_SSDM_DIR <- file.path(OUT_DIR, "ssdm_by_guild")

dir.create(ENV_DIR, showWarnings = FALSE)
dir.create(OUT_DIR, showWarnings = FALSE)
dir.create(SPECIES_SDM_DIR, showWarnings = FALSE)
dir.create(GUILD_SSDM_DIR, showWarnings = FALSE)

## WorldClim resolution:
## 0.5 arc-min ≈ 1 km, good for final analysis but heavier
## 2.5 arc-min ≈ 4–5 km, good for test run
WC_RES <- 2.5

## Japan extent with buffer
JAPAN_EXT <- ext(122, 146.5, 24, 46)

## VIF threshold
VIF_THRESHOLD <- 10

## Maxnet settings
MIN_UNIQUE_COORDS <- 20
BACKGROUND_N <- 10000
BG_MULTIPLIER <- 10
MAX_BG_PER_SPECIES <- 10000
MAXNET_FEATURES <- "lq"   # safer for small n

## Binary richness threshold.
THRESHOLD_QUANTILE <- 0.10

## Main guilds
TARGET_GUILDS <- c("butterfly", "bee", "hoverfly", "hawkmoth")

## =========================
## 2. Helper functions
## =========================

safe_filename <- function(x) {
  x %>%
    str_replace_all("[^A-Za-z0-9_\\-]+", "_") %>%
    str_replace_all("_+", "_") %>%
    str_replace_all("^_|_$", "")
}

make_pa_data <- function(pres_df, env_stack, bg_n) {
  pres_pts <- vect(
    pres_df,
    geom = c("decimalLongitude", "decimalLatitude"),
    crs = "EPSG:4326"
  )

  pres_env <- terra::extract(env_stack, pres_pts, ID = FALSE) %>%
    as_tibble()

  pres_dat <- bind_cols(
    tibble(pa = 1),
    pres_env
  ) %>%
    filter(if_all(everything(), ~ !is.na(.x)))

  bg_pts <- terra::spatSample(
    env_stack[[1]],
    size = bg_n,
    method = "random",
    na.rm = TRUE,
    as.points = TRUE,
    values = FALSE
  )

  bg_env <- terra::extract(env_stack, bg_pts, ID = FALSE) %>%
    as_tibble()

  bg_dat <- bind_cols(
    tibble(pa = 0),
    bg_env
  ) %>%
    filter(if_all(everything(), ~ !is.na(.x)))

  bind_rows(pres_dat, bg_dat)
}

fit_predict_one_species <- function(sp, guild, occ_sp, env_stack) {
  message("\n============================")
  message("Species: ", sp)
  message("Guild: ", guild)
  message("Unique coords: ", nrow(occ_sp))

  if (nrow(occ_sp) < MIN_UNIQUE_COORDS) {
    return(tibble(
      query_name = sp,
      pollinator_group = guild,
      n_unique_coords = nrow(occ_sp),
      status = "skipped_too_few_points",
      n_presence_used = NA_integer_,
      n_background_used = NA_integer_,
      sdm_file = NA_character_,
      threshold_10p = NA_real_
    ))
  }

  bg_n <- min(MAX_BG_PER_SPECIES, max(BACKGROUND_N, nrow(occ_sp) * BG_MULTIPLIER))

  dat <- make_pa_data(occ_sp, env_stack, bg_n)

  n_pres_used <- sum(dat$pa == 1)
  n_bg_used <- sum(dat$pa == 0)

  message("Presence used after env extraction: ", n_pres_used)
  message("Background used: ", n_bg_used)

  if (n_pres_used < MIN_UNIQUE_COORDS || n_bg_used < 100) {
    return(tibble(
      query_name = sp,
      pollinator_group = guild,
      n_unique_coords = nrow(occ_sp),
      status = "skipped_after_env_extraction",
      n_presence_used = n_pres_used,
      n_background_used = n_bg_used,
      sdm_file = NA_character_,
      threshold_10p = NA_real_
    ))
  }

  x <- dat %>% select(-pa) %>% as.data.frame()
  p <- dat$pa

  f <- maxnet::maxnet.formula(
    p = p,
    data = x,
    classes = MAXNET_FEATURES
  )

  mod <- tryCatch(
    maxnet::maxnet(p = p, data = x, f = f),
    error = function(e) {
      message("maxnet error: ", conditionMessage(e))
      NULL
    }
  )

  if (is.null(mod)) {
    return(tibble(
      query_name = sp,
      pollinator_group = guild,
      n_unique_coords = nrow(occ_sp),
      status = "maxnet_error",
      n_presence_used = n_pres_used,
      n_background_used = n_bg_used,
      sdm_file = NA_character_,
      threshold_10p = NA_real_
    ))
  }

  pred <- tryCatch(
    terra::predict(
      env_stack,
      mod,
      type = "cloglog",
      na.rm = TRUE
    ),
    error = function(e) {
      message("predict error: ", conditionMessage(e))
      NULL
    }
  )

  if (is.null(pred)) {
    return(tibble(
      query_name = sp,
      pollinator_group = guild,
      n_unique_coords = nrow(occ_sp),
      status = "predict_error",
      n_presence_used = n_pres_used,
      n_background_used = n_bg_used,
      sdm_file = NA_character_,
      threshold_10p = NA_real_
    ))
  }

  names(pred) <- safe_filename(sp)

  pres_x <- x[p == 1, , drop = FALSE]
  pres_pred <- tryCatch(
    as.numeric(predict(mod, pres_x, type = "cloglog")),
    error = function(e) rep(NA_real_, nrow(pres_x))
  )
  thr <- suppressWarnings(as.numeric(quantile(pres_pred, probs = THRESHOLD_QUANTILE, na.rm = TRUE)))

  out_file <- file.path(
    SPECIES_SDM_DIR,
    paste0(safe_filename(guild), "__", safe_filename(sp), ".tif")
  )

  writeRaster(pred, out_file, overwrite = TRUE)

  tibble(
    query_name = sp,
    pollinator_group = guild,
    n_unique_coords = nrow(occ_sp),
    status = "ok",
    n_presence_used = n_pres_used,
    n_background_used = n_bg_used,
    sdm_file = out_file,
    threshold_10p = thr
  )
}

## =========================
## 3. Download/create environmental rasters
## =========================

env_vif_file <- file.path(ENV_DIR, "env_stack_vif10_japan.tif")

if (!file.exists(env_vif_file)) {

  message("\nDownloading WorldClim BIO1-19...")
  bio <- geodata::worldclim_global(
    var = "bio",
    res = WC_RES,
    path = ENV_DIR
  )

  names(bio) <- paste0("BIO", 1:nlyr(bio))

  message("\nDownloading elevation...")
  elev <- geodata::elevation_global(
    res = WC_RES,
    path = ENV_DIR
  )
  names(elev) <- "elevation"

  bio_jp <- crop(bio, JAPAN_EXT)
  elev_jp <- crop(elev, JAPAN_EXT)

  elev_jp <- resample(elev_jp, bio_jp[[1]], method = "bilinear")

  slope <- terrain(elev_jp, v = "slope", unit = "degrees")
  names(slope) <- "slope"

  roughness <- terrain(elev_jp, v = "roughness")
  names(roughness) <- "roughness"

  env_all <- c(bio_jp, elev_jp, slope, roughness)

  ## -------------------------
  ## FIXED: GADM already returns SpatVector
  ## -------------------------
  message("\nDownloading Japan polygon for masking...")
  jpn <- geodata::gadm(country = "JPN", level = 0, path = ENV_DIR)

  ## geodata::gadm() returns a SpatVector, so vect(jpn) is NOT needed.
  jpn <- project(jpn, crs(env_all))

  env_all <- crop(env_all, jpn)
  env_all <- mask(env_all, jpn)

  all_env_file <- file.path(ENV_DIR, "env_stack_bio1_19_elev_slope_roughness_japan.tif")
  writeRaster(env_all, all_env_file, overwrite = TRUE)

  message("\nSaved full env stack:")
  message(all_env_file)

  ## =========================
  ## 4. VIFSTEP 10
  ## =========================

  message("\nRunning VIFSTEP threshold = ", VIF_THRESHOLD, "...")

  vif_sample <- terra::spatSample(
    env_all,
    size = 10000,
    method = "random",
    na.rm = TRUE,
    as.points = FALSE,
    values = TRUE
  ) %>%
    as.data.frame() %>%
    dplyr::select(-any_of("ID")) %>%
    tidyr::drop_na()

  vif_sample <- vif_sample[complete.cases(vif_sample), ]

  vif_res <- usdm::vifstep(vif_sample, th = VIF_THRESHOLD)

  selected_vars <- vif_res@results$Variables

  selected_tbl <- tibble(
    selected_variables = selected_vars
  )

  write_csv(
    selected_tbl,
    file.path(ENV_DIR, "env_vifstep10_selected_variables.csv")
  )

  tryCatch({
    write_csv(
      as_tibble(vif_res@results),
      file.path(ENV_DIR, "env_vifstep10_results.csv")
    )
  }, error = function(e) {})

  message("\nSelected variables:")
  print(selected_vars)

  env_vif <- env_all[[selected_vars]]

  writeRaster(env_vif, env_vif_file, overwrite = TRUE)

} else {
  message("\nVIF-selected env stack already exists. Loading:")
  message(env_vif_file)
}

env_stack <- rast(env_vif_file)

message("\nFinal env stack variables:")
print(names(env_stack))

## =========================
## 5. Read pollinator occurrences
## =========================

if (!file.exists(OCC_FILE)) {
  stop("Occurrence file not found: ", OCC_FILE)
}

occ <- read_csv(OCC_FILE, show_col_types = FALSE) %>%
  filter(pollinator_group %in% TARGET_GUILDS) %>%
  filter(!is.na(decimalLongitude), !is.na(decimalLatitude)) %>%
  distinct(query_name, pollinator_group, decimalLongitude, decimalLatitude, .keep_all = TRUE)

message("\nOccurrence rows for SSDM:")
print(nrow(occ))
message("\nSpecies by guild:")
print(occ %>% distinct(query_name, pollinator_group) %>% count(pollinator_group, name = "n_species"))

## =========================
## 6. Fit SDMs per species
## =========================

species_list <- occ %>%
  count(query_name, pollinator_group, name = "n_unique_coords") %>%
  filter(n_unique_coords >= MIN_UNIQUE_COORDS) %>%
  arrange(pollinator_group, query_name)

message("\nSpecies models to fit:")
print(species_list %>% count(pollinator_group, name = "n_species"))

model_summary <- pmap_dfr(
  list(
    sp = species_list$query_name,
    guild = species_list$pollinator_group
  ),
  function(sp, guild) {
    occ_sp <- occ %>%
      filter(query_name == sp, pollinator_group == guild)

    fit_predict_one_species(sp, guild, occ_sp, env_stack)
  }
)

write_csv(
  model_summary,
  file.path(OUT_DIR, "sdm_model_summary.csv")
)

message("\nModel status summary:")
print(model_summary %>% count(pollinator_group, status, name = "n_models"))

## =========================
## 7. Make SSDM by guild
## =========================

ok_models <- model_summary %>%
  filter(status == "ok", !is.na(sdm_file), file.exists(sdm_file))

if (nrow(ok_models) == 0) {
  stop("No successful SDM models.")
}

for (g in unique(ok_models$pollinator_group)) {

  message("\nMaking SSDM for guild: ", g)

  files_g <- ok_models %>%
    filter(pollinator_group == g) %>%
    pull(sdm_file)

  s <- rast(files_g)

  ssdm_mean <- mean(s, na.rm = TRUE)
  names(ssdm_mean) <- paste0("ssdm_mean_", g)

  mean_file <- file.path(GUILD_SSDM_DIR, paste0("ssdm_mean_", g, ".tif"))
  writeRaster(ssdm_mean, mean_file, overwrite = TRUE)

  ssdm_sum <- sum(s, na.rm = TRUE)
  names(ssdm_sum) <- paste0("ssdm_sum_", g)

  sum_file <- file.path(GUILD_SSDM_DIR, paste0("ssdm_sum_", g, ".tif"))
  writeRaster(ssdm_sum, sum_file, overwrite = TRUE)

  thresholds <- ok_models %>%
    filter(pollinator_group == g) %>%
    pull(threshold_10p)

  bin_layers <- list()

  for (i in seq_along(files_g)) {
    r <- rast(files_g[i])
    thr <- thresholds[i]

    if (is.na(thr)) next

    rb <- r >= thr
    names(rb) <- names(r)
    bin_layers[[length(bin_layers) + 1]] <- rb
  }

  if (length(bin_layers) > 0) {
    bstack <- rast(bin_layers)
    richness <- sum(bstack, na.rm = TRUE)
    names(richness) <- paste0("ssdm_binary_richness_", g)

    richness_file <- file.path(GUILD_SSDM_DIR, paste0("ssdm_binary_richness_", g, ".tif"))
    writeRaster(richness, richness_file, overwrite = TRUE)
  }
}

## =========================
## 8. Create combined SSDM stacks
## =========================

mean_files <- list.files(
  GUILD_SSDM_DIR,
  pattern = "^ssdm_mean_.*\\.tif$",
  full.names = TRUE
)

if (length(mean_files) > 0) {
  mean_stack <- rast(mean_files)
  writeRaster(
    mean_stack,
    file.path(GUILD_SSDM_DIR, "ssdm_mean_stack_all_guilds.tif"),
    overwrite = TRUE
  )
}

sum_files <- list.files(
  GUILD_SSDM_DIR,
  pattern = "^ssdm_sum_.*\\.tif$",
  full.names = TRUE
)

if (length(sum_files) > 0) {
  sum_stack <- rast(sum_files)
  writeRaster(
    sum_stack,
    file.path(GUILD_SSDM_DIR, "ssdm_sum_stack_all_guilds.tif"),
    overwrite = TRUE
  )
}

rich_files <- list.files(
  GUILD_SSDM_DIR,
  pattern = "^ssdm_binary_richness_.*\\.tif$",
  full.names = TRUE
)

if (length(rich_files) > 0) {
  rich_stack <- rast(rich_files)
  writeRaster(
    rich_stack,
    file.path(GUILD_SSDM_DIR, "ssdm_binary_richness_stack_all_guilds.tif"),
    overwrite = TRUE
  )
}

## =========================
## 9. Done
## =========================

message("\n============================")
message("DONE")
message("============================")
message("Env stack:")
message("  ", env_vif_file)
message("Model summary:")
message("  ", file.path(OUT_DIR, "sdm_model_summary.csv"))
message("Guild SSDMs:")
message("  ", GUILD_SSDM_DIR)
