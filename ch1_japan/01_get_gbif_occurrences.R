############################################################
## GBIF Cirsium occurrences
## Keep species with >=20 unique coordinate points
## Default: Japan only
############################################################

pkgs <- c("rgbif", "dplyr", "readr", "purrr", "stringr", "tibble")

for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
}

library(rgbif)
library(dplyr)
library(readr)
library(purrr)
library(stringr)
library(tibble)

## =========================
## Settings
## =========================

GENUS_NAME <- "Cirsium"
COUNTRY_CODE <- "JP"   # 日本のみ。世界全体なら NULL
MIN_UNIQUE_COORDS <- 20

HAS_COORDINATE <- TRUE
HAS_GEOSPATIAL_ISSUE <- FALSE

PAGE_LIMIT <- 300
SLEEP_SEC <- 0.3

OUT_DIR <- "gbif_cirsium_outputs"
dir.create(OUT_DIR, showWarnings = FALSE)

SPECIES_DIR <- file.path(OUT_DIR, "gbif_cirsium_20plus_unique_coords_by_species")
dir.create(SPECIES_DIR, showWarnings = FALSE)

## =========================
## Get GBIF taxonKey for Cirsium
## =========================

bb <- name_backbone(name = GENUS_NAME, rank = "GENUS")

print(bb)

if (is.null(bb$usageKey) || is.na(bb$usageKey)) {
  stop("GBIF taxonKey for Cirsium was not found.")
}

GENUS_KEY <- bb$usageKey
cat("\nCirsium genus taxonKey:", GENUS_KEY, "\n")

## =========================
## Count occurrences
## =========================

count_args <- list(
  taxonKey = GENUS_KEY,
  hasCoordinate = HAS_COORDINATE,
  hasGeospatialIssue = HAS_GEOSPATIAL_ISSUE
)

if (!is.null(COUNTRY_CODE)) {
  count_args$country <- COUNTRY_CODE
}

total_n <- do.call(occ_count, count_args)

cat("\nTotal GBIF records:", total_n, "\n")

if (total_n == 0) stop("No records found.")

if (total_n > 100000) {
  stop("Too many records for occ_search paging. Use GBIF download API.")
}

## =========================
## Download all pages
## =========================

offsets <- seq(0, total_n - 1, by = PAGE_LIMIT)

download_one_page <- function(start_offset) {
  cat("Downloading offset:", start_offset, "\n")
  
  search_args <- list(
    taxonKey = GENUS_KEY,
    hasCoordinate = HAS_COORDINATE,
    hasGeospatialIssue = HAS_GEOSPATIAL_ISSUE,
    limit = PAGE_LIMIT,
    start = start_offset
  )
  
  if (!is.null(COUNTRY_CODE)) {
    search_args$country <- COUNTRY_CODE
  }
  
  res <- do.call(occ_search, search_args)
  Sys.sleep(SLEEP_SEC)
  
  if (is.null(res$data) || nrow(res$data) == 0) {
    return(tibble())
  }
  
  as_tibble(res$data)
}

occ_all <- map_dfr(offsets, download_one_page)

cat("\nDownloaded rows:", nrow(occ_all), "\n")

## =========================
## Clean columns
## =========================

pick_col <- function(df, colname) {
  if (colname %in% names(df)) df[[colname]] else NA
}

occ_clean <- occ_all %>%
  mutate(
    gbifID = as.character(pick_col(., "key")),
    speciesKey = as.character(pick_col(., "speciesKey")),
    taxonKey = as.character(pick_col(., "taxonKey")),
    scientificName = as.character(pick_col(., "scientificName")),
    acceptedScientificName = as.character(pick_col(., "acceptedScientificName")),
    species = as.character(pick_col(., "species")),
    genus = as.character(pick_col(., "genus")),
    decimalLatitude = as.numeric(pick_col(., "decimalLatitude")),
    decimalLongitude = as.numeric(pick_col(., "decimalLongitude")),
    countryCode = as.character(pick_col(., "countryCode")),
    eventDate = as.character(pick_col(., "eventDate")),
    year = suppressWarnings(as.integer(pick_col(., "year"))),
    basisOfRecord = as.character(pick_col(., "basisOfRecord")),
    datasetKey = as.character(pick_col(., "datasetKey")),
    institutionCode = as.character(pick_col(., "institutionCode")),
    collectionCode = as.character(pick_col(., "collectionCode")),
    catalogNumber = as.character(pick_col(., "catalogNumber")),
    occurrenceStatus = as.character(pick_col(., "occurrenceStatus"))
  ) %>%
  mutate(
    species_from_name = str_extract(scientificName, "^Cirsium\\s+[A-Za-z-]+"),
    species_final = case_when(
      !is.na(species) & species != "" ~ species,
      !is.na(acceptedScientificName) & acceptedScientificName != "" ~
        str_extract(acceptedScientificName, "^Cirsium\\s+[A-Za-z-]+"),
      !is.na(species_from_name) ~ species_from_name,
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(decimalLatitude), !is.na(decimalLongitude)) %>%
  filter(!is.na(species_final)) %>%
  filter(str_detect(species_final, "^Cirsium\\s+")) %>%
  mutate(
    ## 座標重複判定用
    coord_key = paste(decimalLatitude, decimalLongitude, sep = "_")
  )

write_csv(
  occ_clean,
  file.path(OUT_DIR, "gbif_cirsium_all_japan_raw_records.csv")
)

cat("\nSaved raw cleaned records:\n")
cat(file.path(OUT_DIR, "gbif_cirsium_all_japan_raw_records.csv"), "\n")

## =========================
## Remove duplicate coordinates within each species
## =========================

occ_unique_coords <- occ_clean %>%
  group_by(species_final, decimalLatitude, decimalLongitude) %>%
  arrange(
    desc(!is.na(year)),
    desc(year),
    .by_group = TRUE
  ) %>%
  slice(1) %>%
  ungroup()

write_csv(
  occ_unique_coords,
  file.path(OUT_DIR, "gbif_cirsium_all_japan_unique_coords.csv")
)

cat("\nSaved unique-coordinate records:\n")
cat(file.path(OUT_DIR, "gbif_cirsium_all_japan_unique_coords.csv"), "\n")

## =========================
## Count species by raw records and unique coordinates
## =========================

species_counts <- occ_clean %>%
  count(species_final, name = "n_raw_records") %>%
  full_join(
    occ_unique_coords %>%
      count(species_final, name = "n_unique_coords"),
    by = "species_final"
  ) %>%
  mutate(
    n_raw_records = ifelse(is.na(n_raw_records), 0, n_raw_records),
    n_unique_coords = ifelse(is.na(n_unique_coords), 0, n_unique_coords)
  ) %>%
  arrange(desc(n_unique_coords), desc(n_raw_records), species_final)

write_csv(
  species_counts,
  file.path(OUT_DIR, "gbif_cirsium_species_counts_unique_coords_japan.csv")
)

cat("\nSpecies counts by unique coordinates:\n")
print(species_counts, n = Inf)

## =========================
## Keep species with >=20 unique coordinates
## =========================

species_20plus_unique <- species_counts %>%
  filter(n_unique_coords >= MIN_UNIQUE_COORDS)

write_csv(
  species_20plus_unique,
  file.path(OUT_DIR, "gbif_cirsium_species_20plus_unique_coords_japan.csv")
)

cat("\nSpecies with >=", MIN_UNIQUE_COORDS, "unique coordinates:\n")
print(species_20plus_unique, n = Inf)

occ_20plus_unique <- occ_unique_coords %>%
  semi_join(species_20plus_unique, by = "species_final") %>%
  arrange(species_final, decimalLatitude, decimalLongitude)

write_csv(
  occ_20plus_unique,
  file.path(OUT_DIR, "gbif_cirsium_occurrences_20plus_unique_coords_japan.csv")
)

cat("\nSaved occurrences for species with >=20 unique coordinates:\n")
cat(file.path(OUT_DIR, "gbif_cirsium_occurrences_20plus_unique_coords_japan.csv"), "\n")

## =========================
## Save each species separately
## =========================

safe_filename <- function(x) {
  x %>%
    str_replace_all("[^A-Za-z0-9_\\-]+", "_") %>%
    str_replace_all("_+", "_")
}

for (sp in species_20plus_unique$species_final) {
  dat_sp <- occ_20plus_unique %>%
    filter(species_final == sp)
  
  out_file <- file.path(
    SPECIES_DIR,
    paste0(safe_filename(sp), ".csv")
  )
  
  write_csv(dat_sp, out_file)
}

cat("\nSaved per-species unique-coordinate CSVs in:\n")
cat(SPECIES_DIR, "\n")

## =========================
## Done
## =========================

cat("\n============================\n")
cat("DONE\n")
cat("============================\n")
cat("Raw cleaned records:", nrow(occ_clean), "\n")
cat("Unique-coordinate records:", nrow(occ_unique_coords), "\n")
cat("Number of species:", nrow(species_counts), "\n")
cat("Species with >=20 unique coords:", nrow(species_20plus_unique), "\n")

cat("\nOutputs:\n")
cat("  ", file.path(OUT_DIR, "gbif_cirsium_all_japan_raw_records.csv"), "\n")
cat("  ", file.path(OUT_DIR, "gbif_cirsium_all_japan_unique_coords.csv"), "\n")
cat("  ", file.path(OUT_DIR, "gbif_cirsium_species_counts_unique_coords_japan.csv"), "\n")
cat("  ", file.path(OUT_DIR, "gbif_cirsium_species_20plus_unique_coords_japan.csv"), "\n")
cat("  ", file.path(OUT_DIR, "gbif_cirsium_occurrences_20plus_unique_coords_japan.csv"), "\n")
cat("  ", SPECIES_DIR, "\n")