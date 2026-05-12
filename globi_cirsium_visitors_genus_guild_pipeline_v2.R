############################################################
## GloBI Cirsium visitors: clean genus-guild pipeline v2
##
## 修正点:
##   - "Bombus sp.", "Bombus spp.", "cf.", "aff.", "group", "complex" などを
##     SSDMターゲットから除外
##   - GloBI の source_citation は存在しないため study_title を使う
##   - genusベースで bee / butterfly / hoverfly / hawkmoth に分類
##
## Main output:
##   globi_cirsium_outputs/globi_cirsium_pollinator_species_for_ssdm_genus_guild.csv
############################################################

pkgs <- c("httr", "jsonlite", "dplyr", "stringr", "readr", "tibble")

for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

library(httr)
library(jsonlite)
library(dplyr)
library(stringr)
library(readr)
library(tibble)

## =========================
## 1. Settings
## =========================

OUT_DIR <- "globi_cirsium_outputs"
dir.create(OUT_DIR, showWarnings = FALSE)

SOURCE_TAXON <- "Cirsium"
INTERACTION_TYPE <- "visitedBy"
GLOBI_LIMIT <- 10000

KEEP_ONLY_BINOMIAL <- TRUE

## =========================
## 2. Helper functions
## =========================

safe_col <- function(df, col) {
  if (col %in% names(df)) {
    df[[col]]
  } else {
    rep(NA_character_, nrow(df))
  }
}

clean_taxon_name <- function(x) {
  x <- as.character(x)
  x <- str_replace_all(x, "_", " ")
  x <- str_squish(x)
  x
}

is_bad_taxon_name <- function(x) {
  str_detect(
    x,
    regex(
      "\\bsp\\.?$|\\bspp\\.?$|\\bcf\\.?\\b|\\baff\\.?\\b|\\bgroup\\b|\\bcomplex\\b|\\bindet\\.?\\b|\\bunknown\\b",
      ignore_case = TRUE
    )
  )
}

is_binomial_name <- function(x) {
  str_detect(x, "^[A-Z][A-Za-z-]+\\s+[a-z][A-Za-z-]+(?:\\s+[a-z][A-Za-z-]+)?$") &
    !is_bad_taxon_name(x)
}

## =========================
## 3. GloBI API
## =========================

globi_url <- "https://api.globalbioticinteractions.org/interaction"

response <- GET(
  globi_url,
  query = list(
    sourceTaxon = SOURCE_TAXON,
    interactionType = INTERACTION_TYPE,
    limit = GLOBI_LIMIT
  )
)

stop_for_status(response)

json_txt <- content(response, as = "text", encoding = "UTF-8")
json_list <- fromJSON(json_txt)

df_visitors <- as.data.frame(json_list$data, stringsAsFactors = FALSE)
colnames(df_visitors) <- json_list$columns

write_csv(
  df_visitors,
  file.path(OUT_DIR, "globi_cirsium_visitedBy_raw.csv")
)

cat("\n=== GloBI columns ===\n")
print(names(df_visitors))
cat("\nRaw interaction rows:", nrow(df_visitors), "\n")

## =========================
## 4. Extract unique visitors
## =========================

visitors <- tibble(
  source_taxon_name = safe_col(df_visitors, "source_taxon_name"),
  interaction_type  = safe_col(df_visitors, "interaction_type"),
  visitor_name      = safe_col(df_visitors, "target_taxon_name"),
  visitor_id        = safe_col(df_visitors, "target_taxon_external_id"),
  source_reference  = safe_col(df_visitors, "study_title"),
  latitude          = safe_col(df_visitors, "latitude"),
  longitude         = safe_col(df_visitors, "longitude")
) %>%
  mutate(
    visitor_name = clean_taxon_name(visitor_name),
    visitor_genus = str_extract(visitor_name, "^[A-Za-z]+"),
    is_bad_name = is_bad_taxon_name(visitor_name),
    n_words = str_count(visitor_name, "\\s+") + 1,
    is_binomial = is_binomial_name(visitor_name)
  ) %>%
  filter(!is.na(visitor_name), visitor_name != "") %>%
  distinct(visitor_name, .keep_all = TRUE) %>%
  arrange(visitor_name)

write_csv(
  visitors,
  file.path(OUT_DIR, "globi_cirsium_visitors_unique.csv")
)

cat("\nUnique visitor taxa:", nrow(visitors), "\n")
cat("\nBinomial / non-binomial after bad-name filter:\n")
print(table(visitors$is_binomial, useNA = "ifany"))

cat("\nBad-name examples excluded from SSDM:\n")
print(
  visitors %>%
    filter(is_bad_name) %>%
    select(visitor_name, visitor_genus, is_binomial) %>%
    head(30)
)

## =========================
## 5. Genus-based guild
## =========================

bee_genus <- c(
  "Bombus", "Xylocopa", "Apis", "Megachile", "Eucera",
  "Halictus", "Lasioglossum", "Andrena", "Osmia",
  "Ceratina", "Hylaeus", "Anthidium", "Panurgus",
  "Colletes", "Anthophora", "Nomada", "Melitta"
)

butter_genus <- c(
  "Papilio", "Pieris", "Vanessa", "Argynnis", "Aglais",
  "Fabriciana", "Erynnis", "Polygonia", "Nymphalis",
  "Cynthia", "Speyeria", "Aporia", "Erebia", "Araschnia",
  "Colias", "Gonepteryx", "Limenitis", "Lycaena",
  "Lysandra", "Aricia", "Satyrium", "Thecla",
  "Aphantopus", "Coenonympha", "Pararge", "Lasiommata",
  "Melanargia", "Ochlodes", "Thymelicus", "Pyrgus",
  "Hesperia", "Boloria", "Clossiana", "Euphydryas",
  "Hipparchia", "Melitaea", "Maniola", "Polyommatus"
)

hover_genus <- c(
  "Episyrphus", "Eristalis", "Eristalinus", "Melanostoma",
  "Syrphus", "Sphaerophoria", "Eupeodes", "Cheilosia",
  "Volucella", "Platycheirus", "Rhingia", "Sericomyia",
  "Chrysotoxum", "Dasysyrphus", "Didea", "Eriozona",
  "Leucozona", "Myolepta", "Parhelophilus", "Pyrophaena",
  "Neoascia", "Xylota", "Scaeva", "Syritta",
  "Ferdinandea", "Criorhina", "Merodon", "Epistrophe",
  "Meliscaeva", "Myathropa", "Helophilus", "Chrysogaster"
)

macro_genus <- c("Macroglossum", "Hemaris")

visitors_grouped <- visitors %>%
  mutate(
    guild = case_when(
      visitor_genus %in% bee_genus    ~ "bee",
      visitor_genus %in% butter_genus ~ "butterfly",
      visitor_genus %in% hover_genus  ~ "hoverfly",
      visitor_genus %in% macro_genus  ~ "hawkmoth",
      TRUE ~ "other_or_unknown"
    )
  )

write_csv(
  visitors_grouped,
  file.path(OUT_DIR, "globi_cirsium_visitors_genus_guild.csv")
)

guild_summary <- visitors_grouped %>%
  count(guild, sort = TRUE)

write_csv(
  guild_summary,
  file.path(OUT_DIR, "globi_cirsium_visitors_genus_guild_summary.csv")
)

cat("\n=== Genus-based guild summary ===\n")
print(guild_summary)

## =========================
## 6. SSDM target list
## =========================

ssdm_targets <- visitors_grouped %>%
  filter(guild %in% c("bee", "butterfly", "hoverfly", "hawkmoth")) %>%
  filter(if (KEEP_ONLY_BINOMIAL) is_binomial else TRUE) %>%
  filter(!is_bad_name) %>%
  distinct(visitor_name, visitor_genus, guild, .keep_all = TRUE) %>%
  arrange(guild, visitor_name)

write_csv(
  ssdm_targets,
  file.path(OUT_DIR, "globi_cirsium_pollinator_species_for_ssdm_genus_guild.csv")
)

ssdm_summary <- ssdm_targets %>%
  count(guild, name = "n_taxa_for_ssdm") %>%
  arrange(desc(n_taxa_for_ssdm))

write_csv(
  ssdm_summary,
  file.path(OUT_DIR, "globi_cirsium_pollinator_species_for_ssdm_summary.csv")
)

cat("\n=== SSDM target taxa by guild ===\n")
print(ssdm_summary)

for (g in unique(ssdm_targets$guild)) {
  out_g <- ssdm_targets %>% filter(guild == g)

  write_csv(
    out_g,
    file.path(OUT_DIR, paste0("globi_cirsium_pollinator_species_for_ssdm_", g, ".csv"))
  )
}

## =========================
## 7. Check files
## =========================

unknown_check <- visitors_grouped %>%
  filter(guild == "other_or_unknown") %>%
  select(visitor_name, visitor_genus, is_binomial, is_bad_name, visitor_id, source_reference) %>%
  arrange(visitor_genus, visitor_name)

write_csv(
  unknown_check,
  file.path(OUT_DIR, "globi_cirsium_visitors_other_or_unknown_check.csv")
)

bad_name_check <- visitors_grouped %>%
  filter(is_bad_name) %>%
  select(visitor_name, visitor_genus, guild, is_binomial, is_bad_name, visitor_id, source_reference) %>%
  arrange(guild, visitor_name)

write_csv(
  bad_name_check,
  file.path(OUT_DIR, "globi_cirsium_visitors_bad_names_excluded.csv")
)

cat("\nOther/unknown taxa:", nrow(unknown_check), "\n")
cat("Bad names excluded:", nrow(bad_name_check), "\n")

## =========================
## 8. Done
## =========================

cat("\n============================\n")
cat("DONE\n")
cat("============================\n")
cat("Main SSDM target file:\n")
cat("  ", file.path(OUT_DIR, "globi_cirsium_pollinator_species_for_ssdm_genus_guild.csv"), "\n")
