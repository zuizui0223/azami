############################################################
## Cirsium head orientation
## Pollinator grouping sensitivity analysis + model comparison
##
## Purpose:
##   Compare how different pollinator grouping schemes affect models of
##   regional proportion of nodding-flowered Cirsium species.
##
## Response:
##   cbind(n_nodding_species, n_upward_species)
##   = proportion/probability of nodding-flowered species in each grid.
##
## Environmental predictors:
##   Pre-VIF candidate environmental raster stack -> grid-level PCA.
##   This matches the grid-level response.
##
## Pollinator grouping schemes compared:
##   1. broad4_direct:
##        bee + butterfly + fly + hawkmoth
##   2. lepidoptera3_direct:
##        bee + lepidoptera(=butterfly+hawkmoth) + fly
##   3. order3_direct:
##        Hymenoptera + Lepidoptera + Diptera
##        Here Hymenoptera is bee SDM set; Diptera is fly/hoverfly SDM set.
##   4. order3_pca:
##        PCA of Hymenoptera / Lepidoptera / Diptera
##   5. bombus_vs_nonbombus_all:
##        Bombus + all non-Bombus pollinators
##   6. bombus_vs_nonbombus_bees:
##        Bombus + non-Bombus bees only
##   7. total_pollinator:
##        all pollinator SDMs merged
##
## Main model families:
##   - binomial GLM
##   - beta-binomial via glmmTMB
##
## Main outputs:
##   all_scheme_binomial_GLM_model_comparison.csv
##   all_scheme_beta_binomial_model_comparison.csv
##   all_scheme_best_models.csv
##   all_scheme_coefficients.csv
##   all_scheme_predictor_correlations.csv
##   figures/model_comparison_all_schemes_beta_binomial.png
############################################################

## =========================================================
## 0. Packages
## =========================================================

pkgs <- c(
  "terra", "dplyr", "readr", "stringr", "tibble", "tidyr",
  "ggplot2", "forcats", "pROC", "glmmTMB", "purrr"
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
  library(forcats)
  library(pROC)
  library(glmmTMB)
  library(purrr)
})

set.seed(42)

## =========================================================
## 1. Settings
## =========================================================

ROOT_DIR <- "chelsa_pollinator_enmeval_rebuild_no_swe"

SPECIES_SUMMARY_FILE <- file.path(ROOT_DIR, "species_enmeval_run_summary.csv")
ENV_CANDIDATE_TIF <- file.path(ROOT_DIR, "env_candidates", "chelsa_candidates_no_swe_japan.tif")
ENV_CANDIDATE_FILELIST <- file.path(ROOT_DIR, "env_candidates", "chelsa_candidates_no_swe_filelist.csv")

## Candidate occurrence files. The first existing one is used.
OCCURRENCE_CANDIDATES <- c(
  file.path(ROOT_DIR, "cirsium_occurrences_with_guild_sdm_speciesM_no_swe.csv"),
  file.path(ROOT_DIR, "cirsium_candidateEnvPCA_pollinator_metric_reduced", "01_cirsium_occurrences_with_candidateEnv_and_guildSDM.csv"),
  file.path(ROOT_DIR, "cirsium_candidateEnvPCA_pollinator_metric", "01_cirsium_occurrences_with_candidateEnv_and_guildSDM.csv")
)

OUT_DIR <- file.path(ROOT_DIR, "cirsium_pollinator_grouping_sensitivity_model_comparison")
STACK_DIR <- file.path(OUT_DIR, "group_stacks")
FIG_DIR <- file.path(OUT_DIR, "figures")
MAP_DIR <- file.path(FIG_DIR, "maps")
COEF_DIR <- file.path(FIG_DIR, "coefficients")
PCA_DIR <- file.path(FIG_DIR, "pca_diagnostics")
TMP_DIR <- file.path(OUT_DIR, "terra_tmp")

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(STACK_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(MAP_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(COEF_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(PCA_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TMP_DIR, showWarnings = FALSE, recursive = TRUE)

terra::terraOptions(memfrac = 0.45, tempdir = TMP_DIR)

GRID_SIZE_DEG <- 0.5
MIN_TOTAL_SPECIES_PER_GRID <- 2
N_ENV_PC_USE <- 3
N_ENV_PC_SAVE_LOADINGS <- 6
N_ORDER_PC_USE <- 2

## If TRUE, rebuild broad-group raster stacks even if files exist.
FORCE_REBUILD_STACKS <- FALSE

## Correlation thresholds for diagnostics only.
COR_WARN_THRESHOLD <- 0.70
COR_DANGER_THRESHOLD <- 0.85

## =========================================================
## 2. Helper functions
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

to_numeric_coord <- function(x) {
  if (is.list(x)) x <- unlist(x)
  x <- as.character(x)
  x <- stringr::str_replace_all(x, ",", ".")
  suppressWarnings(as.numeric(x))
}

extract_genus <- function(species_name) {
  genus <- stringr::str_extract(as.character(species_name), "^[A-Za-z]+")
  ifelse(is.na(genus) | genus == "", "Unknown", genus)
}

safe_scale <- function(x) {
  x <- as.numeric(x)
  m <- mean(x, na.rm = TRUE)
  s <- stats::sd(x, na.rm = TRUE)
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

brier_prop <- function(obs_prop, p) {
  mean((obs_prop - p)^2, na.rm = TRUE)
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

replace_nonfinite_with_mean <- function(df) {
  df2 <- as.data.frame(df)
  for (nm in names(df2)) {
    x <- as.numeric(df2[[nm]])
    if (all(!is.finite(x))) {
      x[] <- 0
    } else {
      x[!is.finite(x)] <- mean(x, na.rm = TRUE)
    }
    df2[[nm]] <- x
  }
  df2
}

get_rhs <- function(vars) {
  vars <- vars[!is.na(vars) & vars != ""]
  if (length(vars) == 0) return("1")
  paste(vars, collapse = " + ")
}

make_binom_formula <- function(rhs) {
  as.formula(paste("cbind(n_nodding_species, n_upward_species) ~", rhs))
}

empty_summary <- function(scheme, model, method, error = NA_character_) {
  tibble(
    scheme = scheme,
    model = model,
    method = method,
    AIC = NA_real_,
    BIC = NA_real_,
    logLik = NA_real_,
    df = NA_real_,
    n_grid = NA_integer_,
    AUC_presence = NA_real_,
    logloss = NA_real_,
    brier = NA_real_,
    error = error
  )
}

empty_coef <- function(scheme, model, method, error = NA_character_) {
  tibble(
    scheme = scheme,
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

coef_from_glm <- function(m, scheme, model_name, method_name) {
  sm <- summary(m)$coefficients
  as.data.frame(sm) |>
    tibble::rownames_to_column("term") |>
    tibble::as_tibble() |>
    dplyr::rename(
      estimate = Estimate,
      std.error = `Std. Error`,
      statistic = `z value`,
      p.value = `Pr(>|z|)`
    ) |>
    dplyr::mutate(scheme = scheme, model = model_name, method = method_name, .before = 1)
}

coef_from_glmmTMB <- function(m, scheme, model_name, method_name) {
  sm <- summary(m)$coefficients$cond
  as.data.frame(sm) |>
    tibble::rownames_to_column("term") |>
    tibble::as_tibble() |>
    dplyr::rename(
      estimate = Estimate,
      std.error = `Std. Error`,
      statistic = `z value`,
      p.value = `Pr(>|z|)`
    ) |>
    dplyr::mutate(scheme = scheme, model = model_name, method = method_name, .before = 1)
}

fit_glm_binom <- function(formula, data, scheme, model_name) {
  method <- "binomial_GLM"
  tryCatch({
    m <- glm(formula, data = data, family = binomial)
    p <- predict(m, type = "response")
    y_presence <- ifelse(data$n_nodding_species > 0, 1, 0)
    list(
      model = m,
      summary = tibble(
        scheme = scheme,
        model = model_name,
        method = method,
        AIC = AIC(m),
        BIC = BIC(m),
        logLik = as.numeric(logLik(m)),
        df = attr(logLik(m), "df"),
        n_grid = nrow(data),
        AUC_presence = safe_auc(y_presence, p),
        logloss = logloss_binomial(data$n_nodding_species, data$n_total_species, p),
        brier = brier_prop(data$prop_nodding_species, p),
        error = NA_character_
      ),
      coef = coef_from_glm(m, scheme, model_name, method)
    )
  }, error = function(e) {
    list(
      model = NULL,
      summary = empty_summary(scheme, model_name, method, conditionMessage(e)),
      coef = empty_coef(scheme, model_name, method, conditionMessage(e))
    )
  })
}

fit_glmmTMB_betabinom <- function(formula, data, scheme, model_name) {
  method <- "beta_binomial_glmmTMB"
  tryCatch({
    m <- glmmTMB::glmmTMB(formula, data = data, family = glmmTMB::betabinomial(link = "logit"))
    p <- predict(m, type = "response")
    y_presence <- ifelse(data$n_nodding_species > 0, 1, 0)
    list(
      model = m,
      summary = tibble(
        scheme = scheme,
        model = model_name,
        method = method,
        AIC = AIC(m),
        BIC = BIC(m),
        logLik = as.numeric(logLik(m)),
        df = attr(logLik(m), "df"),
        n_grid = nrow(data),
        AUC_presence = safe_auc(y_presence, p),
        logloss = logloss_binomial(data$n_nodding_species, data$n_total_species, p),
        brier = brier_prop(data$prop_nodding_species, p),
        error = NA_character_
      ),
      coef = coef_from_glmmTMB(m, scheme, model_name, method)
    )
  }, error = function(e) {
    list(
      model = NULL,
      summary = empty_summary(scheme, model_name, method, conditionMessage(e)),
      coef = empty_coef(scheme, model_name, method, conditionMessage(e))
    )
  })
}

add_delta_aic <- function(df) {
  df |>
    dplyr::group_by(method) |>
    dplyr::mutate(
      delta_AIC = ifelse(is.finite(AIC), AIC - min(AIC, na.rm = TRUE), NA_real_),
      akaike_weight = ifelse(
        is.finite(delta_AIC),
        exp(-0.5 * delta_AIC) / sum(exp(-0.5 * delta_AIC), na.rm = TRUE),
        NA_real_
      )
    ) |>
    dplyr::ungroup()
}

cor_pairs_table <- function(dat, vars, scheme) {
  vars <- vars[vars %in% names(dat)]
  if (length(vars) < 2) return(tibble())
  cm <- cor(dat[, vars, drop = FALSE], use = "pairwise.complete.obs")
  as.data.frame(as.table(cm)) |>
    as_tibble() |>
    rename(var1 = Var1, var2 = Var2, correlation = Freq) |>
    filter(as.character(var1) < as.character(var2)) |>
    arrange(desc(abs(correlation))) |>
    mutate(
      scheme = scheme,
      abs_correlation = abs(correlation),
      flag = case_when(
        abs_correlation >= COR_DANGER_THRESHOLD ~ "danger_0.85",
        abs_correlation >= COR_WARN_THRESHOLD ~ "warning_0.70",
        TRUE ~ "ok"
      ),
      .before = 1
    )
}

plot_model_compare <- function(df, out_file, title) {
  plot_df <- df |>
    filter(is.finite(AIC)) |>
    arrange(AIC) |>
    mutate(label = paste(scheme, model, sep = " / ")) |>
    slice_head(n = 40) |>
    mutate(label = factor(label, levels = rev(label)))
  if (nrow(plot_df) == 0) return(NULL)
  p <- ggplot(plot_df, aes(x = delta_AIC, y = label)) +
    geom_col() +
    theme_bw(base_size = 11) +
    labs(x = "Delta AIC", y = NULL, title = title)
  ggsave(out_file, p, width = 11, height = 10, dpi = 300)
  p
}

plot_coef <- function(coef_df, out_file, title) {
  df <- coef_df |>
    filter(term != "(Intercept)", is.finite(estimate), is.finite(std.error)) |>
    mutate(
      conf.low = estimate - 1.96 * std.error,
      conf.high = estimate + 1.96 * std.error,
      term = forcats::fct_reorder(term, estimate)
    )
  if (nrow(df) == 0) return(NULL)
  p <- ggplot(df, aes(x = estimate, y = term)) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.15) +
    geom_point(size = 2) +
    theme_bw(base_size = 11) +
    labs(x = "Coefficient", y = NULL, title = title)
  ggsave(out_file, p, width = 9, height = 6.5, dpi = 300)
  p
}

plot_map <- function(dat, value, out_file, title, diverging = FALSE) {
  if (!value %in% names(dat)) return(NULL)
  p <- ggplot(dat, aes(x = grid_lon, y = grid_lat, color = .data[[value]], size = n_total_species)) +
    geom_point(alpha = 0.9) +
    coord_equal() +
    theme_bw(base_size = 11) +
    labs(x = "Longitude", y = "Latitude", color = value, size = "Total species", title = title)
  if (diverging) {
    p <- p + scale_color_gradient2(midpoint = 0)
  } else {
    p <- p + scale_color_viridis_c(limits = c(0, 1))
  }
  ggsave(out_file, p, width = 7.5, height = 7, dpi = 300)
  p
}

## =========================================================
## 3. Pollinator classification
## =========================================================

BEE_GENERA <- c(
  "Bombus", "Xylocopa", "Eucera", "Apis", "Megachile", "Ceratina", "Lasioglossum",
  "Andrena", "Colletes", "Nomada", "Osmia", "Halictus", "Hylaeus", "Sphecodes",
  "Seladonia", "Anthophora"
)

BUTTERFLY_GENERA <- c(
  "Pieris", "Aglais", "Aporia", "Araschnia", "Argynnis", "Fabriciana", "Lycaena",
  "Nymphalis", "Ochlodes", "Papilio", "Polygonia", "Satyrium", "Speyeria",
  "Cynthia", "Vanessa"
)

FLY_GENERA <- c(
  "Syrphus", "Eristalis", "Melanostoma", "Volucella", "Dasysyrphus", "Didea",
  "Episyrphus", "Eupeodes", "Ferdinandea", "Leucozona", "Meliscaeva", "Merodon",
  "Platycheirus", "Rhingia", "Sphaerophoria"
)

HAWKMOTH_GENERA <- c("Macroglossum")

classify_basic_group <- function(genus) {
  case_when(
    genus %in% BEE_GENERA ~ "bee",
    genus %in% BUTTERFLY_GENERA ~ "butterfly",
    genus %in% FLY_GENERA ~ "fly",
    genus %in% HAWKMOTH_GENERA ~ "hawkmoth",
    TRUE ~ "unknown"
  )
}

## Scheme definitions: each named vector gives group name for a basic group/genus logic.
## These are applied to individual SDM rows.
assign_scheme_group <- function(basic_group, genus, scheme) {
  if (scheme == "broad4_direct") {
    return(basic_group)
  }
  if (scheme == "lepidoptera3_direct") {
    return(case_when(
      basic_group == "bee" ~ "bee",
      basic_group %in% c("butterfly", "hawkmoth") ~ "lepidoptera",
      basic_group == "fly" ~ "fly",
      TRUE ~ "unknown"
    ))
  }
  if (scheme == "order3_direct" || scheme == "order3_pca") {
    return(case_when(
      basic_group == "bee" ~ "hymenoptera",
      basic_group %in% c("butterfly", "hawkmoth") ~ "lepidoptera",
      basic_group == "fly" ~ "diptera",
      TRUE ~ "unknown"
    ))
  }
  if (scheme == "bombus_vs_nonbombus_all") {
    return(case_when(
      genus == "Bombus" ~ "bombus",
      basic_group %in% c("bee", "butterfly", "hawkmoth", "fly") ~ "non_bombus",
      TRUE ~ "unknown"
    ))
  }
  if (scheme == "bombus_vs_nonbombus_bees") {
    return(case_when(
      genus == "Bombus" ~ "bombus",
      basic_group == "bee" & genus != "Bombus" ~ "non_bombus_bees",
      TRUE ~ "exclude"
    ))
  }
  if (scheme == "total_pollinator") {
    return(case_when(
      basic_group %in% c("bee", "butterfly", "hawkmoth", "fly") ~ "total_pollinator",
      TRUE ~ "unknown"
    ))
  }
  stop("Unknown scheme: ", scheme)
}

SCHEMES <- c(
  "broad4_direct",
  "lepidoptera3_direct",
  "order3_direct",
  "order3_pca",
  "bombus_vs_nonbombus_all",
  "bombus_vs_nonbombus_bees",
  "total_pollinator"
)

## =========================================================
## 4. Load inputs
## =========================================================

stop_if_missing(SPECIES_SUMMARY_FILE)

OCCURRENCE_FILE <- OCCURRENCE_CANDIDATES[file.exists(OCCURRENCE_CANDIDATES)][1]
if (is.na(OCCURRENCE_FILE) || length(OCCURRENCE_FILE) == 0) {
  stop("No occurrence file found. Checked:\n", paste(OCCURRENCE_CANDIDATES, collapse = "\n"))
}
message("Using occurrence file: ", OCCURRENCE_FILE)

## Environmental candidate stack
if (file.exists(ENV_CANDIDATE_TIF)) {
  message("Loading pre-VIF candidate environment stack: ", ENV_CANDIDATE_TIF)
  env_r <- terra::rast(ENV_CANDIDATE_TIF)
} else if (file.exists(ENV_CANDIDATE_FILELIST)) {
  message("Rebuilding environmental stack from file list: ", ENV_CANDIDATE_FILELIST)
  fl <- read_csv(ENV_CANDIDATE_FILELIST, show_col_types = FALSE)
  if (!all(c("variable", "cropped_file") %in% names(fl))) {
    stop("ENV_CANDIDATE_FILELIST must contain variable and cropped_file columns.")
  }
  env_r <- terra::rast(fl$cropped_file)
  names(env_r) <- safe_name(fl$variable)
} else {
  stop("No pre-VIF environmental candidate stack found.")
}

if (any(tolower(names(env_r)) == "swe")) {
  stop("swe is still in the candidate environmental stack. Use no-SWE candidate stack.")
}

message("Environmental layers for grid PCA: ", terra::nlyr(env_r))
read_csv(SPECIES_SUMMARY_FILE, show_col_types = FALSE) |>
  count(status) |>
  print(n = Inf)

sp_sum <- read_csv(SPECIES_SUMMARY_FILE, show_col_types = FALSE) |>
  filter(status %in% c("success", "skipped_existing")) |>
  filter(!is.na(pred_file), file.exists(pred_file)) |>
  mutate(
    genus = extract_genus(species),
    basic_group = classify_basic_group(genus)
  )

classification_summary <- sp_sum |>
  group_by(basic_group, genus) |>
  summarise(n_species_sdm = n_distinct(species), species_list = paste(sort(unique(species)), collapse = "; "), .groups = "drop") |>
  arrange(basic_group, genus)
write_csv(classification_summary, file.path(OUT_DIR, "00_basic_pollinator_classification.csv"))
message("Basic pollinator classification:")
print(classification_summary, n = Inf, width = Inf)

occ <- read_csv(OCCURRENCE_FILE, show_col_types = FALSE) |>
  mutate(
    decimalLongitude = to_numeric_coord(decimalLongitude),
    decimalLatitude = to_numeric_coord(decimalLatitude)
  ) |>
  filter(
    is.finite(decimalLongitude),
    is.finite(decimalLatitude),
    head_orientation_binary %in% c("upward", "nodding"),
    !is.na(species_final), species_final != ""
  ) |>
  mutate(y_nodding = ifelse(head_orientation_binary == "nodding", 1, 0))

message("Occurrence rows: ", nrow(occ))
print(occ |> count(head_orientation_binary))

## =========================================================
## 5. Build grid response and grid-level pre-VIF environmental PCA
## =========================================================

occ_grid <- occ |>
  mutate(
    grid_lon = floor(decimalLongitude / GRID_SIZE_DEG) * GRID_SIZE_DEG + GRID_SIZE_DEG / 2,
    grid_lat = floor(decimalLatitude / GRID_SIZE_DEG) * GRID_SIZE_DEG + GRID_SIZE_DEG / 2,
    grid_id = paste0(round(grid_lon, 3), "_", round(grid_lat, 3))
  )

grid_species <- occ_grid |>
  distinct(grid_id, grid_lon, grid_lat, species_final, japanese_name, head_orientation_binary, .keep_all = TRUE)

grid_counts <- grid_species |>
  group_by(grid_id, grid_lon, grid_lat) |>
  summarise(
    n_total_species = n_distinct(species_final),
    n_nodding_species = n_distinct(species_final[head_orientation_binary == "nodding"]),
    n_upward_species = n_distinct(species_final[head_orientation_binary == "upward"]),
    prop_nodding_species = n_nodding_species / n_total_species,
    has_nodding = as.integer(n_nodding_species > 0),
    .groups = "drop"
  )

## Extract environmental variables to occurrence points and average by grid.
xy_occ <- cbind(occ_grid$decimalLongitude, occ_grid$decimalLatitude)
colnames(xy_occ) <- c("lon", "lat")

env_vals_occ <- safe_extract_df(env_r, xy_occ)
names(env_vals_occ) <- paste0("preVIF_", safe_name(names(env_r)))

occ_env <- bind_cols(
  occ_grid |> select(grid_id, grid_lon, grid_lat),
  as_tibble(env_vals_occ)
)

env_cols <- names(env_vals_occ)

grid_env <- occ_env |>
  group_by(grid_id, grid_lon, grid_lat) |>
  summarise(across(all_of(env_cols), ~ mean(.x, na.rm = TRUE), .names = "{.col}"), .groups = "drop") |>
  mutate(across(where(is.numeric), ~ ifelse(is.nan(.x), NA_real_, .x)))

env_mat <- grid_env |> select(all_of(env_cols)) |> as.data.frame()

good_env <- names(env_mat)[vapply(env_mat, function(x) {
  sx <- sd(x, na.rm = TRUE)
  sum(is.finite(x)) >= 5 && is.finite(sx) && sx > 0
}, logical(1))]

env_mat <- env_mat[, good_env, drop = FALSE]
if (ncol(env_mat) < 2) stop("Too few environmental variables for grid PCA. Found: ", paste(good_env, collapse = ", "))

env_mat_imp <- replace_nonfinite_with_mean(env_mat)

env_pca <- prcomp(env_mat_imp, center = TRUE, scale. = TRUE)
saveRDS(env_pca, file.path(OUT_DIR, "01_grid_preVIF_env_PCA_model.rds"))

env_scores <- as.data.frame(env_pca$x[, seq_len(min(N_ENV_PC_USE, ncol(env_pca$x))), drop = FALSE]) |>
  as_tibble()
names(env_scores) <- paste0("env_PC", seq_len(ncol(env_scores)))

grid_base <- grid_counts |>
  left_join(grid_env |> select(grid_id, grid_lon, grid_lat), by = c("grid_id", "grid_lon", "grid_lat")) |>
  bind_cols(env_scores) |>
  mutate(
    z_lon = safe_scale(grid_lon),
    z_lat = safe_scale(grid_lat)
  )

for (nm in names(env_scores)) {
  grid_base[[paste0("z_", nm)]] <- safe_scale(grid_base[[nm]])
}

env_var_tbl <- tibble(
  PC = paste0("PC", seq_along(env_pca$sdev)),
  variance = env_pca$sdev^2,
  prop_variance = variance / sum(variance),
  cum_variance = cumsum(prop_variance)
)

env_loading_tbl <- as.data.frame(env_pca$rotation) |>
  rownames_to_column("variable") |>
  as_tibble()

env_top_loading <- env_loading_tbl |>
  pivot_longer(cols = starts_with("PC"), names_to = "PC", values_to = "loading") |>
  filter(PC %in% paste0("PC", seq_len(min(N_ENV_PC_SAVE_LOADINGS, ncol(env_pca$x))))) |>
  group_by(PC) |>
  slice_max(abs(loading), n = 15, with_ties = FALSE) |>
  arrange(PC, desc(abs(loading))) |>
  ungroup()

write_csv(grid_counts, file.path(OUT_DIR, "01_grid_response_counts.csv"))
write_csv(env_var_tbl, file.path(OUT_DIR, "01_grid_preVIF_env_PCA_variance_explained.csv"))
write_csv(env_loading_tbl, file.path(OUT_DIR, "01_grid_preVIF_env_PCA_loadings.csv"))
write_csv(env_top_loading, file.path(OUT_DIR, "01_grid_preVIF_env_PCA_top_loadings.csv"))

p_env_var <- env_var_tbl |>
  slice_head(n = 10) |>
  ggplot(aes(x = PC, y = prop_variance)) +
  geom_col() +
  theme_bw(base_size = 12) +
  labs(x = NULL, y = "Proportion variance explained", title = "Grid-level pre-VIF environmental PCA")
ggsave(file.path(PCA_DIR, "grid_preVIF_env_PCA_variance_explained.png"), p_env_var, width = 7, height = 4.5, dpi = 300)

## =========================================================
## 6. Build group raster stacks for all schemes
## =========================================================

build_group_stack <- function(scheme, group_name, sub_tbl) {
  scheme_dir <- file.path(STACK_DIR, safe_name(scheme))
  dir.create(scheme_dir, showWarnings = FALSE, recursive = TRUE)
  out_file <- file.path(scheme_dir, paste0("group_", safe_name(group_name), "_sum_cloglog.tif"))
  if (!FORCE_REBUILD_STACKS && file.exists(out_file)) return(out_file)
  pred_files <- sub_tbl$pred_file[file.exists(sub_tbl$pred_file)]
  if (length(pred_files) == 0) return(NA_character_)
  r <- terra::rast(pred_files)
  names(r) <- safe_name(tools::file_path_sans_ext(basename(pred_files)))
  sum_r <- terra::app(r, fun = sum, na.rm = TRUE)
  names(sum_r) <- paste0("group_", safe_name(group_name), "_sum_cloglog")
  terra::writeRaster(sum_r, out_file, overwrite = TRUE)
  out_file
}

scheme_manifest_list <- list()

grid_xy <- cbind(grid_base$grid_lon, grid_base$grid_lat)
colnames(grid_xy) <- c("lon", "lat")

for (scheme in SCHEMES) {
  message("Building/extracting pollinator scheme: ", scheme)

  scheme_tbl <- sp_sum |>
    mutate(scheme = .env$scheme, scheme_group = assign_scheme_group(basic_group, genus, .env$scheme)) |>
    filter(!scheme_group %in% c("unknown", "exclude"))

  if (nrow(scheme_tbl) == 0) {
    warning("No SDMs for scheme: ", scheme)
    next
  }

  group_summary <- scheme_tbl |>
    group_by(scheme, scheme_group, basic_group, genus) |>
    summarise(n_species_sdm = n_distinct(species), species_list = paste(sort(unique(species)), collapse = "; "), .groups = "drop") |>
    arrange(scheme_group, basic_group, genus)

  scheme_group_files <- scheme_tbl |>
    group_by(scheme, scheme_group) |>
    group_modify(~ {
      out_file <- build_group_stack(scheme = .y$scheme, group_name = .y$scheme_group, sub_tbl = .x)
      tibble(raster_file = out_file, n_species_sdm = n_distinct(.x$species), genera = paste(sort(unique(.x$genus)), collapse = "; "))
    }) |>
    ungroup()

  write_csv(group_summary, file.path(OUT_DIR, paste0("02_", safe_name(scheme), "_group_classification.csv")))
  write_csv(scheme_group_files, file.path(OUT_DIR, paste0("02_", safe_name(scheme), "_group_raster_files.csv")))

  scheme_manifest_list[[scheme]] <- scheme_group_files
}

scheme_manifest <- bind_rows(scheme_manifest_list)
write_csv(scheme_manifest, file.path(OUT_DIR, "02_all_scheme_group_raster_files.csv"))

## =========================================================
## 7. Analysis for each scheme
## =========================================================

run_scheme_analysis <- function(scheme, scheme_files, grid_base) {
  message("\n==============================")
  message("Running scheme: ", scheme)
  message("==============================")

  scheme_dir <- file.path(OUT_DIR, paste0("scheme_", safe_name(scheme)))
  dir.create(scheme_dir, showWarnings = FALSE, recursive = TRUE)

  files <- scheme_files |>
    filter(.data$scheme == .env$scheme, !is.na(raster_file), file.exists(raster_file)) |>
    arrange(scheme_group)

  if (nrow(files) == 0) {
    warning("No raster files for scheme: ", scheme)
    return(NULL)
  }

  r <- terra::rast(files$raster_file)
  raw_names <- paste0("poll_", safe_name(files$scheme_group), "_sum")
  names(r) <- raw_names

  vals <- safe_extract_df(r, grid_xy)
  names(vals) <- raw_names
  vals <- as_tibble(vals)

  dat <- bind_cols(grid_base, vals)
  for (nm in raw_names) {
    dat[[nm]][is.na(dat[[nm]])] <- 0
    dat[[paste0("z_", nm)]] <- safe_scale(dat[[nm]])
  }

  ## If this is the order PCA scheme, replace raw order variables with order PC variables.
  order_loading_tbl <- NULL
  order_var_tbl <- NULL
  if (scheme == "order3_pca") {
    order_mat <- dat |> select(all_of(raw_names)) |> as.data.frame()
    order_mat_imp <- replace_nonfinite_with_mean(order_mat)
    order_pca <- prcomp(order_mat_imp, center = TRUE, scale. = TRUE)
    saveRDS(order_pca, file.path(scheme_dir, "order_PCA_model.rds"))

    n_pc <- min(N_ORDER_PC_USE, ncol(order_pca$x))
    order_scores <- as.data.frame(order_pca$x[, seq_len(n_pc), drop = FALSE]) |> as_tibble()
    names(order_scores) <- paste0("order_PC", seq_len(n_pc))
    dat <- bind_cols(dat, order_scores)
    for (nm in names(order_scores)) dat[[paste0("z_", nm)]] <- safe_scale(dat[[nm]])

    order_loading_tbl <- as.data.frame(order_pca$rotation) |>
      rownames_to_column("variable") |>
      as_tibble() |>
      pivot_longer(cols = starts_with("PC"), names_to = "PC", values_to = "loading")

    order_var_tbl <- tibble(
      PC = paste0("PC", seq_along(order_pca$sdev)),
      variance = order_pca$sdev^2,
      prop_variance = variance / sum(variance),
      cum_variance = cumsum(prop_variance)
    )

    write_csv(order_loading_tbl, file.path(scheme_dir, "order_PCA_loadings.csv"))
    write_csv(order_var_tbl, file.path(scheme_dir, "order_PCA_variance_explained.csv"))

    p_load <- order_loading_tbl |>
      filter(PC %in% paste0("PC", seq_len(min(3, length(order_pca$sdev))))) |>
      ggplot(aes(x = loading, y = variable)) +
      geom_vline(xintercept = 0, linetype = 2) +
      geom_col() +
      facet_wrap(~ PC, scales = "free_x") +
      theme_bw(base_size = 11) +
      labs(x = "PCA loading", y = NULL, title = "Pollinator order PCA loadings")
    ggsave(file.path(scheme_dir, "order_PCA_loadings.png"), p_load, width = 10, height = 5.5, dpi = 300)
  }

  env_vars <- paste0("z_env_PC", seq_len(N_ENV_PC_USE))
  env_vars <- env_vars[env_vars %in% names(dat)]
  raw_poll_vars <- paste0("z_", raw_names)
  raw_poll_vars <- raw_poll_vars[raw_poll_vars %in% names(dat)]
  if (scheme == "order3_pca") {
    poll_vars <- grep("^z_order_PC", names(dat), value = TRUE)
  } else {
    poll_vars <- raw_poll_vars
  }

  dat <- dat |>
    filter(n_total_species >= MIN_TOTAL_SPECIES_PER_GRID) |>
    drop_na(n_nodding_species, n_upward_species, n_total_species, prop_nodding_species, all_of(env_vars), z_lon, z_lat)

  ## Standardized expected excess for maps.
  global_prop <- sum(dat$n_nodding_species) / sum(dat$n_total_species)
  dat <- dat |>
    mutate(
      expected_nodding_global = n_total_species * global_prop,
      excess_nodding = n_nodding_species - expected_nodding_global,
      standardized_excess_nodding = excess_nodding / sqrt(pmax(expected_nodding_global, 1e-6))
    )

  write_csv(dat, file.path(scheme_dir, "analysis_grid_data.csv"))

  ## Correlations among predictors and response.
  cor_vars <- unique(c(env_vars, poll_vars, "z_lon", "z_lat", "prop_nodding_species"))
  cor_tbl <- cor_pairs_table(dat, cor_vars, scheme)
  write_csv(cor_tbl, file.path(scheme_dir, "predictor_correlation_pairs.csv"))

  dangerous_tbl <- cor_tbl |> filter(abs_correlation >= COR_WARN_THRESHOLD)
  if (nrow(dangerous_tbl) > 0) {
    message("High correlations for ", scheme, ":")
    print(dangerous_tbl, n = Inf, width = Inf)
  }

  rhs_env <- get_rhs(env_vars)
  rhs_poll <- get_rhs(poll_vars)
  rhs_lat <- "z_lat"
  rhs_lonlat <- "z_lon + z_lat"
  rhs_space_poly <- "z_lon + z_lat + I(z_lon^2) + I(z_lat^2) + z_lon:z_lat"
  rhs_env_poll <- get_rhs(c(env_vars, poll_vars))
  rhs_env_poll_lat <- get_rhs(c(env_vars, poll_vars, "z_lat"))
  rhs_env_poll_space <- get_rhs(c(env_vars, poll_vars, "z_lon", "z_lat", "I(z_lon^2)", "I(z_lat^2)", "z_lon:z_lat"))
  rhs_poll_space <- get_rhs(c(poll_vars, "z_lon", "z_lat", "I(z_lon^2)", "I(z_lat^2)", "z_lon:z_lat"))
  rhs_env_space <- get_rhs(c(env_vars, "z_lon", "z_lat", "I(z_lon^2)", "I(z_lat^2)", "z_lon:z_lat"))

  model_specs <- list(
    null = "1",
    lat_only = rhs_lat,
    lonlat = rhs_lonlat,
    space_poly = rhs_space_poly,
    env = rhs_env,
    poll = rhs_poll,
    env_plus_poll = rhs_env_poll,
    env_plus_poll_plus_lat = rhs_env_poll_lat,
    env_plus_space_poly = rhs_env_space,
    poll_plus_space_poly = rhs_poll_space,
    env_plus_poll_plus_space_poly = rhs_env_poll_space
  )

  ## Remove malformed specs.
  model_specs <- model_specs[!vapply(model_specs, function(x) is.na(x) || !nzchar(x), logical(1))]

  glm_results <- list()
  bb_results <- list()

  for (mn in names(model_specs)) {
    f <- make_binom_formula(model_specs[[mn]])
    message("  fitting ", scheme, " / ", mn)
    glm_results[[mn]] <- fit_glm_binom(f, dat, scheme, mn)
    bb_results[[mn]] <- fit_glmmTMB_betabinom(f, dat, scheme, mn)
  }

  glm_summary <- bind_rows(lapply(glm_results, `[[`, "summary")) |> add_delta_aic()
  bb_summary <- bind_rows(lapply(bb_results, `[[`, "summary")) |> add_delta_aic()
  glm_coef <- bind_rows(lapply(glm_results, `[[`, "coef"))
  bb_coef <- bind_rows(lapply(bb_results, `[[`, "coef"))

  write_csv(glm_summary, file.path(scheme_dir, "binomial_GLM_model_comparison.csv"))
  write_csv(bb_summary, file.path(scheme_dir, "beta_binomial_model_comparison.csv"))
  write_csv(glm_coef, file.path(scheme_dir, "binomial_GLM_coefficients.csv"))
  write_csv(bb_coef, file.path(scheme_dir, "beta_binomial_coefficients.csv"))

  ## Best method/model by beta-binomial if available; otherwise binomial GLM.
  best_method <- "beta_binomial_glmmTMB"
  best_tbl <- bb_summary |> filter(is.finite(AIC)) |> arrange(AIC)
  if (nrow(best_tbl) == 0) {
    best_method <- "binomial_GLM"
    best_tbl <- glm_summary |> filter(is.finite(AIC)) |> arrange(AIC)
  }

  best_model_name <- if (nrow(best_tbl) > 0) best_tbl$model[1] else NA_character_
  best_model_obj <- NULL
  best_coef <- tibble()

  if (!is.na(best_model_name)) {
    if (best_method == "beta_binomial_glmmTMB") {
      best_model_obj <- bb_results[[best_model_name]]$model
      best_coef <- bb_coef |> filter(model == best_model_name)
    } else {
      best_model_obj <- glm_results[[best_model_name]]$model
      best_coef <- glm_coef |> filter(model == best_model_name)
    }
  }

  if (!is.null(best_model_obj)) {
    dat$pred_prob_nodding_best <- as.numeric(predict(best_model_obj, newdata = dat, type = "response"))
    dat$resid_prop_nodding_best <- dat$prop_nodding_species - dat$pred_prob_nodding_best
    dat$best_model_name <- best_model_name
    dat$best_model_method <- best_method
    write_csv(dat, file.path(scheme_dir, "predicted_probability_best_model.csv"))
    saveRDS(best_model_obj, file.path(scheme_dir, paste0("best_model_", best_method, "_", best_model_name, ".rds")))

    write_csv(best_coef, file.path(scheme_dir, "best_model_coefficients.csv"))
    plot_coef(best_coef, file.path(scheme_dir, "best_model_coefficients.png"), paste0("Best model: ", scheme, " / ", best_model_name))
    plot_map(dat, "pred_prob_nodding_best", file.path(scheme_dir, "map_predicted_probability_best_model.png"), paste0("Predicted P(nodding): ", scheme))
    plot_map(dat, "resid_prop_nodding_best", file.path(scheme_dir, "map_residual_best_model.png"), paste0("Observed - predicted: ", scheme), diverging = TRUE)
    plot_map(dat, "standardized_excess_nodding", file.path(scheme_dir, "map_standardized_excess_nodding.png"), paste0("Standardized excess nodding: ", scheme), diverging = TRUE)
  }

  list(
    scheme = scheme,
    data = dat,
    glm_summary = glm_summary,
    bb_summary = bb_summary,
    glm_coef = glm_coef,
    bb_coef = bb_coef,
    cor_tbl = cor_tbl,
    best = if (nrow(best_tbl) > 0) best_tbl[1, ] else tibble()
  )
}

scheme_results <- list()
for (scheme in SCHEMES) {
  files <- scheme_manifest |> filter(.data$scheme == .env$scheme)
  if (nrow(files) == 0) next
  scheme_results[[scheme]] <- run_scheme_analysis(scheme, files, grid_base)
}

## =========================================================
## 8. Combined outputs
## =========================================================

all_glm_summary <- bind_rows(lapply(scheme_results, `[[`, "glm_summary")) |>
  group_by(method) |>
  mutate(
    global_delta_AIC = ifelse(is.finite(AIC), AIC - min(AIC, na.rm = TRUE), NA_real_),
    global_akaike_weight = ifelse(
      is.finite(global_delta_AIC),
      exp(-0.5 * global_delta_AIC) / sum(exp(-0.5 * global_delta_AIC), na.rm = TRUE),
      NA_real_
    )
  ) |>
  ungroup() |>
  arrange(AIC)

all_bb_summary <- bind_rows(lapply(scheme_results, `[[`, "bb_summary")) |>
  group_by(method) |>
  mutate(
    global_delta_AIC = ifelse(is.finite(AIC), AIC - min(AIC, na.rm = TRUE), NA_real_),
    global_akaike_weight = ifelse(
      is.finite(global_delta_AIC),
      exp(-0.5 * global_delta_AIC) / sum(exp(-0.5 * global_delta_AIC), na.rm = TRUE),
      NA_real_
    )
  ) |>
  ungroup() |>
  arrange(AIC)

all_glm_coef <- bind_rows(lapply(scheme_results, `[[`, "glm_coef"))
all_bb_coef <- bind_rows(lapply(scheme_results, `[[`, "bb_coef"))
all_cor <- bind_rows(lapply(scheme_results, `[[`, "cor_tbl"))
all_best <- bind_rows(lapply(scheme_results, `[[`, "best")) |>
  arrange(AIC)

write_csv(all_glm_summary, file.path(OUT_DIR, "all_scheme_binomial_GLM_model_comparison.csv"))
write_csv(all_bb_summary, file.path(OUT_DIR, "all_scheme_beta_binomial_model_comparison.csv"))
write_csv(all_glm_coef, file.path(OUT_DIR, "all_scheme_binomial_GLM_coefficients.csv"))
write_csv(all_bb_coef, file.path(OUT_DIR, "all_scheme_beta_binomial_coefficients.csv"))
write_csv(all_cor, file.path(OUT_DIR, "all_scheme_predictor_correlations.csv"))
write_csv(all_best, file.path(OUT_DIR, "all_scheme_best_models.csv"))

## High-correlation summary only.
high_cor <- all_cor |> filter(abs_correlation >= COR_WARN_THRESHOLD)
write_csv(high_cor, file.path(OUT_DIR, "all_scheme_high_predictor_correlations_ge_0p70.csv"))

## Top model tables.
top_bb <- all_bb_summary |> filter(is.finite(AIC)) |> arrange(AIC) |> slice_head(n = 30)
top_glm <- all_glm_summary |> filter(is.finite(AIC)) |> arrange(AIC) |> slice_head(n = 30)
write_csv(top_bb, file.path(OUT_DIR, "top30_beta_binomial_models_across_schemes.csv"))
write_csv(top_glm, file.path(OUT_DIR, "top30_binomial_GLM_models_across_schemes.csv"))

plot_model_compare(
  all_bb_summary |> filter(is.finite(AIC)) |> mutate(delta_AIC = global_delta_AIC),
  file.path(FIG_DIR, "model_comparison_all_schemes_beta_binomial.png"),
  "All pollinator grouping schemes: beta-binomial model comparison"
)

plot_model_compare(
  all_glm_summary |> filter(is.finite(AIC)) |> mutate(delta_AIC = global_delta_AIC),
  file.path(FIG_DIR, "model_comparison_all_schemes_binomial_GLM.png"),
  "All pollinator grouping schemes: binomial GLM model comparison"
)

## Compact bar: best model per scheme.
if (nrow(all_best) > 0) {
  p_best <- all_best |>
    filter(is.finite(AIC)) |>
    mutate(label = paste(scheme, model, sep = " / ")) |>
    arrange(AIC) |>
    mutate(label = factor(label, levels = rev(label))) |>
    ggplot(aes(x = AIC - min(AIC, na.rm = TRUE), y = label)) +
    geom_col() +
    theme_bw(base_size = 11) +
    labs(x = "Delta AIC among best models", y = NULL, title = "Best model from each pollinator grouping scheme")
  ggsave(file.path(FIG_DIR, "best_model_each_scheme_delta_AIC.png"), p_best, width = 10, height = 6, dpi = 300)
}

## =========================================================
## 9. Console summary
## =========================================================

message("\n🎉 DONE — pollinator grouping sensitivity model comparison")
message("Output dir: ", OUT_DIR)
message("\nTop beta-binomial models across schemes:")
print(top_bb, n = 20, width = Inf)
message("\nHigh correlations >= ", COR_WARN_THRESHOLD, ":")
print(high_cor, n = Inf, width = Inf)
message("\nKey outputs:")
message("  all_scheme_beta_binomial_model_comparison.csv")
message("  all_scheme_binomial_GLM_model_comparison.csv")
message("  all_scheme_best_models.csv")
message("  all_scheme_predictor_correlations.csv")
message("  all_scheme_high_predictor_correlations_ge_0p70.csv")
message("  top30_beta_binomial_models_across_schemes.csv")
message("  figures/model_comparison_all_schemes_beta_binomial.png")
message("  figures/best_model_each_scheme_delta_AIC.png")
message("\nEach scheme has its own folder: scheme_<scheme_name>/")
message("  analysis_grid_data.csv")
message("  beta_binomial_model_comparison.csv")
message("  best_model_coefficients.csv")
message("  predicted_probability_best_model.csv")
message("  map_predicted_probability_best_model.png")
message("  map_standardized_excess_nodding.png")
