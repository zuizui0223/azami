############################################################
## Global Cirsium head-orientation analysis
## FULL rerun + publication-level figure export
##
## Purpose:
##   Re-run the full analysis from the existing enriched table:
##     1. RF without space
##     2. RF with space
##     3. RF-top variable selection
##     4. correlation clustering
##     5. VIF pruning
##     6. common-N weighted GLMs
##     7. incremental likelihood-ratio tests
##     8. GAMs
##     9. publication-level figures in PNG/PDF/TIFF
##
## Important:
##   This script does NOT download CHELSA, topography, or SoilGrids.
##   It assumes those values were already extracted into:
##
##   C:/Users/zuizui/cirsium_inat/
##     global_cirsium_CHELSA_SOILGRIDS_TOPO_RF_GLM_GAM/
##       01_orientation_points_CHELSA_topography_soilgrids.csv
##
## Run:
##   setwd("C:/Users/zuizui/cirsium_inat")
##   source("global_cirsium_full_rerun_RF_GLM_GAM_publication_figures.R")
############################################################

## =========================================================
## 0. Packages
## =========================================================

pkgs <- c(
  "dplyr", "readr", "stringr", "tibble", "tidyr",
  "ggplot2", "forcats", "patchwork", "broom",
  "randomForest", "pROC", "car", "mgcv", "viridis",
  "scales", "purrr"
)

for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, dependencies = TRUE)
  }
}

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
  library(tidyr)
  library(ggplot2)
  library(forcats)
  library(patchwork)
  library(broom)
  library(randomForest)
  library(pROC)
  library(car)
  library(mgcv)
  library(viridis)
  library(scales)
  library(purrr)
})

## Avoid common conflicts
select <- dplyr::select
filter <- dplyr::filter
mutate <- dplyr::mutate
summarise <- dplyr::summarise
arrange <- dplyr::arrange
group_by <- dplyr::group_by
ungroup <- dplyr::ungroup
left_join <- dplyr::left_join
bind_rows <- dplyr::bind_rows
slice <- dplyr::slice
slice_max <- dplyr::slice_max
slice_sample <- dplyr::slice_sample
all_of <- dplyr::all_of
any_of <- dplyr::any_of
across <- dplyr::across

set.seed(42)

## =========================================================
## 1. Paths and settings
## =========================================================

ROOT_DIR <- "C:/Users/zuizui/cirsium_inat"

BASE_OUT_DIR <- file.path(
  ROOT_DIR,
  "global_cirsium_CHELSA_SOILGRIDS_TOPO_RF_GLM_GAM"
)

ENRICHED_FILE <- file.path(
  BASE_OUT_DIR,
  "01_orientation_points_CHELSA_topography_soilgrids.csv"
)

NO_RANGE_VAR_FILE <- file.path(
  ROOT_DIR,
  "global_cirsium_ALL_CHELSA_RF_no_minmaxrange",
  "02_variables_used_no_minmaxrange.csv"
)

## Fresh full rerun output directory
OUT_DIR <- file.path(
  ROOT_DIR,
  "global_cirsium_FULL_RERUN_RF_GLM_GAM_PUBLICATION"
)

RF_DIR <- file.path(OUT_DIR, "01_RF")
MODEL_DIR <- file.path(OUT_DIR, "02_GLM_GAM")
DIAG_DIR <- file.path(OUT_DIR, "03_diagnostics")
PUB_DIR <- file.path(OUT_DIR, "04_publication_figures")
PUB_PNG_DIR <- file.path(PUB_DIR, "png")
PUB_PDF_DIR <- file.path(PUB_DIR, "pdf")
PUB_TIFF_DIR <- file.path(PUB_DIR, "tiff")
PUB_DATA_DIR <- file.path(PUB_DIR, "figure_data")

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(RF_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(MODEL_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(DIAG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(PUB_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(PUB_PNG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(PUB_PDF_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(PUB_TIFF_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(PUB_DATA_DIR, showWarnings = FALSE, recursive = TRUE)

## Topography variables expected in enriched file
TOPO_VARS <- c("elevation", "slope", "roughness")

## Soil variables as 0-30 cm composite columns expected in enriched file
SOIL_VARS <- c(
  "soil_bdod_0_30cm",
  "soil_cec_0_30cm",
  "soil_cfvo_0_30cm",
  "soil_clay_0_30cm",
  "soil_sand_0_30cm",
  "soil_silt_0_30cm",
  "soil_nitrogen_0_30cm",
  "soil_phh2o_0_30cm",
  "soil_soc_0_30cm",
  "soil_ocd_0_30cm"
)

## RF settings
RF_NTREE_MAIN <- 1500
RF_NTREE_REPEAT <- 800
RF_REPEATS <- 20
RF_MTRY <- NULL
BALANCE_CLASSES_FOR_RF <- TRUE
MAX_OCC_RF_ROWS <- 50000
MAX_PROP_NA_PREDICTOR <- 0.20

## Variable selection / GLM / GAM settings
TOP_N_RF_FOR_MODELS <- 30
COR_CLUSTER_THRESHOLD <- 0.70
MAX_REP_VARS <- 14
VIF_THRESHOLD <- 5
USE_CLASS_BALANCED_WEIGHTS <- TRUE
FORCE_TOPOGRAPHY_IN_GLM <- TRUE
FORCE_SOIL_SUMMARY_IN_GLM <- TRUE

GAM_MAX_N <- 25000
GAM_ENV_K <- 5

SPACE_VARS <- c("z_lon", "z_lat")
SPACE_POLY_VARS <- c("z_lon", "z_lat", "z_lon2", "z_lat2", "z_lon_lat")

## Publication style
SHOW_PLOT_TITLES <- FALSE
BASE_SIZE <- 10.5
DPI <- 450

## =========================================================
## 2. Helper functions
## =========================================================

msg <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", ..., "\n")
}

stop_if_missing <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)
}

safe_scale <- function(x) {
  x <- as.numeric(x)
  s <- stats::sd(x, na.rm = TRUE)
  m <- mean(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(0, length(x)))
  (x - m) / s
}

auc_safe <- function(y, p) {
  tryCatch(as.numeric(pROC::auc(y, p, quiet = TRUE)), error = function(e) NA_real_)
}

choose_mtry <- function(p) {
  if (is.null(RF_MTRY) || length(RF_MTRY) == 0 || is.na(RF_MTRY)) {
    return(max(1L, floor(sqrt(p))))
  }
  m <- as.integer(RF_MTRY)
  if (!is.finite(m) || is.na(m)) return(max(1L, floor(sqrt(p))))
  max(1L, min(m, p))
}

make_class_weights <- function(y01) {
  tab <- table(y01)
  w <- rep(1, length(y01))
  if (length(tab) < 2) return(w)
  n0 <- as.numeric(tab[["0"]])
  n1 <- as.numeric(tab[["1"]])
  n <- length(y01)
  w[y01 == 0] <- n / (2 * n0)
  w[y01 == 1] <- n / (2 * n1)
  w
}

balance_binary_df <- function(df, y_col = "orientation_factor") {
  tab <- table(df[[y_col]])
  if (length(tab) < 2) return(df)
  n_min <- min(tab)

  df |>
    group_by(.data[[y_col]]) |>
    slice_sample(n = n_min) |>
    ungroup()
}

weighted_logloss <- function(y, p, w = NULL, eps = 1e-8) {
  p <- pmin(pmax(p, eps), 1 - eps)
  if (is.null(w)) w <- rep(1, length(y))
  -sum(w * (y * log(p) + (1 - y) * log(1 - p)), na.rm = TRUE) / sum(w, na.rm = TRUE)
}

safe_vif_table <- function(model) {
  tryCatch({
    v <- car::vif(model)
    tibble(variable = names(v), VIF = as.numeric(v)) |>
      arrange(desc(VIF))
  }, error = function(e) {
    tibble(variable = NA_character_, VIF = NA_real_, error = conditionMessage(e))
  })
}

safe_auc_col <- function(y, p) {
  if (length(unique(y[is.finite(y)])) < 2) return(NA_real_)
  auc_safe(y, p)
}

theme_paper <- function(base_size = BASE_SIZE) {
  theme_classic(base_size = base_size) +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      panel.grid.major.x = element_line(colour = "grey88", linewidth = 0.25),
      panel.grid.major.y = element_blank(),
      legend.title = element_text(face = "bold"),
      legend.key = element_blank(),
      strip.background = element_rect(fill = "grey92", colour = "grey60", linewidth = 0.35),
      strip.text = element_text(face = "bold"),
      plot.title = if (SHOW_PLOT_TITLES) element_text(face = "bold", hjust = 0) else element_blank(),
      plot.subtitle = if (SHOW_PLOT_TITLES) element_text(hjust = 0) else element_blank()
    )
}

paper_labs <- function(...) {
  args <- list(...)
  if (!isTRUE(SHOW_PLOT_TITLES)) {
    args$title <- NULL
    args$subtitle <- NULL
  }
  do.call(ggplot2::labs, args)
}

save_pub <- function(plot, name, width = 7, height = 5) {
  png_file <- file.path(PUB_PNG_DIR, paste0(name, ".png"))
  pdf_file <- file.path(PUB_PDF_DIR, paste0(name, ".pdf"))
  tiff_file <- file.path(PUB_TIFF_DIR, paste0(name, ".tiff"))

  ggsave(png_file, plot, width = width, height = height, dpi = DPI, bg = "white")
  ggsave(pdf_file, plot, width = width, height = height, bg = "white")
  ggsave(tiff_file, plot, width = width, height = height, dpi = DPI, bg = "white", compression = "lzw")

  tibble(
    figure = name,
    png = png_file,
    pdf = pdf_file,
    tiff = tiff_file,
    exists_png = file.exists(png_file),
    exists_pdf = file.exists(pdf_file),
    exists_tiff = file.exists(tiff_file)
  )
}

save_data <- function(x, name) {
  readr::write_csv(x, file.path(PUB_DATA_DIR, paste0(name, ".csv")))
}

pretty_var <- function(x) {
  x0 <- stringr::str_remove(as.character(x), "^z_")
  dplyr::case_when(
    x0 == "lon" ~ "Longitude",
    x0 == "lat" ~ "Latitude",
    x0 == "z_lon" ~ "Longitude",
    x0 == "z_lat" ~ "Latitude",
    x0 == "z_lon2" ~ "Longitude²",
    x0 == "z_lat2" ~ "Latitude²",
    x0 == "z_lon_lat" ~ "Longitude × latitude",
    x0 == "longitude" ~ "Longitude",
    x0 == "latitude" ~ "Latitude",
    x0 == "elevation" ~ "Elevation",
    x0 == "slope" ~ "Slope",
    x0 == "roughness" ~ "Roughness",
    x0 == "soil_bdod_0_30cm" ~ "Soil bulk density",
    x0 == "soil_cec_0_30cm" ~ "Soil CEC",
    x0 == "soil_cfvo_0_30cm" ~ "Soil coarse fragments",
    x0 == "soil_clay_0_30cm" ~ "Soil clay",
    x0 == "soil_sand_0_30cm" ~ "Soil sand",
    x0 == "soil_silt_0_30cm" ~ "Soil silt",
    x0 == "soil_nitrogen_0_30cm" ~ "Soil nitrogen",
    x0 == "soil_phh2o_0_30cm" ~ "Soil pH",
    x0 == "soil_soc_0_30cm" ~ "Soil organic carbon",
    x0 == "soil_ocd_0_30cm" ~ "Soil organic carbon density",
    TRUE ~ x0 |>
      stringr::str_replace_all("_", " ") |>
      stringr::str_replace_all("BIO", "BIO")
  )
}

class_of_var <- function(v) {
  vv <- stringr::str_remove(as.character(v), "^z_")
  dplyr::case_when(
    vv %in% TOPO_VARS ~ "Topography",
    vv %in% SOIL_VARS ~ "Soil",
    vv %in% c("lon", "lat", "longitude", "latitude", "z_lon", "z_lat", "z_lon2", "z_lat2", "z_lon_lat") ~ "Space",
    vv %in% c("z_lon", "z_lat", "z_lon2", "z_lat2", "z_lon_lat") ~ "Space",
    TRUE ~ "Climate"
  )
}

make_safe_x <- function(df, predictors) {
  predictors <- unique(as.character(predictors))
  predictors <- predictors[predictors %in% names(df)]

  if (length(predictors) == 0) {
    stop("No predictors remain after matching to data columns.")
  }

  x <- as.data.frame(df[, predictors, drop = FALSE])
  original_names <- names(x)
  safe_names <- make.names(original_names, unique = TRUE)
  names(x) <- safe_names

  list(
    x = x,
    predictors_original = original_names,
    predictors_safe = safe_names
  )
}

geom_ci_horizontal <- function(mapping = NULL, data = NULL, height = 0.15, linewidth = 0.4, ...) {
  ## Avoids geom_errorbarh() warning: "`height` was translated to `width`."
  ggplot2::geom_errorbar(
    mapping = mapping,
    data = data,
    orientation = "y",
    width = height,
    linewidth = linewidth,
    ...
  )
}

## =========================================================
## 3. Random forest function
## =========================================================

run_rf_classification <- function(df, predictors, y_col, prefix, out_dir,
                                  ntree = RF_NTREE_MAIN,
                                  repeats = RF_REPEATS) {
  predictors <- unique(as.character(predictors))
  predictors <- predictors[predictors %in% names(df)]

  keep_cols <- unique(c(y_col, predictors))
  rf_df <- as.data.frame(df[, keep_cols, drop = FALSE])
  rf_df <- tidyr::drop_na(rf_df)

  rf_df[[y_col]] <- factor(rf_df[[y_col]], levels = c("upward", "nodding"))

  if (length(unique(rf_df[[y_col]])) < 2) {
    warning("RF skipped for ", prefix, ": only one class.")
    return(NULL)
  }

  if (nrow(rf_df) > MAX_OCC_RF_ROWS && stringr::str_detect(prefix, "^occurrence")) {
    n_each <- floor(MAX_OCC_RF_ROWS / 2)
    rf_df <- rf_df |>
      group_by(.data[[y_col]]) |>
      group_modify(function(.x, .y) {
        slice_sample(.x, n = min(n_each, nrow(.x)))
      }) |>
      ungroup()
  }

  rf_train <- if (BALANCE_CLASSES_FOR_RF) balance_binary_df(rf_df, y_col) else rf_df

  xx <- make_safe_x(rf_train, predictors)
  x_train <- xx$x
  y_train <- rf_train[[y_col]]
  mtry_used <- choose_mtry(ncol(x_train))

  msg("RF ", prefix)
  msg("  original n = ", nrow(rf_df))
  msg("  used n     = ", nrow(rf_train))
  msg("  predictors = ", ncol(x_train))
  msg("  mtry       = ", mtry_used)
  msg("  balanced   = ", BALANCE_CLASSES_FOR_RF)
  print(table(y_train))

  rf <- randomForest::randomForest(
    x = x_train,
    y = y_train,
    ntree = ntree,
    mtry = mtry_used,
    importance = TRUE,
    proximity = FALSE
  )

  pred_prob <- predict(rf, x_train, type = "prob")[, "nodding"]
  pred_class <- predict(rf, x_train, type = "response")

  conf <- as.data.frame.matrix(table(observed = y_train, predicted = pred_class)) |>
    tibble::rownames_to_column("observed") |>
    as_tibble()

  metrics <- tibble(
    analysis = prefix,
    n_original_before_RF = nrow(rf_df),
    n_used_for_RF = nrow(rf_train),
    class_balanced = BALANCE_CLASSES_FOR_RF,
    n_predictors = ncol(x_train),
    mtry = mtry_used,
    OOB_error = rf$err.rate[nrow(rf$err.rate), "OOB"],
    AUC_train = auc_safe(y_train, pred_prob)
  )

  imp <- as.data.frame(randomForest::importance(rf)) |>
    tibble::rownames_to_column("variable_safe") |>
    as_tibble() |>
    mutate(
      variable = xx$predictors_original[match(variable_safe, xx$predictors_safe)],
      .after = variable_safe
    ) |>
    arrange(desc(MeanDecreaseAccuracy))

  write_csv(metrics, file.path(out_dir, paste0(prefix, "_RF_metrics.csv")))
  write_csv(conf, file.path(out_dir, paste0(prefix, "_RF_confusion_matrix.csv")))
  write_csv(imp, file.path(out_dir, paste0(prefix, "_RF_importance.csv")))
  saveRDS(rf, file.path(out_dir, paste0(prefix, "_RF_model.rds")))

  rep_list <- list()

  for (i in seq_len(repeats)) {
    set.seed(1000 + i)

    train_i <- if (BALANCE_CLASSES_FOR_RF) balance_binary_df(rf_df, y_col) else rf_df

    xxi <- make_safe_x(train_i, predictors)
    x_i <- xxi$x
    y_i <- train_i[[y_col]]

    fit_i <- randomForest::randomForest(
      x = x_i,
      y = y_i,
      ntree = RF_NTREE_REPEAT,
      mtry = choose_mtry(ncol(x_i)),
      importance = TRUE,
      proximity = FALSE
    )

    imp_i <- as.data.frame(randomForest::importance(fit_i)) |>
      tibble::rownames_to_column("variable_safe") |>
      as_tibble() |>
      mutate(
        variable = xxi$predictors_original[match(variable_safe, xxi$predictors_safe)],
        repeat_id = i,
        .before = 1
      )

    rep_list[[i]] <- imp_i
  }

  imp_rep <- bind_rows(rep_list)

  imp_summary <- imp_rep |>
    group_by(variable) |>
    summarise(
      mean_MDA = mean(MeanDecreaseAccuracy, na.rm = TRUE),
      sd_MDA = sd(MeanDecreaseAccuracy, na.rm = TRUE),
      prop_positive_MDA = mean(MeanDecreaseAccuracy > 0, na.rm = TRUE),
      mean_MDG = mean(MeanDecreaseGini, na.rm = TRUE),
      .groups = "drop"
    ) |>
    arrange(desc(mean_MDA))

  write_csv(imp_rep, file.path(out_dir, paste0(prefix, "_RF_importance_repeats_raw.csv")))
  write_csv(imp_summary, file.path(out_dir, paste0(prefix, "_RF_importance_repeats_summary.csv")))

  list(model = rf, metrics = metrics, importance = imp_summary)
}

## =========================================================
## 4. Variable selection helpers
## =========================================================

select_cluster_representatives <- function(dat, candidate_vars, importance_tbl,
                                           threshold = 0.70,
                                           max_vars = 14,
                                           force_keep = character()) {
  candidate_vars <- candidate_vars[candidate_vars %in% names(dat)]

  if (length(candidate_vars) <= 1) {
    return(list(
      reps = candidate_vars,
      clusters = tibble(variable = candidate_vars, cluster = 1L, mean_MDA = NA_real_),
      cor = matrix(NA_real_, nrow = length(candidate_vars), ncol = length(candidate_vars))
    ))
  }

  x <- dat[, candidate_vars, drop = FALSE] |>
    mutate(across(everything(), as.numeric)) |>
    drop_na()

  keep <- names(x)[vapply(x, function(z) stats::sd(z, na.rm = TRUE) > 0, logical(1))]
  x <- x[, keep, drop = FALSE]
  candidate_vars <- keep

  cor_mat <- stats::cor(x, use = "pairwise.complete.obs")
  cor_mat[is.na(cor_mat)] <- 0

  hc <- stats::hclust(stats::as.dist(1 - abs(cor_mat)), method = "average")
  cl <- stats::cutree(hc, h = 1 - threshold)

  cl_tbl <- tibble(variable = names(cl), cluster = as.integer(cl)) |>
    left_join(importance_tbl |> select(variable, mean_MDA), by = "variable") |>
    mutate(
      mean_MDA = ifelse(is.na(mean_MDA), 0, mean_MDA),
      is_forced = variable %in% force_keep
    )

  reps_by_cluster <- cl_tbl |>
    group_by(cluster) |>
    slice_max(mean_MDA, n = 1, with_ties = FALSE) |>
    ungroup() |>
    arrange(desc(mean_MDA)) |>
    slice(seq_len(min(max_vars, n()))) |>
    pull(variable)

  forced <- intersect(force_keep, candidate_vars)
  reps <- unique(c(reps_by_cluster, forced))

  list(reps = reps, clusters = cl_tbl, cor = cor_mat, hc = hc)
}

vif_prune <- function(dat, response, predictors, weights_col = NULL, threshold = 5) {
  predictors <- predictors[predictors %in% names(dat)]
  removed <- tibble(step = integer(), removed_variable = character(), max_VIF = numeric())
  step <- 0

  repeat {
    if (length(predictors) <= 1) break

    step <- step + 1
    cols <- unique(c(response, predictors, weights_col))
    d <- dat[, cols, drop = FALSE] |> drop_na()

    f <- as.formula(paste(response, "~", paste(predictors, collapse = " + ")))

    m <- tryCatch({
      if (!is.null(weights_col) && weights_col %in% names(d)) {
        glm(f, data = d, family = binomial, weights = d[[weights_col]])
      } else {
        glm(f, data = d, family = binomial)
      }
    }, error = function(e) NULL)

    if (is.null(m)) break

    vt <- safe_vif_table(m)
    if (!"VIF" %in% names(vt) || all(is.na(vt$VIF))) break

    max_vif <- max(vt$VIF, na.rm = TRUE)

    if (!is.finite(max_vif) || max_vif <= threshold) {
      return(list(kept = predictors, removed = removed, final_vif = vt))
    }

    drop_var <- vt$variable[which.max(vt$VIF)]

    removed <- bind_rows(
      removed,
      tibble(step = step, removed_variable = drop_var, max_VIF = max_vif)
    )

    predictors <- setdiff(predictors, drop_var)
  }

  final_vif <- if (length(predictors) > 1) {
    cols <- unique(c(response, predictors, weights_col))
    d <- dat[, cols, drop = FALSE] |> drop_na()
    f <- as.formula(paste(response, "~", paste(predictors, collapse = " + ")))

    m <- tryCatch({
      if (!is.null(weights_col) && weights_col %in% names(d)) {
        glm(f, data = d, family = binomial, weights = d[[weights_col]])
      } else {
        glm(f, data = d, family = binomial)
      }
    }, error = function(e) NULL)

    if (is.null(m)) tibble(variable = predictors, VIF = NA_real_) else safe_vif_table(m)
  } else {
    tibble(variable = predictors, VIF = NA_real_)
  }

  list(kept = predictors, removed = removed, final_vif = final_vif)
}

fit_glm_set <- function(dat, response, vars, prefix, weights_col = NULL) {
  vars <- vars[vars %in% names(dat)]
  if (length(vars) == 0) stop("No vars for model: ", prefix)

  cols <- unique(c(response, vars, weights_col))
  dat2 <- dat[, cols, drop = FALSE] |> drop_na()

  f <- as.formula(paste(response, "~", paste(vars, collapse = " + ")))

  if (!is.null(weights_col) && weights_col %in% names(dat2)) {
    m <- glm(f, data = dat2, family = binomial, weights = dat2[[weights_col]])
    ww <- dat2[[weights_col]]
  } else {
    m <- glm(f, data = dat2, family = binomial)
    ww <- rep(1, nrow(dat2))
  }

  pred <- predict(m, type = "response")

  list(
    model = m,
    data = dat2,
    pred = pred,
    summary = tibble(
      model = prefix,
      n = nobs(m),
      AIC = AIC(m),
      BIC = BIC(m),
      AUC = auc_safe(dat2[[response]], pred),
      weighted_logloss = weighted_logloss(dat2[[response]], pred, ww)
    ),
    coef = broom::tidy(m, conf.int = TRUE) |>
      mutate(model = prefix, .before = 1),
    vif = safe_vif_table(m) |>
      mutate(model = prefix, .before = 1)
  )
}

extract_gam_smooth <- function(m, model_name) {
  st <- summary(m)$s.table
  if (is.null(st) || nrow(st) == 0) {
    return(tibble(model = model_name, term = character()))
  }

  st_df <- as.data.frame(st) |>
    tibble::rownames_to_column("term") |>
    as_tibble()

  stat_col <- intersect(c("Chi.sq", "F"), names(st_df))[1]
  p_col <- intersect(c("p-value", "p.value"), names(st_df))[1]

  st_df |>
    rename(
      statistic = all_of(stat_col),
      p.value = all_of(p_col)
    ) |>
    mutate(model = model_name, .before = 1)
}

## =========================================================
## 5. Load data
## =========================================================

stop_if_missing(ENRICHED_FILE)
stop_if_missing(NO_RANGE_VAR_FILE)

msg("Loading enriched table...")
dat <- readr::read_csv(ENRICHED_FILE, show_col_types = FALSE)

if (!"orientation_factor" %in% names(dat)) {
  if ("y_nodding" %in% names(dat)) {
    dat$orientation_factor <- factor(
      ifelse(dat$y_nodding == 1, "nodding", "upward"),
      levels = c("upward", "nodding")
    )
  } else {
    stop("No orientation_factor or y_nodding column found.")
  }
} else {
  dat$orientation_factor <- factor(dat$orientation_factor, levels = c("upward", "nodding"))
}

dat <- dat |>
  mutate(
    y_nodding01 = ifelse(orientation_factor == "nodding", 1, 0),
    z_lon = if ("z_lon" %in% names(dat)) z_lon else safe_scale(longitude),
    z_lat = if ("z_lat" %in% names(dat)) z_lat else safe_scale(latitude),
    z_lon2 = z_lon^2,
    z_lat2 = z_lat^2,
    z_lon_lat = z_lon * z_lat
  )

dat$case_weight <- if (USE_CLASS_BALANCED_WEIGHTS) make_class_weights(dat$y_nodding01) else 1

msg("Loaded rows: ", nrow(dat))
print(dat |> count(orientation_factor))

topo_vars_available <- TOPO_VARS[TOPO_VARS %in% names(dat)]
soil_vars_available <- SOIL_VARS[SOIL_VARS %in% names(dat)]

msg("Available topography variables: ", paste(topo_vars_available, collapse = ", "))
msg("Available soil variables: ", paste(soil_vars_available, collapse = ", "))

write_csv(
  tibble(
    variable = c(topo_vars_available, soil_vars_available),
    class = c(
      rep("topography", length(topo_vars_available)),
      rep("soilgrids_0_30cm", length(soil_vars_available))
    ),
    n_NA = vapply(c(topo_vars_available, soil_vars_available), function(v) sum(is.na(dat[[v]])), integer(1)),
    prop_NA = vapply(c(topo_vars_available, soil_vars_available), function(v) mean(is.na(dat[[v]])), numeric(1))
  ),
  file.path(OUT_DIR, "00_available_topography_soil_summary.csv")
)

base_vars <- readr::read_csv(NO_RANGE_VAR_FILE, show_col_types = FALSE)$variable
base_vars <- unique(as.character(base_vars))
base_vars <- base_vars[base_vars %in% names(dat)]

raw_predictors <- unique(c(base_vars, topo_vars_available, soil_vars_available))
raw_predictors <- raw_predictors[raw_predictors %in% names(dat)]

na_tbl <- tibble(
  variable = raw_predictors,
  n_NA = vapply(raw_predictors, function(v) sum(is.na(dat[[v]])), integer(1)),
  prop_NA = vapply(raw_predictors, function(v) mean(is.na(dat[[v]])), numeric(1))
) |>
  arrange(desc(prop_NA), variable)

write_csv(na_tbl, file.path(RF_DIR, "01_predictor_NA_summary.csv"))

rf_predictors <- na_tbl |>
  filter(prop_NA <= MAX_PROP_NA_PREDICTOR) |>
  pull(variable)

write_csv(
  tibble(
    variable = rf_predictors,
    class = case_when(
      variable %in% base_vars ~ "CHELSA_no_minmaxrange",
      variable %in% topo_vars_available ~ "topography",
      variable %in% soil_vars_available ~ "soilgrids_0_30cm",
      TRUE ~ "other"
    )
  ),
  file.path(RF_DIR, "02_variables_used_CHELSA_topography_soilgrids.csv")
)

msg("RF predictors after NA screen: ", length(rf_predictors))
print(
  tibble(variable = rf_predictors) |>
    mutate(
      class = case_when(
        variable %in% base_vars ~ "CHELSA_no_minmaxrange",
        variable %in% topo_vars_available ~ "topography",
        variable %in% soil_vars_available ~ "soilgrids_0_30cm",
        TRUE ~ "other"
      )
    ) |>
    count(class)
)

## =========================================================
## 6. Random forests
## =========================================================

rf_occ <- run_rf_classification(
  df = dat,
  predictors = rf_predictors,
  y_col = "orientation_factor",
  prefix = "occurrence_CHELSA_TOPO_SOIL",
  out_dir = RF_DIR,
  ntree = RF_NTREE_MAIN,
  repeats = RF_REPEATS
)

rf_occ_space <- run_rf_classification(
  df = dat,
  predictors = c(rf_predictors, SPACE_VARS),
  y_col = "orientation_factor",
  prefix = "occurrence_CHELSA_TOPO_SOIL_with_space",
  out_dir = RF_DIR,
  ntree = RF_NTREE_MAIN,
  repeats = RF_REPEATS
)

imp_rf <- readr::read_csv(
  file.path(RF_DIR, "occurrence_CHELSA_TOPO_SOIL_RF_importance_repeats_summary.csv"),
  show_col_types = FALSE
)

## =========================================================
## 7. Correlation clustering + VIF
## =========================================================

imp_primary <- imp_rf |>
  filter(variable %in% rf_predictors) |>
  arrange(desc(mean_MDA))

top_rf_vars <- imp_primary |>
  slice(seq_len(min(TOP_N_RF_FOR_MODELS, n()))) |>
  pull(variable)

forced_vars <- character()
if (FORCE_TOPOGRAPHY_IN_GLM) forced_vars <- unique(c(forced_vars, topo_vars_available))
if (FORCE_SOIL_SUMMARY_IN_GLM) forced_vars <- unique(c(forced_vars, soil_vars_available))

candidate_vars <- unique(c(top_rf_vars, forced_vars))
candidate_vars <- candidate_vars[candidate_vars %in% names(dat)]

forced_imp <- tibble(
  variable = setdiff(forced_vars, imp_primary$variable),
  mean_MDA = 0,
  sd_MDA = NA_real_,
  prop_positive_MDA = NA_real_,
  mean_MDG = NA_real_
)

imp_for_selection <- bind_rows(
  imp_primary |> select(any_of(c("variable", "mean_MDA", "sd_MDA", "prop_positive_MDA", "mean_MDG"))),
  forced_imp
)

write_csv(
  tibble(
    candidate_variable = candidate_vars,
    class = case_when(
      candidate_variable %in% base_vars ~ "CHELSA_no_minmaxrange",
      candidate_variable %in% topo_vars_available ~ "topography",
      candidate_variable %in% soil_vars_available ~ "soilgrids_0_30cm",
      TRUE ~ "other"
    ),
    forced = candidate_variable %in% forced_vars
  ),
  file.path(MODEL_DIR, "01_candidate_variables_new_RFtop_topography_soilgrids.csv")
)

rep_obj <- select_cluster_representatives(
  dat = dat,
  candidate_vars = candidate_vars,
  importance_tbl = imp_for_selection,
  threshold = COR_CLUSTER_THRESHOLD,
  max_vars = MAX_REP_VARS,
  force_keep = forced_vars
)

cluster_reps <- rep_obj$reps

cluster_tbl <- rep_obj$clusters |>
  mutate(
    selected_after_correlation_clustering = variable %in% cluster_reps,
    class = case_when(
      variable %in% base_vars ~ "CHELSA_no_minmaxrange",
      variable %in% topo_vars_available ~ "topography",
      variable %in% soil_vars_available ~ "soilgrids_0_30cm",
      TRUE ~ "other"
    )
  )

write_csv(cluster_tbl, file.path(MODEL_DIR, "02_correlation_clusters_new_RFtop_topography_soilgrids.csv"))

cor_long <- as.data.frame(as.table(rep_obj$cor)) |>
  as_tibble() |>
  rename(var1 = Var1, var2 = Var2, correlation = Freq) |>
  filter(as.character(var1) < as.character(var2)) |>
  mutate(abs_correlation = abs(correlation)) |>
  arrange(desc(abs_correlation))

write_csv(cor_long, file.path(MODEL_DIR, "02_correlation_pairs_new_RFtop_topography_soilgrids.csv"))

msg("Representatives after correlation clustering:")
print(cluster_reps)

dat_model <- dat

for (v in unique(c(candidate_vars, cluster_reps))) {
  dat_model[[paste0("z_", v)]] <- safe_scale(dat_model[[v]])
}

z_cluster_reps <- paste0("z_", cluster_reps)

vif_obj <- vif_prune(
  dat = dat_model,
  response = "y_nodding01",
  predictors = z_cluster_reps,
  weights_col = "case_weight",
  threshold = VIF_THRESHOLD
)

z_final_vars <- vif_obj$kept
final_vars <- str_remove(z_final_vars, "^z_")

write_csv(vif_obj$removed, file.path(MODEL_DIR, "03_VIF_removed_variables.csv"))
write_csv(vif_obj$final_vif, file.path(MODEL_DIR, "03_VIF_final_table.csv"))

write_csv(
  tibble(
    final_variable = final_vars,
    final_scaled_variable = z_final_vars,
    class = case_when(
      final_variable %in% base_vars ~ "CHELSA_no_minmaxrange",
      final_variable %in% topo_vars_available ~ "topography",
      final_variable %in% soil_vars_available ~ "soilgrids_0_30cm",
      TRUE ~ "other"
    )
  ),
  file.path(MODEL_DIR, "04_final_variables_for_GLM_GAM.csv")
)

msg("Final variables after VIF pruning:")
print(final_vars)

## =========================================================
## 8. Common-N weighted GLMs
## =========================================================

glm_common_cols <- unique(c("y_nodding01", "case_weight", z_final_vars, SPACE_POLY_VARS, "longitude", "latitude", "orientation_factor"))
glm_common_cols <- glm_common_cols[glm_common_cols %in% names(dat_model)]

glm_dat_common <- dat_model[, glm_common_cols, drop = FALSE] |>
  drop_na()

write_csv(
  tibble(
    n_original = nrow(dat_model),
    n_common_complete_case = nrow(glm_dat_common),
    n_removed_by_complete_case = nrow(dat_model) - nrow(glm_dat_common)
  ),
  file.path(MODEL_DIR, "05_GLM_commonN_summary.csv")
)

msg("GLM common complete-case n = ", nrow(glm_dat_common))

glm_env <- fit_glm_set(
  dat = glm_dat_common,
  response = "y_nodding01",
  vars = z_final_vars,
  prefix = "env_topography_soil_new_RFtop_VIF",
  weights_col = "case_weight"
)

glm_space <- fit_glm_set(
  dat = glm_dat_common,
  response = "y_nodding01",
  vars = SPACE_POLY_VARS,
  prefix = "space_only_commonN",
  weights_col = "case_weight"
)

glm_env_space <- fit_glm_set(
  dat = glm_dat_common,
  response = "y_nodding01",
  vars = c(z_final_vars, SPACE_POLY_VARS),
  prefix = "env_topography_soil_plus_space_new_RFtop_VIF",
  weights_col = "case_weight"
)

glm_summaries <- bind_rows(glm_env$summary, glm_space$summary, glm_env_space$summary) |>
  arrange(AIC) |>
  mutate(delta_AIC = AIC - min(AIC))

glm_coefs <- bind_rows(glm_env$coef, glm_space$coef, glm_env_space$coef)
glm_vif <- bind_rows(glm_env$vif, glm_space$vif, glm_env_space$vif)

write_csv(glm_summaries, file.path(MODEL_DIR, "05_GLM_model_comparison.csv"))
write_csv(glm_coefs, file.path(MODEL_DIR, "05_GLM_coefficients.csv"))
write_csv(glm_vif, file.path(MODEL_DIR, "05_GLM_VIF.csv"))

msg("GLM model comparison, common N:")
print(glm_summaries, n = Inf, width = Inf)

lrt_env_after_space <- anova(glm_space$model, glm_env_space$model, test = "Chisq") |>
  as.data.frame() |>
  rownames_to_column("step") |>
  as_tibble() |>
  mutate(test = "env_topography_soil_after_space_commonN")

lrt_space_after_env <- anova(glm_env$model, glm_env_space$model, test = "Chisq") |>
  as.data.frame() |>
  rownames_to_column("step") |>
  as_tibble() |>
  mutate(test = "space_after_env_topography_soil_commonN")

lrt_tbl <- bind_rows(lrt_env_after_space, lrt_space_after_env)
write_csv(lrt_tbl, file.path(MODEL_DIR, "06_GLM_incremental_tests.csv"))

msg("GLM incremental tests, common N:")
print(lrt_tbl, n = Inf, width = Inf)

## Add GLM predictions to common data
glm_pred_dat <- glm_dat_common |>
  mutate(
    pred_env = predict(glm_env$model, newdata = glm_dat_common, type = "response"),
    pred_space = predict(glm_space$model, newdata = glm_dat_common, type = "response"),
    pred_env_space = predict(glm_env_space$model, newdata = glm_dat_common, type = "response"),
    resid_env_space = y_nodding01 - pred_env_space
  )

write_csv(glm_pred_dat, file.path(MODEL_DIR, "05_GLM_commonN_predictions.csv"))

## =========================================================
## 9. GAMs
## =========================================================

gam_cols <- unique(c("y_nodding01", "case_weight", "longitude", "latitude", z_final_vars))
gam_cols <- gam_cols[gam_cols %in% names(dat_model)]

gam_dat <- dat_model[, gam_cols, drop = FALSE] |>
  drop_na()

if (nrow(gam_dat) > GAM_MAX_N) {
  set.seed(42)
  n_each <- floor(GAM_MAX_N / 2)

  gam_dat <- gam_dat |>
    group_by(y_nodding01) |>
    group_modify(function(.x, .y) {
      slice_sample(.x, n = min(n_each, nrow(.x)))
    }) |>
    ungroup()
}

k_space <- min(80, max(20, floor(sqrt(nrow(gam_dat)))))
smooth_terms <- paste0("s(", z_final_vars, ", k = ", GAM_ENV_K, ")", collapse = " + ")

gam_formula_env <- as.formula(paste("y_nodding01 ~", smooth_terms))

gam_formula_env_space <- as.formula(
  paste(
    "y_nodding01 ~",
    smooth_terms,
    "+ s(longitude, latitude, k =", k_space, ")"
  )
)

msg("Fitting GAM env/topography/soil")
gam_env <- mgcv::gam(
  gam_formula_env,
  data = gam_dat,
  family = binomial,
  weights = case_weight,
  method = "REML"
)

msg("Fitting GAM env/topography/soil + space")
gam_env_space <- mgcv::gam(
  gam_formula_env_space,
  data = gam_dat,
  family = binomial,
  weights = case_weight,
  method = "REML"
)

gam_models <- list(
  GAM_env_topography_soil = gam_env,
  GAM_env_topography_soil_plus_space = gam_env_space
)

gam_summary <- bind_rows(lapply(names(gam_models), function(nm) {
  m <- gam_models[[nm]]
  pred <- predict(m, type = "response")
  tibble(
    model = nm,
    n = nobs(m),
    AIC = AIC(m),
    BIC = BIC(m),
    AUC = auc_safe(m$y, pred),
    deviance_explained = summary(m)$dev.expl
  )
})) |>
  arrange(AIC) |>
  mutate(delta_AIC = AIC - min(AIC))

gam_smooth <- bind_rows(lapply(names(gam_models), function(nm) {
  extract_gam_smooth(gam_models[[nm]], nm)
}))

write_csv(gam_summary, file.path(MODEL_DIR, "07_GAM_model_comparison.csv"))
write_csv(gam_smooth, file.path(MODEL_DIR, "07_GAM_smooth_terms.csv"))
saveRDS(gam_env, file.path(MODEL_DIR, "07_GAM_env_topography_soil.rds"))
saveRDS(gam_env_space, file.path(MODEL_DIR, "07_GAM_env_topography_soil_plus_space.rds"))

msg("GAM comparison:")
print(gam_summary, n = Inf, width = Inf)

## =========================================================
## 10. Publication figures
## =========================================================

fig_index <- list()

## FigA: RF importance
rf_imp_plot <- imp_rf |>
  slice_max(mean_MDA, n = 25, with_ties = FALSE) |>
  mutate(
    variable_label = pretty_var(variable),
    variable_label = fct_reorder(variable_label, mean_MDA),
    class = class_of_var(variable)
  )

save_data(rf_imp_plot, "FigA_RF_importance_top25_data")

p_rf <- ggplot(rf_imp_plot, aes(x = mean_MDA, y = variable_label)) +
  geom_col(width = 0.72) +
  geom_ci_horizontal(aes(xmin = mean_MDA - sd_MDA, xmax = mean_MDA + sd_MDA), height = 0.15, linewidth = 0.35) +
  theme_paper(10) +
  paper_labs(
    x = "Random forest importance\n(mean decrease accuracy)",
    y = NULL,
    title = "Random forest variable importance"
  )

fig_index[[length(fig_index) + 1]] <- save_pub(p_rf, "FigA_RF_importance_top25", width = 7.4, height = 7.2)

## FigA2: RF with space importance
imp_rf_space <- read_csv(file.path(RF_DIR, "occurrence_CHELSA_TOPO_SOIL_with_space_RF_importance_repeats_summary.csv"), show_col_types = FALSE)

rf_imp_space_plot <- imp_rf_space |>
  slice_max(mean_MDA, n = 25, with_ties = FALSE) |>
  mutate(
    variable_label = pretty_var(variable),
    variable_label = fct_reorder(variable_label, mean_MDA),
    class = class_of_var(variable)
  )

save_data(rf_imp_space_plot, "FigA2_RF_with_space_importance_top25_data")

p_rf_space <- ggplot(rf_imp_space_plot, aes(x = mean_MDA, y = variable_label)) +
  geom_col(width = 0.72) +
  geom_ci_horizontal(aes(xmin = mean_MDA - sd_MDA, xmax = mean_MDA + sd_MDA), height = 0.15, linewidth = 0.35) +
  theme_paper(10) +
  paper_labs(
    x = "Random forest importance\n(mean decrease accuracy)",
    y = NULL,
    title = "Random forest variable importance including space"
  )

fig_index[[length(fig_index) + 1]] <- save_pub(p_rf_space, "FigA2_RF_with_space_importance_top25", width = 7.4, height = 7.2)

## FigB: GLM comparison
glm_comp_plot <- glm_summaries |>
  arrange(AIC) |>
  mutate(
    model_label = case_when(
      model == "env_topography_soil_new_RFtop_VIF" ~ "Environment + topography + soil",
      model == "space_only_commonN" ~ "Space only",
      model == "env_topography_soil_plus_space_new_RFtop_VIF" ~ "Environment + topography + soil + space",
      TRUE ~ model
    ),
    model_label = fct_reorder(model_label, delta_AIC)
  )

save_data(glm_comp_plot, "FigB_GLM_model_comparison_data")

p_glm_comp <- ggplot(glm_comp_plot, aes(x = delta_AIC, y = model_label)) +
  geom_col(width = 0.72) +
  theme_paper(11) +
  paper_labs(
    x = expression(Delta*AIC),
    y = NULL,
    title = "Weighted GLM model comparison"
  )

fig_index[[length(fig_index) + 1]] <- save_pub(p_glm_comp, "FigB_GLM_model_comparison_deltaAIC", width = 6.4, height = 3.4)

## FigC: GLM coefficients
coef_src <- glm_coefs |>
  filter(model == "env_topography_soil_plus_space_new_RFtop_VIF", term != "(Intercept)")

if (!"conf.low" %in% names(coef_src)) {
  coef_src$conf.low <- coef_src$estimate - 1.96 * coef_src$std.error
}
if (!"conf.high" %in% names(coef_src)) {
  coef_src$conf.high <- coef_src$estimate + 1.96 * coef_src$std.error
}

coef_plot <- coef_src |>
  mutate(
    term_label = pretty_var(term),
    term_label = fct_reorder(term_label, estimate),
    class = class_of_var(term)
  )

save_data(coef_plot, "FigC_GLM_coefficients_data")

p_coef <- ggplot(coef_plot, aes(x = estimate, y = term_label)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.35) +
  geom_ci_horizontal(aes(xmin = conf.low, xmax = conf.high), height = 0.15, linewidth = 0.42) +
  geom_point(size = 2.2) +
  theme_paper(10.5) +
  paper_labs(
    x = "Weighted logistic coefficient",
    y = NULL,
    title = "GLM coefficients"
  )

fig_index[[length(fig_index) + 1]] <- save_pub(p_coef, "FigC_GLM_coefficients_env_topo_soil_space", width = 7.2, height = 6.4)

## FigD1: top correlated pairs
cor_top <- cor_long |>
  filter(is.finite(abs_correlation)) |>
  slice_max(abs_correlation, n = 30, with_ties = FALSE) |>
  mutate(
    pair = paste(pretty_var(as.character(var1)), pretty_var(as.character(var2)), sep = " — "),
    pair = fct_reorder(pair, abs_correlation)
  )

save_data(cor_top, "FigD1_top_correlated_pairs_data")

p_cor_pairs <- ggplot(cor_top, aes(x = abs_correlation, y = pair)) +
  geom_col(width = 0.72) +
  geom_vline(xintercept = COR_CLUSTER_THRESHOLD, linetype = 2, linewidth = 0.35) +
  theme_paper(8.5) +
  paper_labs(
    x = "|Pearson correlation|",
    y = NULL,
    title = "Top correlated candidate predictors"
  )

fig_index[[length(fig_index) + 1]] <- save_pub(p_cor_pairs, "FigD1_top_correlated_pairs", width = 7.6, height = 7.0)

## FigD2: final VIF
vif_plot <- vif_obj$final_vif |>
  filter(!is.na(variable)) |>
  mutate(
    variable_label = pretty_var(variable),
    variable_label = fct_reorder(variable_label, VIF)
  )

save_data(vif_plot, "FigD2_final_VIF_data")

p_vif <- ggplot(vif_plot, aes(x = VIF, y = variable_label)) +
  geom_col(width = 0.72) +
  geom_vline(xintercept = VIF_THRESHOLD, linetype = 2, linewidth = 0.35) +
  theme_paper(10) +
  paper_labs(
    x = "Variance inflation factor",
    y = NULL,
    title = "VIF after pruning"
  )

fig_index[[length(fig_index) + 1]] <- save_pub(p_vif, "FigD2_final_VIF_after_pruning", width = 5.8, height = 4.8)

## FigE: GLM response curves
response_curve_data <- list()

for (v in z_final_vars) {
  newdat <- glm_dat_common[rep(1, 140), , drop = FALSE]

  for (vv in c(z_final_vars, SPACE_POLY_VARS)) {
    if (vv %in% names(newdat)) newdat[[vv]] <- median(glm_dat_common[[vv]], na.rm = TRUE)
  }

  grid <- seq(
    quantile(glm_dat_common[[v]], 0.02, na.rm = TRUE),
    quantile(glm_dat_common[[v]], 0.98, na.rm = TRUE),
    length.out = 140
  )

  newdat[[v]] <- grid

  pr <- predict(glm_env_space$model, newdata = newdat, type = "link", se.fit = TRUE)

  dd <- tibble(
    variable = v,
    variable_label = pretty_var(v),
    x = grid,
    fit_link = pr$fit,
    se_link = pr$se.fit,
    fit = plogis(pr$fit),
    lo = plogis(pr$fit - 1.96 * pr$se.fit),
    hi = plogis(pr$fit + 1.96 * pr$se.fit)
  )

  response_curve_data[[v]] <- dd
}

rc <- bind_rows(response_curve_data)
save_data(rc, "FigE_GLM_response_curves_data")

p_curves <- rc |>
  mutate(variable_label = fct_reorder(variable_label, fit, .fun = mean)) |>
  ggplot(aes(x = x, y = fit)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.18) +
  geom_line(linewidth = 0.75) +
  facet_wrap(~ variable_label, scales = "free_x", ncol = 3) +
  theme_paper(9) +
  theme(panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.22)) +
  paper_labs(
    x = "Standardized predictor",
    y = "Predicted probability of nodding",
    title = "GLM response curves"
  )

fig_index[[length(fig_index) + 1]] <- save_pub(p_curves, "FigE_GLM_response_curves", width = 8.8, height = 7.8)

## FigF: GLM prediction map
save_data(glm_pred_dat, "FigF_GLM_prediction_map_data")

p_map <- ggplot(glm_pred_dat, aes(x = longitude, y = latitude, color = pred_env_space)) +
  geom_point(size = 0.55, alpha = 0.60) +
  coord_equal() +
  scale_color_viridis_c(option = "C", limits = c(0, 1)) +
  theme_paper(10) +
  theme(
    panel.grid.major = element_line(colour = "grey90", linewidth = 0.2),
    panel.grid.minor = element_blank()
  ) +
  paper_labs(
    x = "Longitude",
    y = "Latitude",
    color = "Predicted\nprobability",
    title = "Predicted distribution of nodding head orientation"
  )

fig_index[[length(fig_index) + 1]] <- save_pub(p_map, "FigF_GLM_prediction_map", width = 7.2, height = 4.8)

## FigF2: GLM residual map
p_resid_map <- ggplot(glm_pred_dat, aes(x = longitude, y = latitude, color = resid_env_space)) +
  geom_point(size = 0.55, alpha = 0.60) +
  coord_equal() +
  scale_color_gradient2(midpoint = 0) +
  theme_paper(10) +
  theme(
    panel.grid.major = element_line(colour = "grey90", linewidth = 0.2),
    panel.grid.minor = element_blank()
  ) +
  paper_labs(
    x = "Longitude",
    y = "Latitude",
    color = "Observed -\npredicted",
    title = "GLM residual map"
  )

fig_index[[length(fig_index) + 1]] <- save_pub(p_resid_map, "FigF2_GLM_residual_map", width = 7.2, height = 4.8)

## FigG: GAM comparison
gam_comp_plot <- gam_summary |>
  mutate(
    model_label = case_when(
      model == "GAM_env_topography_soil" ~ "GAM: environment",
      model == "GAM_env_topography_soil_plus_space" ~ "GAM: environment + space",
      TRUE ~ model
    ),
    model_label = fct_reorder(model_label, delta_AIC)
  )

save_data(gam_comp_plot, "FigG_GAM_model_comparison_data")

p_gam_comp <- ggplot(gam_comp_plot, aes(x = delta_AIC, y = model_label)) +
  geom_col(width = 0.72) +
  theme_paper(11) +
  paper_labs(
    x = expression(Delta*AIC),
    y = NULL,
    title = "GAM model comparison"
  )

fig_index[[length(fig_index) + 1]] <- save_pub(p_gam_comp, "FigG_GAM_model_comparison_deltaAIC", width = 6.2, height = 3.0)

## FigG2: GAM smooth term statistics
gam_smooth_plot <- gam_smooth |>
  filter(model == "GAM_env_topography_soil_plus_space") |>
  mutate(
    term_label = term |>
      str_replace_all("^s\\(", "") |>
      str_replace_all("\\)$", "") |>
      str_replace_all(",.*$", "") |>
      pretty_var(),
    term_label = fct_reorder(term_label, statistic)
  )

if (nrow(gam_smooth_plot) > 0) {
  save_data(gam_smooth_plot, "FigG2_GAM_smooth_terms_data")

  p_gam_terms <- ggplot(gam_smooth_plot, aes(x = statistic, y = term_label)) +
    geom_col(width = 0.72) +
    theme_paper(10) +
    paper_labs(
      x = "Smooth-term test statistic",
      y = NULL,
      title = "GAM smooth-term importance"
    )

  fig_index[[length(fig_index) + 1]] <- save_pub(p_gam_terms, "FigG2_GAM_smooth_term_statistics", width = 6.4, height = 5.2)
}

## FigH: GAM smooths using mgcv base plot
png(file.path(PUB_PNG_DIR, "FigH_GAM_smooths_mgcv.png"), width = 2400, height = 1900, res = 260, bg = "white")
par(mfrow = c(4, 4), mar = c(4, 4, 2, 1), cex = 0.75)
plot(gam_env_space, pages = 1, shade = TRUE, residuals = FALSE, seWithMean = TRUE)
dev.off()

pdf(file.path(PUB_PDF_DIR, "FigH_GAM_smooths_mgcv.pdf"), width = 11, height = 8.5)
par(mfrow = c(4, 4), mar = c(4, 4, 2, 1), cex = 0.75)
plot(gam_env_space, pages = 1, shade = TRUE, residuals = FALSE, seWithMean = TRUE)
dev.off()

tiff(file.path(PUB_TIFF_DIR, "FigH_GAM_smooths_mgcv.tiff"), width = 11, height = 8.5, units = "in", res = DPI, compression = "lzw", bg = "white")
par(mfrow = c(4, 4), mar = c(4, 4, 2, 1), cex = 0.75)
plot(gam_env_space, pages = 1, shade = TRUE, residuals = FALSE, seWithMean = TRUE)
dev.off()

fig_index[[length(fig_index) + 1]] <- tibble(
  figure = "FigH_GAM_smooths_mgcv",
  png = file.path(PUB_PNG_DIR, "FigH_GAM_smooths_mgcv.png"),
  pdf = file.path(PUB_PDF_DIR, "FigH_GAM_smooths_mgcv.pdf"),
  tiff = file.path(PUB_TIFF_DIR, "FigH_GAM_smooths_mgcv.tiff"),
  exists_png = file.exists(file.path(PUB_PNG_DIR, "FigH_GAM_smooths_mgcv.png")),
  exists_pdf = file.exists(file.path(PUB_PDF_DIR, "FigH_GAM_smooths_mgcv.pdf")),
  exists_tiff = file.exists(file.path(PUB_TIFF_DIR, "FigH_GAM_smooths_mgcv.tiff"))
)

## FigI: orientation occurrence map
p_obs <- ggplot(glm_pred_dat, aes(x = longitude, y = latitude, shape = orientation_factor)) +
  geom_point(size = 0.60, alpha = 0.45) +
  coord_equal() +
  theme_paper(10) +
  theme(
    panel.grid.major = element_line(colour = "grey90", linewidth = 0.2),
    panel.grid.minor = element_blank()
  ) +
  paper_labs(
    x = "Longitude",
    y = "Latitude",
    shape = "Orientation",
    title = "Observed/predicted head orientation records"
  )

fig_index[[length(fig_index) + 1]] <- save_pub(p_obs, "FigI_orientation_occurrence_map", width = 7.2, height = 4.8)

## FigJ: raw final predictor boxplots
box_dat <- dat_model |>
  select(any_of(c("orientation_factor", final_vars))) |>
  drop_na(orientation_factor) |>
  pivot_longer(-orientation_factor, names_to = "variable", values_to = "value") |>
  mutate(variable_label = pretty_var(variable))

save_data(box_dat, "FigJ_raw_final_predictor_boxplots_data")

p_box <- ggplot(box_dat, aes(x = orientation_factor, y = value)) +
  geom_boxplot(outlier.alpha = 0.06, linewidth = 0.30) +
  facet_wrap(~ variable_label, scales = "free_y", ncol = 3) +
  theme_paper(8.5) +
  theme(panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.22)) +
  paper_labs(
    x = NULL,
    y = "Raw predictor value",
    title = "Raw environmental differences between orientation classes"
  )

fig_index[[length(fig_index) + 1]] <- save_pub(p_box, "FigJ_raw_final_predictor_boxplots", width = 8.8, height = 7.6)

## Combined main figures
p_main1 <- (p_rf | p_glm_comp) / p_coef +
  plot_layout(heights = c(1, 1.15)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 12))

fig_index[[length(fig_index) + 1]] <- save_pub(p_main1, "MainFig1_RF_GLM_model_selection_coefficients", width = 12, height = 10)

p_main2 <- p_curves / p_map +
  plot_layout(heights = c(1.4, 0.8)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 12))

fig_index[[length(fig_index) + 1]] <- save_pub(p_main2, "MainFig2_response_curves_prediction_map", width = 11, height = 12)

if (exists("p_gam_terms")) {
  p_main3 <- p_gam_comp | p_gam_terms +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(face = "bold", size = 12))

  fig_index[[length(fig_index) + 1]] <- save_pub(p_main3, "MainFig3_GAM_summary", width = 10.5, height = 5.5)
}

## =========================================================
## 11. Diagnostics and README
## =========================================================

fig_index_tbl <- bind_rows(fig_index)
write_csv(fig_index_tbl, file.path(PUB_DIR, "00_publication_figure_index.csv"))

summary_tbl <- tibble(
  output_dir = OUT_DIR,
  publication_figure_dir = PUB_DIR,
  n_rows_loaded = nrow(dat),
  n_rf_predictors = length(rf_predictors),
  n_candidate_variables = length(candidate_vars),
  n_final_variables = length(final_vars),
  best_glm_model = glm_summaries$model[which.min(glm_summaries$AIC)],
  best_glm_AIC = min(glm_summaries$AIC, na.rm = TRUE),
  best_glm_AUC = glm_summaries$AUC[which.min(glm_summaries$AIC)],
  best_gam_model = gam_summary$model[which.min(gam_summary$AIC)],
  best_gam_AIC = min(gam_summary$AIC, na.rm = TRUE),
  best_gam_deviance_explained = gam_summary$deviance_explained[which.min(gam_summary$AIC)],
  RF_NTREE_MAIN = RF_NTREE_MAIN,
  RF_REPEATS = RF_REPEATS,
  COR_CLUSTER_THRESHOLD = COR_CLUSTER_THRESHOLD,
  VIF_THRESHOLD = VIF_THRESHOLD
)

write_csv(summary_tbl, file.path(OUT_DIR, "00_analysis_summary.csv"))
write_csv(summary_tbl, file.path(PUB_DIR, "00_publication_figure_summary.csv"))

summary_lines <- c(
  "Global Cirsium full rerun: RF + GLM/GAM + publication figures",
  "==============================================================",
  paste0("Input enriched table: ", ENRICHED_FILE),
  paste0("No-min/max/range variable file: ", NO_RANGE_VAR_FILE),
  paste0("Output dir: ", OUT_DIR),
  "",
  "Analysis:",
  paste0("- RF predictors after NA screen: ", length(rf_predictors)),
  paste0("- RF ntree main: ", RF_NTREE_MAIN),
  paste0("- RF repeats: ", RF_REPEATS),
  paste0("- RF top variables for candidate set: ", TOP_N_RF_FOR_MODELS),
  paste0("- Correlation clustering threshold: |r| >= ", COR_CLUSTER_THRESHOLD),
  paste0("- VIF threshold: ", VIF_THRESHOLD),
  paste0("- Final variables after VIF pruning: ", paste(final_vars, collapse = ", ")),
  "",
  "Key output folders:",
  paste0("- RF: ", RF_DIR),
  paste0("- GLM/GAM: ", MODEL_DIR),
  paste0("- Publication figures: ", PUB_DIR),
  "",
  "Publication figure formats:",
  "- PNG: 04_publication_figures/png",
  "- PDF: 04_publication_figures/pdf",
  "- TIFF: 04_publication_figures/tiff",
  "- Figure data: 04_publication_figures/figure_data",
  "",
  "Suggested main figure panels:",
  "- MainFig1_RF_GLM_model_selection_coefficients",
  "- MainFig2_response_curves_prediction_map",
  "- MainFig3_GAM_summary"
)

writeLines(summary_lines, file.path(OUT_DIR, "README_outputs.txt"))
writeLines(summary_lines, file.path(PUB_DIR, "README_publication_figures.txt"))

msg("DONE — full rerun + publication-level figures")
msg("Output dir: ", OUT_DIR)
msg("Publication figure dir: ", PUB_DIR)
msg("Final variables:")
print(final_vars)
print(fig_index_tbl, n = Inf, width = Inf)
