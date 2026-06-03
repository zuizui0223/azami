############################################################
## Resume RF only from extracted CHELSA table
##
## Use this when extraction already finished but RF crashed with:
##   all_of(predictors): elements swb, gdd10, cmimax... don't exist
##
## It reads:
##   global_cirsium_ALL_CHELSA_random_forest_reuse_existing/
##     05_orientation_points_with_all_CHELSA.csv
##     05_CHELSA_variables_kept_for_RF.csv
##
## Then runs:
##   occurrence RF
##   occurrence RF + lon/lat
##   species RF
##   species RF + lon/lat
##
## It uses base data.frame indexing instead of dplyr::select(all_of()),
## so it is much harder to crash from select/name conflicts.
############################################################

pkgs <- c(
  "dplyr", "readr", "stringr", "tibble", "tidyr",
  "ggplot2", "randomForest", "pROC", "forcats", "patchwork"
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
  library(randomForest)
  library(pROC)
  library(forcats)
  library(patchwork)
})

## avoid select conflicts
select <- dplyr::select
filter <- dplyr::filter
mutate <- dplyr::mutate
summarise <- dplyr::summarise
arrange <- dplyr::arrange
group_by <- dplyr::group_by
ungroup <- dplyr::ungroup
slice_sample <- dplyr::slice_sample
all_of <- dplyr::all_of
bind_rows <- dplyr::bind_rows

set.seed(42)

ROOT_DIR <- "C:/Users/zuizui/cirsium_inat"

OUT_DIR <- file.path(
  ROOT_DIR,
  "global_cirsium_ALL_CHELSA_random_forest_reuse_existing"
)

FIG_DIR <- file.path(OUT_DIR, "figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

DATA_FILE <- file.path(OUT_DIR, "05_orientation_points_with_all_CHELSA.csv")
VAR_FILE  <- file.path(OUT_DIR, "05_CHELSA_variables_kept_for_RF.csv")

if (!file.exists(DATA_FILE)) {
  stop("Cannot find extracted table: ", DATA_FILE,
       "\nRun the extraction script first, or check output folder name.")
}
if (!file.exists(VAR_FILE)) {
  stop("Cannot find variable list: ", VAR_FILE)
}

MAX_OCC_RF_ROWS <- 50000
MIN_SPECIES_POINTS <- 3
RF_NTREE <- 1500
RF_REPEATS <- 20
RF_MTRY <- NULL

## TRUE caused upward to be downsampled to the nodding count.
## This is usually good for imbalanced classification.
BALANCE_CLASSES_FOR_RF <- TRUE

msg <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", ..., "\n")
}

auc_safe <- function(y_factor, prob_nodding) {
  tryCatch({
    yy <- ifelse(y_factor == "nodding", 1, 0)
    as.numeric(pROC::auc(yy, prob_nodding, quiet = TRUE))
  }, error = function(e) NA_real_)
}

balance_binary_df <- function(df, y_col = "orientation_factor") {
  tab <- table(df[[y_col]])
  if (length(tab) < 2) return(df)
  n_min <- min(tab)
  df |>
    dplyr::group_by(.data[[y_col]]) |>
    dplyr::slice_sample(n = n_min) |>
    dplyr::ungroup()
}

make_safe_x <- function(df, predictors) {
  ## Critical: use base R indexing and intersect at the last second.
  ## This avoids all_of()/select() failures.
  predictors <- unique(as.character(predictors))
  predictors <- predictors[predictors %in% names(df)]

  if (length(predictors) == 0) {
    stop("No predictors remain after matching to data columns.")
  }

  x <- as.data.frame(df[, predictors, drop = FALSE])

  ## randomForest handles syntactic names best
  original_names <- names(x)
  safe_names <- make.names(original_names, unique = TRUE)
  names(x) <- safe_names

  list(
    x = x,
    predictors_original = original_names,
    predictors_safe = safe_names
  )
}


choose_mtry <- function(p) {
  ## randomForest x/y interface can choke when mtry = NULL.
  ## For classification, default is floor(sqrt(p)).
  if (is.null(RF_MTRY) || length(RF_MTRY) == 0 || is.na(RF_MTRY)) {
    return(max(1L, floor(sqrt(p))))
  }
  m <- as.integer(RF_MTRY)
  if (!is.finite(m) || is.na(m)) {
    return(max(1L, floor(sqrt(p))))
  }
  m <- max(1L, min(m, p))
  m
}

run_rf_classification <- function(df, predictors, y_col, prefix, out_dir, ntree = 1500, repeats = 20) {
  predictors <- unique(as.character(predictors))
  predictors <- predictors[predictors %in% names(df)]

  missing0 <- setdiff(predictors, names(df))
  if (length(missing0) > 0) {
    warning(prefix, ": predictors missing before RF: ", paste(missing0, collapse = ", "))
  }

  keep_cols <- unique(c(y_col, predictors))
  rf_df <- as.data.frame(df[, keep_cols, drop = FALSE])
  rf_df <- tidyr::drop_na(rf_df)

  rf_df[[y_col]] <- factor(rf_df[[y_col]], levels = c("upward", "nodding"))

  if (length(unique(rf_df[[y_col]])) < 2) {
    warning("RF skipped for ", prefix, ": only one class.")
    return(NULL)
  }

  if (nrow(rf_df) > MAX_OCC_RF_ROWS && stringr::str_detect(prefix, "^occurrence")) {
    rf_df <- rf_df |>
      dplyr::group_by(.data[[y_col]]) |>
      dplyr::slice_sample(n = min(MAX_OCC_RF_ROWS / 2, dplyr::n()), replace = FALSE) |>
      dplyr::ungroup()
  }

  if (BALANCE_CLASSES_FOR_RF) {
    rf_train <- balance_binary_df(rf_df, y_col)
  } else {
    rf_train <- rf_df
  }

  xx <- make_safe_x(rf_train, predictors)
  x_train <- xx$x
  y_train <- rf_train[[y_col]]

  msg("RF ", prefix, ": original n=", nrow(rf_df), " used n=", nrow(rf_train),
      " predictors=", ncol(x_train))
  msg("Class balancing for RF: ", BALANCE_CLASSES_FOR_RF)
  msg("mtry used: ", choose_mtry(ncol(x_train)))
  print(table(y_train))

  rf <- randomForest::randomForest(
    x = x_train,
    y = y_train,
    ntree = ntree,
    mtry = choose_mtry(ncol(x_train)),
    importance = TRUE,
    proximity = FALSE
  )

  pred_prob <- predict(rf, x_train, type = "prob")[, "nodding"]
  pred_class <- predict(rf, x_train, type = "response")

  conf <- as.data.frame.matrix(table(observed = y_train, predicted = pred_class)) |>
    tibble::rownames_to_column("observed") |>
    tibble::as_tibble()

  metrics <- tibble::tibble(
    analysis = prefix,
    n_original_before_RF = nrow(rf_df),
    n_used_for_RF = nrow(rf_train),
    class_balanced = BALANCE_CLASSES_FOR_RF,
    n_predictors = ncol(x_train),
    OOB_error = rf$err.rate[nrow(rf$err.rate), "OOB"],
    AUC_train = auc_safe(y_train, pred_prob)
  )

  imp <- as.data.frame(randomForest::importance(rf)) |>
    tibble::rownames_to_column("variable_safe") |>
    tibble::as_tibble() |>
    dplyr::mutate(
      variable = xx$predictors_original[match(variable_safe, xx$predictors_safe)],
      .after = variable_safe
    ) |>
    dplyr::arrange(dplyr::desc(MeanDecreaseAccuracy))

  readr::write_csv(metrics, file.path(out_dir, paste0(prefix, "_RF_metrics.csv")))
  readr::write_csv(conf, file.path(out_dir, paste0(prefix, "_RF_confusion_matrix.csv")))
  readr::write_csv(imp, file.path(out_dir, paste0(prefix, "_RF_importance.csv")))
  saveRDS(rf, file.path(out_dir, paste0(prefix, "_RF_model.rds")))

  rep_list <- list()

  for (i in seq_len(repeats)) {
    set.seed(1000 + i)

    if (BALANCE_CLASSES_FOR_RF) {
      train_i <- balance_binary_df(rf_df, y_col)
    } else {
      train_i <- rf_df
    }

    if (nrow(train_i) > MAX_OCC_RF_ROWS && stringr::str_detect(prefix, "^occurrence")) {
      train_i <- train_i |>
        dplyr::group_by(.data[[y_col]]) |>
        dplyr::slice_sample(n = min(MAX_OCC_RF_ROWS / 2, dplyr::n()), replace = FALSE) |>
        dplyr::ungroup()
    }

    xxi <- make_safe_x(train_i, predictors)
    x_i <- xxi$x
    y_i <- train_i[[y_col]]

    fit_i <- randomForest::randomForest(
      x = x_i,
      y = y_i,
      ntree = 800,
      mtry = choose_mtry(ncol(x_i)),
      importance = TRUE,
      proximity = FALSE
    )

    imp_i <- as.data.frame(randomForest::importance(fit_i)) |>
      tibble::rownames_to_column("variable_safe") |>
      tibble::as_tibble() |>
      dplyr::mutate(
        variable = xxi$predictors_original[match(variable_safe, xxi$predictors_safe)],
        repeat_id = i,
        .before = 1
      )

    rep_list[[i]] <- imp_i
  }

  imp_rep <- dplyr::bind_rows(rep_list)

  imp_summary <- imp_rep |>
    dplyr::group_by(variable) |>
    dplyr::summarise(
      mean_MDA = mean(MeanDecreaseAccuracy, na.rm = TRUE),
      sd_MDA = sd(MeanDecreaseAccuracy, na.rm = TRUE),
      prop_positive_MDA = mean(MeanDecreaseAccuracy > 0, na.rm = TRUE),
      mean_MDG = mean(MeanDecreaseGini, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::arrange(dplyr::desc(mean_MDA))

  readr::write_csv(imp_rep, file.path(out_dir, paste0(prefix, "_RF_importance_repeats_raw.csv")))
  readr::write_csv(imp_summary, file.path(out_dir, paste0(prefix, "_RF_importance_repeats_summary.csv")))

  p_imp <- imp_summary |>
    dplyr::slice_max(mean_MDA, n = 30) |>
    dplyr::mutate(variable = forcats::fct_reorder(variable, mean_MDA)) |>
    ggplot2::ggplot(ggplot2::aes(mean_MDA, variable)) +
    ggplot2::geom_col() +
    ggplot2::geom_errorbarh(
      ggplot2::aes(xmin = mean_MDA - sd_MDA, xmax = mean_MDA + sd_MDA),
      height = 0.15
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::labs(
      x = "Mean decrease accuracy",
      y = NULL,
      title = paste0(prefix, " RF variable importance")
    )

  ggplot2::ggsave(file.path(FIG_DIR, paste0(prefix, "_RF_importance_top30.png")),
                  p_imp, width = 8, height = 8, dpi = 300)

  ## Partial dependence-like median curves for top 8 predictors
  top_vars <- head(imp_summary$variable, 8)
  pdp_list <- list()

  for (v in top_vars) {
    if (is.na(v) || !v %in% xx$predictors_original) next

    v_safe <- xx$predictors_safe[match(v, xx$predictors_original)]
    if (is.na(v_safe) || !v_safe %in% names(x_train)) next

    grid <- seq(
      stats::quantile(x_train[[v_safe]], 0.02, na.rm = TRUE),
      stats::quantile(x_train[[v_safe]], 0.98, na.rm = TRUE),
      length.out = 50
    )

    base <- x_train[rep(1, length(grid)), , drop = FALSE]

    for (pv in names(base)) {
      base[[pv]] <- stats::median(x_train[[pv]], na.rm = TRUE)
    }

    base[[v_safe]] <- grid
    pp <- predict(rf, base, type = "prob")[, "nodding"]

    pd <- tibble::tibble(variable = v, x = grid, pred_prob_nodding = pp)
    readr::write_csv(pd, file.path(out_dir, paste0(prefix, "_PDP_", v, ".csv")))

    pdp_list[[v]] <- ggplot2::ggplot(pd, ggplot2::aes(x, pred_prob_nodding)) +
      ggplot2::geom_line(linewidth = 1) +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::labs(
        x = v,
        y = "Predicted P(nodding)",
        title = v
      )
  }

  if (length(pdp_list) > 0) {
    p_pdp <- patchwork::wrap_plots(pdp_list, ncol = 2)
    ggplot2::ggsave(file.path(FIG_DIR, paste0(prefix, "_RF_partial_dependence_top8.png")),
                    p_pdp, width = 10, height = 10, dpi = 300)
  }

  list(model = rf, metrics = metrics, importance = imp_summary)
}

## =========================================================
## Load data and variable list
## =========================================================

dat_env_cc <- readr::read_csv(DATA_FILE, show_col_types = FALSE)
var_tbl <- readr::read_csv(VAR_FILE, show_col_types = FALSE)

if (!"variable" %in% names(var_tbl)) {
  stop("Variable file must have a column named 'variable'.")
}

env_vars_keep <- unique(as.character(var_tbl$variable))
env_vars_keep <- env_vars_keep[env_vars_keep %in% names(dat_env_cc)]

if (!"orientation_factor" %in% names(dat_env_cc)) {
  if ("y_nodding" %in% names(dat_env_cc)) {
    dat_env_cc$orientation_factor <- factor(
      ifelse(dat_env_cc$y_nodding == 1, "nodding", "upward"),
      levels = c("upward", "nodding")
    )
  } else {
    stop("No orientation_factor or y_nodding column found.")
  }
} else {
  dat_env_cc$orientation_factor <- factor(dat_env_cc$orientation_factor, levels = c("upward", "nodding"))
}

## Ensure space vars exist
if (!"z_lon" %in% names(dat_env_cc)) {
  dat_env_cc$z_lon <- as.numeric(scale(dat_env_cc$longitude))
}
if (!"z_lat" %in% names(dat_env_cc)) {
  dat_env_cc$z_lat <- as.numeric(scale(dat_env_cc$latitude))
}

msg("Loaded extracted CHELSA table:")
msg("n rows: ", nrow(dat_env_cc))
print(dat_env_cc |> dplyr::count(orientation_factor))
msg("n env predictors found in data: ", length(env_vars_keep))

readr::write_csv(tibble::tibble(variable = env_vars_keep),
                 file.path(OUT_DIR, "RF_resume_variables_used.csv"))

## =========================================================
## RF occurrence-level
## =========================================================

rf_occ <- run_rf_classification(
  df = dat_env_cc,
  predictors = env_vars_keep,
  y_col = "orientation_factor",
  prefix = "occurrence_all_CHELSA",
  out_dir = OUT_DIR,
  ntree = RF_NTREE,
  repeats = RF_REPEATS
)

rf_occ_space <- run_rf_classification(
  df = dat_env_cc,
  predictors = c(env_vars_keep, "z_lon", "z_lat"),
  y_col = "orientation_factor",
  prefix = "occurrence_all_CHELSA_with_space",
  out_dir = OUT_DIR,
  ntree = RF_NTREE,
  repeats = RF_REPEATS
)

## =========================================================
## Species-level RF
## =========================================================

if (!"species_name" %in% names(dat_env_cc)) {
  warning("No species_name column; species-level RF skipped.")
} else {
  species_dat <- dat_env_cc |>
    dplyr::filter(!is.na(species_name), species_name != "") |>
    dplyr::group_by(species_name) |>
    dplyr::summarise(
      n_points = dplyr::n(),
      prop_nodding = mean(ifelse(orientation_factor == "nodding", 1, 0), na.rm = TRUE),
      orientation_factor = factor(
        ifelse(prop_nodding >= 0.5, "nodding", "upward"),
        levels = c("upward", "nodding")
      ),
      mean_lon = mean(longitude, na.rm = TRUE),
      mean_lat = mean(latitude, na.rm = TRUE),
      dplyr::across(dplyr::all_of(env_vars_keep), ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    ) |>
    dplyr::filter(n_points >= MIN_SPECIES_POINTS) |>
    dplyr::mutate(
      z_lon = as.numeric(scale(mean_lon)),
      z_lat = as.numeric(scale(mean_lat))
    )

  readr::write_csv(species_dat, file.path(OUT_DIR, "06_species_level_all_CHELSA_summary.csv"))

  msg("Species-level counts:")
  print(species_dat |> dplyr::count(orientation_factor))
  msg("n species used: ", nrow(species_dat))

  if (nrow(species_dat) >= 20 && length(unique(species_dat$orientation_factor)) == 2) {
    rf_species <- run_rf_classification(
      df = species_dat,
      predictors = env_vars_keep,
      y_col = "orientation_factor",
      prefix = "species_all_CHELSA",
      out_dir = OUT_DIR,
      ntree = RF_NTREE,
      repeats = RF_REPEATS
    )

    rf_species_space <- run_rf_classification(
      df = species_dat,
      predictors = c(env_vars_keep, "z_lon", "z_lat"),
      y_col = "orientation_factor",
      prefix = "species_all_CHELSA_with_space",
      out_dir = OUT_DIR,
      ntree = RF_NTREE,
      repeats = RF_REPEATS
    )
  } else {
    warning("Species-level RF skipped: not enough species or only one class.")
  }
}

msg("DONE — RF resume from extracted CHELSA table.")
msg("Key outputs:")
msg("  occurrence_all_CHELSA_RF_importance_repeats_summary.csv")
msg("  occurrence_all_CHELSA_with_space_RF_importance_repeats_summary.csv")
msg("  species_all_CHELSA_RF_importance_repeats_summary.csv")
msg("  species_all_CHELSA_with_space_RF_importance_repeats_summary.csv")
