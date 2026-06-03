############################################################
## FULL SCRIPT
## Cirsium head orientation:
##   Pre-VIF candidate environmental PCA
##   + reduced pollinator SDM predictors
##   + multicollinearity-safe GLM
##
## Key design:
##   1. PCA uses PRE-VIF candidate environmental variables, not VIF-selected stack.
##      Input:
##        chelsa_pollinator_enmeval_rebuild_no_swe/env_candidates/
##          chelsa_candidates_no_swe_japan.tif
##
##   2. Pollinator SDM predictors use ONE metric at a time:
##        "sum_cloglog"
##        "mean_cloglog"
##        "binary_richness"
##
##   3. Main analysis uses sum_cloglog by default.
##      Sensitivity analysis can also run mean_cloglog and binary_richness.
##
##   4. Multicollinearity control:
##      - Do not mix sum / mean / binary in one model.
##      - Drop hawkmoth by default.
##      - Reduce guild predictors by pairwise correlation.
##      - Screen pollinator predictors highly correlated with env PCs.
##      - Report VIF for final candidate sets.
##
## Output:
##   chelsa_pollinator_enmeval_rebuild_no_swe/
##     cirsium_candidateEnvPCA_pollinator_metric_reduced/
############################################################

## =========================================================
## 0. Packages
## =========================================================

pkgs <- c(
  "terra", "dplyr", "readr", "stringr", "tibble", "tidyr",
  "ggplot2", "broom", "pROC", "usdm", "patchwork", "forcats",
  "randomForest"
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
  library(broom)
  library(pROC)
  library(usdm)
  library(patchwork)
  library(forcats)
  library(randomForest)
})

set.seed(42)

## =========================================================
## 1. Settings
## =========================================================

ROOT_DIR <- "chelsa_pollinator_enmeval_rebuild_no_swe"

ENV_CANDIDATE_TIF <- file.path(
  ROOT_DIR,
  "env_candidates",
  "chelsa_candidates_no_swe_japan.tif"
)

ENV_CANDIDATE_FILELIST <- file.path(
  ROOT_DIR,
  "env_candidates",
  "chelsa_candidates_no_swe_filelist.csv"
)

CIRSIUM_OCC_WITH_GUILD_FILE <- file.path(
  ROOT_DIR,
  "cirsium_occurrences_with_guild_sdm_speciesM_no_swe.csv"
)

OUT_ROOT <- file.path(
  ROOT_DIR,
  "cirsium_candidateEnvPCA_pollinator_metric_reduced"
)

FIG_ROOT <- file.path(OUT_ROOT, "figures")

dir.create(OUT_ROOT, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_ROOT, showWarnings = FALSE, recursive = TRUE)

## Environmental PCs used in GLM
N_ENV_PC_USE <- 3
N_ENV_PC_SAVE_LOADINGS <- 6

## Pollinator metric:
## Main recommendation = "sum_cloglog".
## Other options = "binary_richness", "mean_cloglog".
MAIN_POLLINATOR_METRIC <- "sum_cloglog"

## If TRUE, run all three metrics separately for sensitivity analysis.
RUN_SENSITIVITY_METRICS <- TRUE
METRICS_TO_RUN <- if (RUN_SENSITIVITY_METRICS) {
  c("sum_cloglog", "binary_richness", "mean_cloglog")
} else {
  MAIN_POLLINATOR_METRIC
}

## Multicollinearity thresholds
DROP_HAWKMOTH <- TRUE
POLL_GUILD_COR_THRESHOLD <- 0.70
POLL_ENVPC_COR_THRESHOLD <- 0.70

## For final model, if VIF is high, this script reports it.
## It does not automatically delete env PCs because PCs are orthogonal by construction.
VIF_THRESHOLD_REPORT <- 5

## Pollinator priority for keeping variables when guild predictors are correlated.
## Keep biologically central variables first.
POLL_PRIORITY_REGEX <- c(
  "bee_hanging_effective",
  "bee_non_hanging",
  "butterfly",
  "hoverfly"
)

## =========================================================
## 2. Helper functions
## =========================================================

stop_if_missing <- function(path) {
  if (!file.exists(path)) {
    stop("File not found: ", path)
  }
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

  if (!is.finite(s) || s == 0) {
    return(rep(0, length(x)))
  }

  (x - m) / s
}

calc_auc <- function(y, pred) {
  tryCatch(
    as.numeric(pROC::auc(y, pred, quiet = TRUE)),
    error = function(e) NA_real_
  )
}

fit_glm_safe <- function(formula, data, model_name) {
  tryCatch({
    m <- stats::glm(formula, data = data, family = stats::binomial)
    pred <- stats::predict(m, type = "response")

    list(
      model = m,
      aic = tibble::tibble(
        model = model_name,
        AIC = stats::AIC(m),
        BIC = stats::BIC(m),
        n = stats::nobs(m),
        logLik = as.numeric(stats::logLik(m)),
        df = attr(stats::logLik(m), "df"),
        AUC = calc_auc(data$y_nodding, pred)
      ),
      coef = broom::tidy(m) |>
        dplyr::mutate(model = model_name, .before = 1)
    )
  }, error = function(e) {
    list(
      model = NULL,
      aic = tibble::tibble(
        model = model_name,
        AIC = NA_real_,
        BIC = NA_real_,
        n = NA_integer_,
        logLik = NA_real_,
        df = NA_real_,
        AUC = NA_real_,
        error = conditionMessage(e)
      ),
      coef = tibble::tibble(
        model = model_name,
        term = NA_character_,
        estimate = NA_real_,
        std.error = NA_real_,
        statistic = NA_real_,
        p.value = NA_real_,
        error = conditionMessage(e)
      )
    )
  })
}

select_by_correlation <- function(data, vars, priority, threshold = 0.70) {
  vars <- vars[vars %in% names(data)]
  priority <- unique(priority[priority %in% vars])
  priority <- c(priority, setdiff(vars, priority))

  selected <- character(0)
  dropped <- tibble::tibble(
    variable = character(),
    dropped_because = character(),
    correlation = numeric()
  )

  for (v in priority) {
    if (length(selected) == 0) {
      selected <- c(selected, v)
      next
    }

    cors <- sapply(selected, function(s) {
      suppressWarnings(stats::cor(data[[v]], data[[s]], use = "pairwise.complete.obs"))
    })

    if (all(!is.finite(cors))) {
      selected <- c(selected, v)
      next
    }

    max_abs <- max(abs(cors), na.rm = TRUE)

    if (is.finite(max_abs) && max_abs >= threshold) {
      s_max <- names(which.max(abs(cors)))
      dropped <- dplyr::bind_rows(
        dropped,
        tibble::tibble(
          variable = v,
          dropped_because = s_max,
          correlation = unname(cors[s_max])
        )
      )
    } else {
      selected <- c(selected, v)
    }
  }

  list(selected = selected, dropped = dropped)
}

vif_check <- function(data, vars, label) {
  vars <- vars[vars %in% names(data)]

  if (length(vars) < 2) {
    return(tibble::tibble(model_set = label, variable = vars, VIF = NA_real_))
  }

  x <- data |>
    dplyr::select(dplyr::all_of(vars)) |>
    tidyr::drop_na()

  x <- x[, vapply(x, function(z) stats::sd(z, na.rm = TRUE) > 0, logical(1)), drop = FALSE]

  if (ncol(x) < 2) {
    return(tibble::tibble(model_set = label, variable = names(x), VIF = NA_real_))
  }

  v <- tryCatch(usdm::vif(x), error = function(e) NULL)

  if (is.null(v)) {
    return(tibble::tibble(model_set = label, variable = names(x), VIF = NA_real_))
  }

  as.data.frame(v) |>
    tibble::as_tibble() |>
    dplyr::rename(variable = Variables, VIF = VIF) |>
    dplyr::mutate(model_set = label, .before = 1)
}

metric_to_regex <- function(metric) {
  if (metric == "sum_cloglog") return("^z_guild_.*sum_cloglog")
  if (metric == "mean_cloglog") return("^z_guild_.*mean_cloglog")
  if (metric == "binary_richness") return("^z_guild_.*binary_richness")
  stop("Unknown metric: ", metric)
}

## =========================================================
## 3. Load pre-VIF candidate environmental stack
## =========================================================

stop_if_missing(CIRSIUM_OCC_WITH_GUILD_FILE)

if (file.exists(ENV_CANDIDATE_TIF)) {
  message("Loading pre-VIF candidate environmental stack: ", ENV_CANDIDATE_TIF)
  env_candidate <- terra::rast(ENV_CANDIDATE_TIF)
} else if (file.exists(ENV_CANDIDATE_FILELIST)) {
  message("Candidate stack tif not found; rebuilding from file list.")
  fl <- readr::read_csv(ENV_CANDIDATE_FILELIST, show_col_types = FALSE)

  if (!all(c("variable", "cropped_file") %in% names(fl))) {
    stop("Candidate file list must contain variable and cropped_file columns.")
  }

  if (!all(file.exists(fl$cropped_file))) {
    missing_files <- fl$cropped_file[!file.exists(fl$cropped_file)]
    stop("Missing candidate cropped files:\n", paste(missing_files, collapse = "\n"))
  }

  env_candidate <- terra::rast(fl$cropped_file)
  names(env_candidate) <- fl$variable
} else {
  stop("No pre-VIF candidate environmental stack found.")
}

if ("swe" %in% names(env_candidate)) {
  stop("swe is still in candidate environment stack. Use no-SWE rebuild output.")
}

message("Candidate environmental layers used for PCA:")
print(names(env_candidate))

## =========================================================
## 4. Load Cirsium occurrence + pollinator guild SDM data
## =========================================================

cir_occ <- readr::read_csv(CIRSIUM_OCC_WITH_GUILD_FILE, show_col_types = FALSE) |>
  dplyr::mutate(
    decimalLongitude = to_numeric_coord(decimalLongitude),
    decimalLatitude = to_numeric_coord(decimalLatitude)
  ) |>
  dplyr::filter(
    is.finite(decimalLongitude),
    is.finite(decimalLatitude),
    head_orientation_binary %in% c("upward", "nodding")
  ) |>
  dplyr::mutate(
    y_nodding = ifelse(head_orientation_binary == "nodding", 1, 0)
  )

message("Cirsium occurrence rows: ", nrow(cir_occ))
print(cir_occ |> dplyr::count(head_orientation_binary))

## =========================================================
## 5. Extract PRE-VIF candidate environmental variables
## =========================================================

xy <- cbind(cir_occ$decimalLongitude, cir_occ$decimalLatitude)
colnames(xy) <- c("lon", "lat")

env_vals <- safe_extract_df(env_candidate, xy)

## Use a unique prefix to avoid collision with old env_* columns.
names(env_vals) <- paste0("pcaenvCand_", names(env_vals))
env_vars <- names(env_vals)

cir_occ_env <- dplyr::bind_cols(cir_occ, env_vals, .name_repair = "unique")

readr::write_csv(
  cir_occ_env,
  file.path(OUT_ROOT, "01_cirsium_occurrences_with_candidateEnv_and_guildSDM.csv")
)

guild_cols <- names(cir_occ_env)[
  stringr::str_detect(
    names(cir_occ_env),
    "guild_.*(sum_cloglog|mean_cloglog|binary_richness)"
  )
]

message("Detected guild SDM columns:")
print(guild_cols)

numeric_mean_cols <- c(env_vars, guild_cols)
numeric_mean_cols <- numeric_mean_cols[numeric_mean_cols %in% names(cir_occ_env)]

species_dat <- cir_occ_env |>
  dplyr::group_by(species_final, japanese_name, head_orientation_binary) |>
  dplyr::summarise(
    y_nodding = dplyr::first(y_nodding),
    n_occ = dplyr::n(),
    dplyr::across(
      dplyr::all_of(numeric_mean_cols),
      ~ mean(.x, na.rm = TRUE),
      .names = "{.col}"
    ),
    mean_lon = mean(decimalLongitude, na.rm = TRUE),
    mean_lat = mean(decimalLatitude, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    dplyr::across(where(is.numeric), ~ ifelse(is.nan(.x), NA_real_, .x))
  )

readr::write_csv(
  species_dat,
  file.path(OUT_ROOT, "02_species_mean_candidateEnv_guildSDM.csv")
)

message("Species-level counts:")
print(species_dat |> dplyr::count(head_orientation_binary))

## =========================================================
## 6. PCA using PRE-VIF candidate variables
## =========================================================

env_vars <- env_vars[env_vars %in% names(species_dat)]

env_mat <- species_dat |>
  dplyr::select(dplyr::all_of(env_vars)) |>
  dplyr::select(where(is.numeric))

good_env_cols <- names(env_mat)[
  vapply(env_mat, function(x) {
    sx <- stats::sd(x, na.rm = TRUE)
    sum(is.finite(x)) >= 5 && is.finite(sx) && sx > 0
  }, logical(1))
]

env_mat <- env_mat[, good_env_cols, drop = FALSE]

if (ncol(env_mat) < 2) {
  stop("Too few pre-VIF candidate environmental variables for PCA.")
}

env_mat_imp <- env_mat
for (nm in names(env_mat_imp)) {
  x <- env_mat_imp[[nm]]
  x[!is.finite(x)] <- mean(x, na.rm = TRUE)
  env_mat_imp[[nm]] <- x
}

pca <- stats::prcomp(env_mat_imp, center = TRUE, scale. = TRUE)
saveRDS(pca, file.path(OUT_ROOT, "03_candidateEnv_PCA_model.rds"))

pca_scores <- as.data.frame(pca$x) |>
  tibble::as_tibble() |>
  dplyr::bind_cols(
    species_dat |>
      dplyr::select(
        species_final, japanese_name, head_orientation_binary,
        y_nodding, n_occ, mean_lon, mean_lat
      )
  ) |>
  dplyr::relocate(
    species_final, japanese_name, head_orientation_binary,
    y_nodding, n_occ, mean_lon, mean_lat
  )

pca_loadings <- as.data.frame(pca$rotation) |>
  tibble::rownames_to_column("variable") |>
  tibble::as_tibble()

pca_var <- tibble::tibble(
  PC = paste0("PC", seq_along(pca$sdev)),
  variance = pca$sdev^2,
  prop_variance = variance / sum(variance),
  cum_variance = cumsum(prop_variance)
)

top_loadings <- pca_loadings |>
  tidyr::pivot_longer(cols = dplyr::starts_with("PC"), names_to = "PC", values_to = "loading") |>
  dplyr::filter(PC %in% paste0("PC", 1:min(6, ncol(pca$x)))) |>
  dplyr::group_by(PC) |>
  dplyr::slice_max(order_by = abs(loading), n = 12, with_ties = FALSE) |>
  dplyr::arrange(PC, dplyr::desc(abs(loading))) |>
  dplyr::ungroup()

readr::write_csv(pca_scores, file.path(OUT_ROOT, "03_candidateEnv_PCA_scores.csv"))
readr::write_csv(pca_loadings, file.path(OUT_ROOT, "03_candidateEnv_PCA_loadings.csv"))
readr::write_csv(pca_var, file.path(OUT_ROOT, "03_candidateEnv_PCA_variance_explained.csv"))
readr::write_csv(top_loadings, file.path(OUT_ROOT, "03_candidateEnv_PCA_top_loadings.csv"))

p_pca <- ggplot2::ggplot(
  pca_scores,
  ggplot2::aes(PC1, PC2, color = head_orientation_binary)
) +
  ggplot2::geom_point(ggplot2::aes(size = n_occ), alpha = 0.85) +
  ggplot2::theme_bw(base_size = 13) +
  ggplot2::labs(
    x = paste0("PC1 (", round(pca_var$prop_variance[1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(pca_var$prop_variance[2] * 100, 1), "%)"),
    color = "Head orientation",
    size = "n occurrence",
    title = "PCA from pre-VIF candidate environmental variables"
  )

ggplot2::ggsave(
  file.path(FIG_ROOT, "candidateEnv_PCA_PC1_PC2.png"),
  p_pca, width = 6, height = 5, dpi = 300
)

## =========================================================
## 7. Build analysis table with PCs and z-scored guild variables
## =========================================================

analysis_dat <- species_dat |>
  dplyr::select(
    species_final, japanese_name, head_orientation_binary,
    y_nodding, n_occ, mean_lon, mean_lat,
    dplyr::all_of(guild_cols)
  ) |>
  dplyr::left_join(
    pca_scores |>
      dplyr::select(species_final, dplyr::starts_with("PC")),
    by = "species_final"
  )

## For M-restricted maps:
## NA in sum/richness means no predicted accessible guild species there.
sum_cols <- guild_cols[stringr::str_detect(guild_cols, "sum_cloglog")]
rich_cols <- guild_cols[stringr::str_detect(guild_cols, "binary_richness")]

analysis_dat <- analysis_dat |>
  dplyr::mutate(
    dplyr::across(
      dplyr::all_of(c(sum_cols, rich_cols)),
      ~ ifelse(is.na(.x), 0, .x)
    )
  )

pc_cols <- paste0("PC", 1:min(N_ENV_PC_USE, ncol(pca$x)))
pc_cols <- pc_cols[pc_cols %in% names(analysis_dat)]

for (nm in pc_cols) {
  analysis_dat[[paste0("z_", nm)]] <- safe_scale(analysis_dat[[nm]])
}

for (nm in guild_cols) {
  if (nm %in% names(analysis_dat)) {
    analysis_dat[[paste0("z_", nm)]] <- safe_scale(analysis_dat[[nm]])
  }
}

readr::write_csv(
  analysis_dat,
  file.path(OUT_ROOT, "04_species_analysis_data_candidateEnvPC_allGuildMetrics.csv")
)

## =========================================================
## 8. Function: run reduced pollinator analysis for one metric
## =========================================================

run_metric_analysis <- function(metric, analysis_dat, pc_cols, guild_cols) {
  message("\n=======================================")
  message("Running pollinator metric analysis: ", metric)
  message("=======================================")

  metric_safe <- safe_name(metric)
  metric_out <- file.path(OUT_ROOT, paste0("metric_", metric_safe))
  metric_fig <- file.path(metric_out, "figures")
  dir.create(metric_out, showWarnings = FALSE, recursive = TRUE)
  dir.create(metric_fig, showWarnings = FALSE, recursive = TRUE)

  z_pc_vars <- paste0("z_", pc_cols)
  z_pc_vars <- z_pc_vars[z_pc_vars %in% names(analysis_dat)]

  poll_regex <- metric_to_regex(metric)
  poll_primary <- names(analysis_dat)[stringr::str_detect(names(analysis_dat), poll_regex)]

  if (DROP_HAWKMOTH) {
    poll_primary <- poll_primary[!stringr::str_detect(poll_primary, "hawkmoth")]
  }

  if (length(poll_primary) == 0) {
    warning("No pollinator variables found for metric: ", metric)
    return(NULL)
  }

  poll_priority <- character(0)
  for (rx in POLL_PRIORITY_REGEX) {
    poll_priority <- c(poll_priority, poll_primary[stringr::str_detect(poll_primary, rx)])
  }
  poll_priority <- unique(c(poll_priority, poll_primary))

  ## Guild-guild correlation reduction
  corr_sel <- select_by_correlation(
    data = analysis_dat,
    vars = poll_primary,
    priority = poll_priority,
    threshold = POLL_GUILD_COR_THRESHOLD
  )

  poll_selected_guildcorr <- corr_sel$selected
  poll_dropped_guildcorr <- corr_sel$dropped

  readr::write_csv(
    tibble::tibble(metric = metric, variable = poll_primary),
    file.path(metric_out, "00_pollinator_primary_candidates.csv")
  )

  readr::write_csv(
    poll_dropped_guildcorr,
    file.path(metric_out, "01_pollinator_dropped_by_guild_correlation.csv")
  )

  ## Correlation with env PCs
  env_corr_records <- list()

  for (pv in poll_selected_guildcorr) {
    for (ev in z_pc_vars) {
      if (pv %in% names(analysis_dat) && ev %in% names(analysis_dat)) {
        env_corr_records[[paste(pv, ev)]] <- tibble::tibble(
          pollinator_var = pv,
          env_PC = ev,
          correlation = suppressWarnings(stats::cor(
            analysis_dat[[pv]], analysis_dat[[ev]],
            use = "pairwise.complete.obs"
          ))
        )
      }
    }
  }

  env_corr_tbl <- dplyr::bind_rows(env_corr_records) |>
    dplyr::arrange(dplyr::desc(abs(correlation)))

  readr::write_csv(
    env_corr_tbl,
    file.path(metric_out, "02_pollinator_vs_candidateEnvPC_correlations.csv")
  )

  high_env_corr_poll <- env_corr_tbl |>
    dplyr::filter(abs(correlation) > POLL_ENVPC_COR_THRESHOLD) |>
    dplyr::pull(pollinator_var) |>
    unique()

  poll_selected_independent <- setdiff(poll_selected_guildcorr, high_env_corr_poll)

  if (length(poll_selected_independent) == 0 && length(poll_selected_guildcorr) > 0) {
    poll_selected_independent <- poll_selected_guildcorr[1]
    warning(
      "All pollinator variables for metric ", metric,
      " are highly correlated with env PCs. Keeping first-priority variable: ",
      poll_selected_independent
    )
  }

  ## Wide selected predictor summary without duplicate rows
  selected_summary <- tibble::tibble(
    variable = unique(c(z_pc_vars, poll_primary, poll_selected_guildcorr, poll_selected_independent))
  ) |>
    dplyr::mutate(
      metric = metric,
      is_env_PC = variable %in% z_pc_vars,
      is_pollinator_primary_candidate = variable %in% poll_primary,
      retained_after_guild_correlation = variable %in% poll_selected_guildcorr,
      high_corr_with_envPC = variable %in% high_env_corr_poll,
      retained_for_env_plus_pollinator = variable %in% poll_selected_independent
    )

  readr::write_csv(
    selected_summary,
    file.path(metric_out, "00_selected_predictors_summary_wide.csv")
  )

  ## VIF and correlations
  candidate_sets <- list(
    env_PC = z_pc_vars,
    pollinator_primary_metric_only = poll_primary,
    pollinator_after_guild_correlation = poll_selected_guildcorr,
    envPC_plus_pollinator_independent = c(z_pc_vars, poll_selected_independent)
  )

  vif_tbl <- dplyr::bind_rows(lapply(names(candidate_sets), function(nm) {
    vif_check(analysis_dat, candidate_sets[[nm]], nm)
  }))

  readr::write_csv(
    vif_tbl,
    file.path(metric_out, "03_vif_candidate_sets.csv")
  )

  corr_vars <- unique(c(z_pc_vars, poll_selected_guildcorr))
  corr_vars <- corr_vars[corr_vars %in% names(analysis_dat)]

  corr_mat <- analysis_dat |>
    dplyr::select(dplyr::all_of(corr_vars)) |>
    as.data.frame()

  corr_mat <- corr_mat[
    ,
    vapply(corr_mat, function(x) stats::sd(x, na.rm = TRUE) > 0, logical(1)),
    drop = FALSE
  ]

  corr <- stats::cor(corr_mat, use = "pairwise.complete.obs")

  readr::write_csv(
    as.data.frame(corr) |> tibble::rownames_to_column("variable"),
    file.path(metric_out, "03_correlation_matrix_reduced_predictors.csv")
  )

  corr_pairs <- as.data.frame(as.table(corr)) |>
    tibble::as_tibble() |>
    dplyr::rename(var1 = Var1, var2 = Var2, correlation = Freq) |>
    dplyr::filter(as.character(var1) < as.character(var2)) |>
    dplyr::arrange(dplyr::desc(abs(correlation)))

  readr::write_csv(
    corr_pairs,
    file.path(metric_out, "03_correlation_pairs_reduced_predictors.csv")
  )

  ## GLM model set
  fml <- function(rhs) stats::as.formula(paste("y_nodding ~", rhs))

  model_specs <- list(null = "1")

  if ("z_PC1" %in% z_pc_vars) {
    model_specs$env_PC1 <- "z_PC1"
  }

  if (all(c("z_PC1", "z_PC2") %in% z_pc_vars)) {
    model_specs$env_PC1_PC2 <- "z_PC1 + z_PC2"
  }

  if (all(c("z_PC1", "z_PC2", "z_PC3") %in% z_pc_vars)) {
    model_specs$env_PC1_PC2_PC3 <- "z_PC1 + z_PC2 + z_PC3"
  }

  ## Pollinator variables individually
  for (pv in poll_selected_guildcorr) {
    clean_nm <- stringr::str_remove(pv, "^z_guild_")
    model_specs[[paste0("poll_", clean_nm)]] <- pv
  }

  ## Pollinator reduced combined model
  if (length(poll_selected_guildcorr) >= 2) {
    model_specs$pollinator_reduced <- paste(poll_selected_guildcorr, collapse = " + ")
  } else if (length(poll_selected_guildcorr) == 1) {
    model_specs$pollinator_reduced <- poll_selected_guildcorr
  }

  ## Env + independent pollinator
  if (length(poll_selected_independent) >= 1) {
    if ("z_PC1" %in% z_pc_vars) {
      model_specs$env_PC1_plus_pollinator_independent <- paste(
        c("z_PC1", poll_selected_independent),
        collapse = " + "
      )
    }

    if (all(c("z_PC1", "z_PC2") %in% z_pc_vars)) {
      model_specs$env_PC1_PC2_plus_pollinator_independent <- paste(
        c("z_PC1", "z_PC2", poll_selected_independent),
        collapse = " + "
      )
    }
  }

  ## Specific bee contrast
  bee_hang <- poll_selected_guildcorr[stringr::str_detect(poll_selected_guildcorr, "bee_hanging_effective")]
  bee_non <- poll_selected_guildcorr[stringr::str_detect(poll_selected_guildcorr, "bee_non_hanging")]

  if (length(bee_hang) >= 1 && length(bee_non) >= 1) {
    model_specs$bee_function_reduced <- paste(c(bee_hang[1], bee_non[1]), collapse = " + ")

    if ("z_PC1" %in% z_pc_vars) {
      model_specs$env_PC1_plus_bee_function_reduced <- paste(
        c("z_PC1", bee_hang[1], bee_non[1]),
        collapse = " + "
      )
    }
  }

  fits <- list()

  for (nm in names(model_specs)) {
    rhs <- model_specs[[nm]]

    formula_i <- if (rhs == "1") {
      stats::as.formula("y_nodding ~ 1")
    } else {
      fml(rhs)
    }

    fits[[nm]] <- fit_glm_safe(formula_i, analysis_dat, nm)
  }

  aic_tbl <- dplyr::bind_rows(lapply(fits, `[[`, "aic")) |>
    dplyr::arrange(AIC) |>
    dplyr::mutate(
      metric = metric,
      delta_AIC = AIC - min(AIC, na.rm = TRUE),
      akaike_weight = exp(-0.5 * delta_AIC) / sum(exp(-0.5 * delta_AIC), na.rm = TRUE),
      .before = 1
    )

  coef_tbl <- dplyr::bind_rows(lapply(fits, `[[`, "coef")) |>
    dplyr::mutate(metric = metric, .before = 1)

  readr::write_csv(
    aic_tbl,
    file.path(metric_out, "04_glm_metric_reduced_aic.csv")
  )

  readr::write_csv(
    coef_tbl,
    file.path(metric_out, "04_glm_metric_reduced_coefficients.csv")
  )

  ## Response plots
  plot_response <- function(dat, xvar, xlab = xvar) {
    if (!xvar %in% names(dat)) return(NULL)

    dd <- dat |>
      dplyr::filter(is.finite(.data[[xvar]]), is.finite(y_nodding))

    if (nrow(dd) < 5) return(NULL)

    ggplot2::ggplot(dd, ggplot2::aes(x = .data[[xvar]], y = y_nodding)) +
      ggplot2::geom_jitter(height = 0.05, width = 0, alpha = 0.65) +
      ggplot2::geom_smooth(
        method = "glm",
        method.args = list(family = stats::binomial),
        se = TRUE
      ) +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::labs(x = xlab, y = "Probability of nodding", title = xlab)
  }

  plot_vars <- c("z_PC1", "z_PC2", poll_selected_guildcorr)
  plot_vars <- plot_vars[plot_vars %in% names(analysis_dat)]

  plots <- lapply(plot_vars, function(v) plot_response(analysis_dat, v, v))
  plots <- plots[!vapply(plots, is.null, logical(1))]

  if (length(plots) > 0) {
    p_all <- patchwork::wrap_plots(plots, ncol = 2)
    ggplot2::ggsave(
      file.path(metric_fig, paste0("response_plots_", metric_safe, ".png")),
      p_all, width = 9, height = 8, dpi = 300
    )
  }

  ## Top model coefficient plot
  top_model_name <- aic_tbl$model[1]

  top_coef <- coef_tbl |>
    dplyr::filter(model == top_model_name, term != "(Intercept)") |>
    dplyr::mutate(
      conf.low = estimate - 1.96 * std.error,
      conf.high = estimate + 1.96 * std.error,
      term = forcats::fct_reorder(term, estimate)
    )

  if (nrow(top_coef) > 0) {
    p_coef <- ggplot2::ggplot(top_coef, ggplot2::aes(x = estimate, y = term)) +
      ggplot2::geom_vline(xintercept = 0, linetype = 2) +
      ggplot2::geom_errorbarh(ggplot2::aes(xmin = conf.low, xmax = conf.high), height = 0.15) +
      ggplot2::geom_point(size = 2) +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::labs(
        x = "Logit coefficient",
        y = NULL,
        title = paste0("Top GLM: ", top_model_name, " / ", metric)
      )

    ggplot2::ggsave(
      file.path(metric_fig, paste0("top_model_coefficients_", metric_safe, ".png")),
      p_coef, width = 7, height = 4, dpi = 300
    )
  }

  ## RF exploratory
  rf_vars <- unique(c(z_pc_vars, poll_selected_guildcorr))
  rf_vars <- rf_vars[rf_vars %in% names(analysis_dat)]

  rf_dat <- analysis_dat |>
    dplyr::select(head_orientation_binary, dplyr::all_of(rf_vars)) |>
    tidyr::drop_na() |>
    dplyr::mutate(head_orientation_binary = factor(head_orientation_binary))

  rf_importance <- NULL

  if (nrow(rf_dat) >= 10 && length(unique(rf_dat$head_orientation_binary)) == 2) {
    rf <- randomForest::randomForest(
      head_orientation_binary ~ .,
      data = rf_dat,
      importance = TRUE,
      ntree = 1000
    )

    saveRDS(rf, file.path(metric_out, paste0("05_random_forest_", metric_safe, ".rds")))

    rf_importance <- as.data.frame(randomForest::importance(rf)) |>
      tibble::rownames_to_column("variable") |>
      tibble::as_tibble() |>
      dplyr::mutate(metric = metric, .before = 1)

    readr::write_csv(
      rf_importance,
      file.path(metric_out, paste0("05_random_forest_importance_", metric_safe, ".csv"))
    )
  }

  message("\nAIC table for metric: ", metric)
  print(aic_tbl)

  message("\nTop coefficients for metric: ", metric)
  print(
    coef_tbl |>
      dplyr::filter(term != "(Intercept)") |>
      dplyr::arrange(p.value),
    n = 100
  )

  list(
    metric = metric,
    selected_summary = selected_summary,
    dropped_guildcorr = poll_dropped_guildcorr,
    env_corr = env_corr_tbl,
    vif = vif_tbl,
    aic = aic_tbl,
    coef = coef_tbl,
    rf_importance = rf_importance
  )
}

## =========================================================
## 9. Run selected metrics
## =========================================================

metric_results <- list()

for (metric in METRICS_TO_RUN) {
  metric_results[[metric]] <- run_metric_analysis(
    metric = metric,
    analysis_dat = analysis_dat,
    pc_cols = pc_cols,
    guild_cols = guild_cols
  )
}

## Combined summaries across metrics
combined_aic <- dplyr::bind_rows(lapply(metric_results, `[[`, "aic"))
combined_coef <- dplyr::bind_rows(lapply(metric_results, `[[`, "coef"))
combined_vif <- dplyr::bind_rows(lapply(metric_results, `[[`, "vif"))
combined_selected <- dplyr::bind_rows(lapply(metric_results, `[[`, "selected_summary"))
combined_env_corr <- dplyr::bind_rows(lapply(metric_results, `[[`, "env_corr"))

readr::write_csv(combined_aic, file.path(OUT_ROOT, "ALL_metrics_glm_aic.csv"))
readr::write_csv(combined_coef, file.path(OUT_ROOT, "ALL_metrics_glm_coefficients.csv"))
readr::write_csv(combined_vif, file.path(OUT_ROOT, "ALL_metrics_vif.csv"))
readr::write_csv(combined_selected, file.path(OUT_ROOT, "ALL_metrics_selected_predictors_summary.csv"))
readr::write_csv(combined_env_corr, file.path(OUT_ROOT, "ALL_metrics_pollinator_vs_envPC_correlations.csv"))

saveRDS(metric_results, file.path(OUT_ROOT, "ALL_metric_results.rds"))

message("\n🎉 DONE — candidate-env PCA + one-metric-at-a-time reduced pollinator analysis")
message("Output root: ", OUT_ROOT)
message("Main outputs:")
message("  03_candidateEnv_PCA_top_loadings.csv")
message("  04_species_analysis_data_candidateEnvPC_allGuildMetrics.csv")
message("  ALL_metrics_glm_aic.csv")
message("  ALL_metrics_glm_coefficients.csv")
message("  ALL_metrics_vif.csv")
message("  ALL_metrics_selected_predictors_summary.csv")
message("Metric-specific outputs:")
message("  metric_sum_cloglog/")
message("  metric_binary_richness/")
message("  metric_mean_cloglog/")
