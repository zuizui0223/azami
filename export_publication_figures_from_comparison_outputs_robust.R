############################################################
## Publication figures from pollinator grouping comparison outputs
## Robust version:
##   - Does NOT require scheme_<best_scheme>/best_model_coefficients.csv
##   - Uses top-level all_scheme_* CSVs first
##   - Re-fits the best model from scheme analysis_grid_data.csv if predictions
##     or coefficients are missing
##   - Uses top-level 01_grid_preVIF_env_PCA_loadings.csv for environmental PCA
##
## Run:
##   setwd("C:/Users/zuizui/cirsium_inat")
##   source("export_publication_figures_from_comparison_outputs_robust.R")
############################################################

pkgs <- c(
  "readr", "dplyr", "stringr", "tidyr", "ggplot2",
  "forcats", "tibble", "patchwork", "glmmTMB", "pROC"
)

for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p, dependencies = TRUE)
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(ggplot2)
  library(forcats)
  library(tibble)
  library(patchwork)
  library(glmmTMB)
  library(pROC)
})

ROOT_DIR <- "chelsa_pollinator_enmeval_rebuild_no_swe"

MODEL_COMP_DIR <- file.path(
  ROOT_DIR,
  "cirsium_pollinator_grouping_sensitivity_model_comparison"
)

PUB_FIG_DIR <- file.path(
  MODEL_COMP_DIR,
  "publication_figures_best_model"
)

dir.create(PUB_FIG_DIR, showWarnings = FALSE, recursive = TRUE)

BEST_METHOD <- "beta_binomial_glmmTMB"
N_TOP_MODELS <- 20

## Clean manuscript panels: no plot titles/subtitles.
SHOW_PLOT_TITLES <- FALSE

## Optionally force a model. Leave NA to auto-detect global best.
FORCE_BEST_SCHEME <- NA_character_
FORCE_BEST_MODEL <- NA_character_

safe_read_csv <- function(path, required = FALSE) {
  if (!file.exists(path)) {
    msg <- paste0("Missing file: ", path)
    if (required) stop(msg) else {
      warning(msg)
      return(NULL)
    }
  }
  readr::read_csv(path, show_col_types = FALSE)
}

safe_name <- function(x) {
  x |>
    as.character() |>
    stringr::str_replace_all("[^A-Za-z0-9_\\-]+", "_") |>
    stringr::str_replace_all("_+", "_") |>
    stringr::str_replace_all("^_|_$", "")
}

theme_paper <- function(base_size = 11) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      axis.line = ggplot2::element_line(linewidth = 0.35),
      axis.ticks = ggplot2::element_line(linewidth = 0.35),
      panel.grid.major.x = ggplot2::element_line(colour = "grey88", linewidth = 0.25),
      panel.grid.major.y = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(face = "bold"),
      strip.background = ggplot2::element_rect(fill = "grey92", colour = "grey55", linewidth = 0.35),
      strip.text = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_blank(),
      plot.subtitle = ggplot2::element_blank()
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

save_plot <- function(plot, filename, width = 8, height = 6) {
  png_file <- file.path(PUB_FIG_DIR, filename)
  pdf_file <- file.path(PUB_FIG_DIR, stringr::str_replace(filename, "\\.png$", ".pdf"))

  ggplot2::ggsave(png_file, plot, width = width, height = height, dpi = 400, bg = "white")
  ggplot2::ggsave(pdf_file, plot, width = width, height = height, bg = "white")

  invisible(c(png_file, pdf_file))
}

pretty_scheme <- function(x) {
  dplyr::case_when(
    x == "order3_pca" ~ "Order PCA\n(Hymenoptera/Lepidoptera/Diptera)",
    x == "lepidoptera3_direct" ~ "Direct 3 groups\n(bee/Lepidoptera/fly)",
    x == "broad4_direct" ~ "Direct 4 groups\n(bee/butterfly/fly/hawkmoth)",
    x == "order3_direct" ~ "Direct orders\n(Hymenoptera/Lepidoptera/Diptera)",
    x == "bombus_vs_nonbombus_all" ~ "Bombus vs\nnon-Bombus all",
    x == "bombus_vs_nonbombus_bees" ~ "Bombus vs\nnon-Bombus bees",
    x == "total_pollinator" ~ "Total pollinator",
    TRUE ~ x
  )
}

pretty_model <- function(x) {
  dplyr::case_when(
    x == "null" ~ "Null",
    x == "lat_only" ~ "Latitude only",
    x == "lonlat" ~ "Longitude + latitude",
    x == "space_poly" ~ "Space polynomial",
    x == "env" ~ "Environment",
    x == "poll" ~ "Pollinator",
    x == "env_plus_poll" ~ "Environment + pollinator",
    x == "env_plus_poll_plus_lat" ~ "Environment + pollinator + latitude",
    x == "env_plus_space_poly" ~ "Environment + space",
    x == "poll_plus_space_poly" ~ "Pollinator + space",
    x == "env_plus_poll_plus_space_poly" ~ "Environment + pollinator + space",
    TRUE ~ x
  )
}

pretty_term <- function(x) {
  dplyr::case_when(
    x == "z_env_PC1" ~ "Env PC1",
    x == "z_env_PC2" ~ "Env PC2",
    x == "z_env_PC3" ~ "Env PC3",
    x == "z_order_PC1" ~ "Order PC1",
    x == "z_order_PC2" ~ "Order PC2",
    x == "z_order_PC3" ~ "Order PC3",
    x == "z_poll_bombus_sum" ~ "Bombus suitability",
    x == "z_poll_non_bombus_sum" ~ "Non-Bombus suitability",
    x == "z_poll_non_bombus_bees_sum" ~ "Non-Bombus bee suitability",
    x == "z_poll_total_pollinator_sum" ~ "Total pollinator suitability",
    x == "z_poll_bee_sum" ~ "Bee suitability",
    x == "z_poll_lepidoptera_sum" ~ "Lepidoptera suitability",
    x == "z_poll_butterfly_sum" ~ "Butterfly suitability",
    x == "z_poll_hawkmoth_sum" ~ "Hawkmoth suitability",
    x == "z_poll_fly_sum" ~ "Fly suitability",
    x == "z_lat" ~ "Latitude",
    x == "z_lon" ~ "Longitude",
    x == "I(z_lat^2)" ~ "Latitude²",
    x == "I(z_lon^2)" ~ "Longitude²",
    x == "z_lon:z_lat" ~ "Longitude × Latitude",
    TRUE ~ x
  )
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

detect_pred_col <- function(dat) {
  cand <- c("pred_prob_nodding_best", "predicted_probability", "pred_prob_nodding", "pred")
  hit <- cand[cand %in% names(dat)]
  if (length(hit) == 0) NA_character_ else hit[1]
}

rhs_from_model <- function(model_name, dat) {
  env_vars <- paste0("z_env_PC", 1:3)
  env_vars <- env_vars[env_vars %in% names(dat)]

  if (any(grepl("^z_order_PC", names(dat)))) {
    poll_vars <- grep("^z_order_PC", names(dat), value = TRUE)
  } else {
    poll_vars <- grep("^z_poll_", names(dat), value = TRUE)
  }

  ## Keep only actual first few pollinator PCs when PCA-based.
  poll_vars <- poll_vars[!grepl("\\.\\.\\.", poll_vars)]

  space_poly <- c("z_lon", "z_lat", "I(z_lon^2)", "I(z_lat^2)", "z_lon:z_lat")

  vars <- switch(
    model_name,
    null = character(0),
    lat_only = "z_lat",
    lonlat = c("z_lon", "z_lat"),
    space_poly = space_poly,
    env = env_vars,
    poll = poll_vars,
    env_plus_poll = c(env_vars, poll_vars),
    env_plus_poll_plus_lat = c(env_vars, poll_vars, "z_lat"),
    env_plus_space_poly = c(env_vars, space_poly),
    poll_plus_space_poly = c(poll_vars, space_poly),
    env_plus_poll_plus_space_poly = c(env_vars, poll_vars, space_poly),
    c(env_vars, poll_vars, space_poly)
  )

  vars <- vars[!is.na(vars) & vars != ""]
  if (length(vars) == 0) "1" else paste(vars, collapse = " + ")
}

make_formula <- function(rhs) {
  stats::as.formula(paste("cbind(n_nodding_species, n_upward_species) ~", rhs))
}

coef_from_model <- function(m, scheme, model_name) {
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
    dplyr::mutate(
      scheme = scheme,
      model = model_name,
      method = "beta_binomial_glmmTMB",
      .before = 1
    )
}

## =========================================================
## 1. Detect best model from top-level comparison outputs
## =========================================================

bb_all <- safe_read_csv(
  file.path(MODEL_COMP_DIR, "all_scheme_beta_binomial_model_comparison.csv"),
  required = TRUE
)

if (!"global_delta_AIC" %in% names(bb_all)) {
  bb_all <- bb_all |>
    dplyr::mutate(global_delta_AIC = AIC - min(AIC, na.rm = TRUE))
}

if (!"global_akaike_weight" %in% names(bb_all)) {
  bb_all <- bb_all |>
    dplyr::mutate(
      global_akaike_weight = exp(-0.5 * global_delta_AIC) /
        sum(exp(-0.5 * global_delta_AIC), na.rm = TRUE)
    )
}

bb_ok <- bb_all |>
  dplyr::filter(method == BEST_METHOD, is.finite(AIC)) |>
  dplyr::arrange(AIC)

if (nrow(bb_ok) == 0) stop("No valid beta-binomial models found.")

if (!is.na(FORCE_BEST_SCHEME) && !is.na(FORCE_BEST_MODEL)) {
  best_scheme <- FORCE_BEST_SCHEME
  best_model <- FORCE_BEST_MODEL
  best_row <- bb_ok |>
    dplyr::filter(scheme == best_scheme, model == best_model) |>
    dplyr::slice(1)
} else {
  best_row <- bb_ok |> dplyr::slice(1)
  best_scheme <- best_row$scheme[1]
  best_model <- best_row$model[1]
}

if (nrow(best_row) == 0) stop("Best model not found in comparison table.")

SCHEME_DIR <- file.path(MODEL_COMP_DIR, paste0("scheme_", safe_name(best_scheme)))

message("Best scheme: ", best_scheme)
message("Best model : ", best_model)
message("Scheme dir : ", SCHEME_DIR)

write_csv(best_row, file.path(PUB_FIG_DIR, "00_global_best_model_row.csv"))

## =========================================================
## 2. Figure 1: global top model comparison
## =========================================================

top_models <- bb_ok |>
  dplyr::arrange(AIC) |>
  dplyr::slice_head(n = N_TOP_MODELS) |>
  dplyr::mutate(
    label = paste0(pretty_scheme(scheme), " / ", pretty_model(model)),
    label = forcats::fct_reorder(label, global_delta_AIC)
  )

p_global <- ggplot(top_models, aes(x = global_delta_AIC, y = label)) +
  geom_col(width = 0.72) +
  theme_paper(10) +
  paper_labs(
    x = expression(Global~Delta*AIC),
    y = NULL,
    title = "Model comparison across pollinator grouping schemes"
  )

save_plot(p_global, "Fig1_global_model_comparison_top20.png", width = 10.5, height = 7.5)

## =========================================================
## 3. Figure 2: best model per scheme
## =========================================================

best_each_scheme <- bb_ok |>
  dplyr::group_by(scheme) |>
  dplyr::slice_min(AIC, n = 1, with_ties = FALSE) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    best_scheme_delta_AIC = AIC - min(AIC, na.rm = TRUE),
    scheme_pretty = pretty_scheme(scheme),
    scheme_pretty = forcats::fct_reorder(scheme_pretty, best_scheme_delta_AIC)
  )

write_csv(best_each_scheme, file.path(PUB_FIG_DIR, "Fig2_best_model_each_scheme_data.csv"))

p_scheme <- ggplot(best_each_scheme, aes(x = best_scheme_delta_AIC, y = scheme_pretty)) +
  geom_col(width = 0.72) +
  theme_paper(11) +
  paper_labs(
    x = expression(Delta*AIC~"of best model in each scheme"),
    y = NULL,
    title = "Sensitivity to pollinator grouping scheme"
  )

save_plot(p_scheme, "Fig2_best_model_each_scheme.png", width = 8.2, height = 5.8)

## =========================================================
## 4. Best-model coefficients:
##    First use top-level all_scheme_beta_binomial_coefficients.csv.
##    If missing, refit from analysis_grid_data.csv.
## =========================================================

all_coef <- safe_read_csv(file.path(MODEL_COMP_DIR, "all_scheme_beta_binomial_coefficients.csv"))

coef <- NULL
if (!is.null(all_coef)) {
  coef <- all_coef |>
    dplyr::filter(
      scheme == best_scheme,
      model == best_model,
      method == BEST_METHOD
    )
  if (nrow(coef) == 0) coef <- NULL
}

analysis_file <- file.path(SCHEME_DIR, "analysis_grid_data.csv")
dat <- safe_read_csv(analysis_file, required = TRUE)

best_model_object <- NULL

if (is.null(coef)) {
  message("Best-model coefficients not found in top-level coefficient table. Re-fitting best model.")
  rhs <- rhs_from_model(best_model, dat)
  f <- make_formula(rhs)

  best_model_object <- glmmTMB::glmmTMB(
    f,
    data = dat,
    family = glmmTMB::betabinomial(link = "logit")
  )

  coef <- coef_from_model(best_model_object, best_scheme, best_model)
}

write_csv(coef, file.path(PUB_FIG_DIR, "Fig3_best_model_coefficients_data.csv"))

coef_plot <- coef |>
  dplyr::filter(term != "(Intercept)", is.finite(estimate), is.finite(std.error)) |>
  dplyr::mutate(
    conf.low = estimate - 1.96 * std.error,
    conf.high = estimate + 1.96 * std.error,
    term_pretty = pretty_term(term),
    term_pretty = forcats::fct_reorder(term_pretty, estimate)
  )

p_coef <- ggplot(coef_plot, aes(x = estimate, y = term_pretty)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.35) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.18, linewidth = 0.45) +
  geom_point(size = 2.4) +
  theme_paper(11) +
  paper_labs(
    x = "Coefficient on logit scale",
    y = NULL,
    title = "Predictors of nodding-flowered species proportion"
  )

save_plot(p_coef, "Fig3_best_model_coefficients.png", width = 8.8, height = 6.2)

## =========================================================
## 5. Figure 4: pollinator order PCA loadings, if available
## =========================================================

order_load_candidates <- c(
  file.path(SCHEME_DIR, "order_PCA_loadings.csv"),
  file.path(SCHEME_DIR, "pollinator_order_PCA_loadings.csv"),
  file.path(SCHEME_DIR, "03_order_PCA_loadings.csv")
)

order_load_file <- order_load_candidates[file.exists(order_load_candidates)][1]

if (!is.na(order_load_file)) {
  order_load <- readr::read_csv(order_load_file, show_col_types = FALSE)

  if (!"variable" %in% names(order_load)) {
    var_cand <- intersect(c("pollinator_variable", "order", "group", "var"), names(order_load))
    if (length(var_cand) > 0) names(order_load)[names(order_load) == var_cand[1]] <- "variable"
  }
  if (!"PC" %in% names(order_load)) {
    pc_cand <- intersect(c("component", "axis", "Dim"), names(order_load))
    if (length(pc_cand) > 0) names(order_load)[names(order_load) == pc_cand[1]] <- "PC"
  }

  if (any(stringr::str_detect(names(order_load), "^PC[0-9]+$"))) {
    order_long <- order_load |>
      dplyr::select(variable, dplyr::matches("^PC[0-9]+$")) |>
      tidyr::pivot_longer(cols = dplyr::matches("^PC[0-9]+$"), names_to = "PC", values_to = "loading")
  } else {
    order_long <- order_load |>
      dplyr::mutate(
        PC = as.character(PC),
        PC = ifelse(stringr::str_detect(PC, "^PC"), PC, paste0("PC", stringr::str_extract(PC, "[0-9]+")))
      ) |>
      dplyr::select(variable, PC, loading)
  }

  order_long <- order_long |>
    dplyr::mutate(
      variable_pretty = dplyr::case_when(
        stringr::str_detect(variable, regex("hymenoptera", ignore_case = TRUE)) ~ "Hymenoptera",
        stringr::str_detect(variable, regex("lepidoptera", ignore_case = TRUE)) ~ "Lepidoptera",
        stringr::str_detect(variable, regex("diptera", ignore_case = TRUE)) ~ "Diptera",
        TRUE ~ variable
      )
    ) |>
    dplyr::filter(PC %in% c("PC1", "PC2", "PC3"))

  write_csv(order_long, file.path(PUB_FIG_DIR, "Fig4_pollinator_order_PCA_loadings_data.csv"))

  p_order <- ggplot(order_long, aes(x = loading, y = variable_pretty)) +
    geom_vline(xintercept = 0, linetype = 2, linewidth = 0.35) +
    geom_col(width = 0.72) +
    facet_wrap(~ PC, scales = "free_x") +
    theme_paper(11) +
    paper_labs(
      x = "PCA loading",
      y = NULL,
      title = "Pollinator order composition axes"
    )

  save_plot(p_order, "Fig4_pollinator_order_PCA_loadings.png", width = 9.5, height = 5.2)
} else {
  message("No order PCA loading file found. Skipping Fig4.")
}

## =========================================================
## 6. Figure 5: environmental PCA top loadings
##    Use exactly the top-level grid environmental PCA loadings produced by
##    the comparison script. This avoids accidentally reading older outputs.
## =========================================================

env_load_file <- file.path(MODEL_COMP_DIR, "01_grid_preVIF_env_PCA_loadings.csv")

if (file.exists(env_load_file)) {
  env_wide <- readr::read_csv(env_load_file, show_col_types = FALSE)

  if (!"variable" %in% names(env_wide)) {
    var_cand <- intersect(c("var", "layer", "name"), names(env_wide))
    if (length(var_cand) > 0) names(env_wide)[names(env_wide) == var_cand[1]] <- "variable"
  }

  env_long <- env_wide |>
    dplyr::select(variable, dplyr::matches("^PC[0-9]+$")) |>
    tidyr::pivot_longer(
      cols = dplyr::matches("^PC[0-9]+$"),
      names_to = "PC",
      values_to = "loading"
    ) |>
    dplyr::mutate(
      variable = as.character(variable),
      variable = stringr::str_remove(variable, "^preVIF_"),
      variable = stringr::str_remove(variable, "^pcaenvCand_"),
      variable = stringr::str_remove(variable, "^pcaenv_"),
      variable = stringr::str_remove(variable, "^env_"),
      variable = stringr::str_replace_all(variable, "_", " "),
      loading = as.numeric(loading),
      abs_loading = abs(loading)
    ) |>
    dplyr::filter(is.finite(loading), PC %in% c("PC1", "PC2", "PC3"))

  env_top <- env_long |>
    dplyr::group_by(PC) |>
    dplyr::slice_max(abs_loading, n = 12, with_ties = FALSE) |>
    dplyr::arrange(PC, dplyr::desc(abs_loading)) |>
    dplyr::ungroup()

  write_csv(env_top, file.path(PUB_FIG_DIR, "Fig5_environmental_PCA_top_loadings_data.csv"))

  p_env <- env_top |>
    dplyr::mutate(variable = forcats::fct_reorder(variable, loading)) |>
    ggplot(aes(x = loading, y = variable)) +
    geom_vline(xintercept = 0, linetype = 2, linewidth = 0.35) +
    geom_col(width = 0.72) +
    facet_wrap(~ PC, scales = "free_y") +
    theme_paper(9.5) +
    paper_labs(
      x = "PCA loading",
      y = NULL,
      title = "Environmental PCA axes"
    )

  save_plot(p_env, "Fig5_environmental_PCA_top_loadings.png", width = 10.5, height = 7.8)
} else {
  warning("Environmental PCA loading file not found: ", env_load_file)
  writeLines(
    c(
      "Environmental PCA loading file not found.",
      paste0("Expected: ", env_load_file)
    ),
    con = file.path(PUB_FIG_DIR, "Fig5_environmental_PCA_loading_file_not_found.txt")
  )
}

## =========================================================
## 7. Prediction maps / residuals
##    Prefer existing predicted_probability_best_model.csv.
##    If missing, predict using re-fitted model.
## =========================================================

pred_file <- file.path(SCHEME_DIR, "predicted_probability_best_model.csv")

if (file.exists(pred_file)) {
  pred <- readr::read_csv(pred_file, show_col_types = FALSE)
} else {
  message("Prediction file not found. Creating predictions from analysis_grid_data.csv.")

  pred <- dat

  if (is.null(best_model_object)) {
    rhs <- rhs_from_model(best_model, dat)
    f <- make_formula(rhs)
    best_model_object <- glmmTMB::glmmTMB(
      f,
      data = dat,
      family = glmmTMB::betabinomial(link = "logit")
    )
  }

  pred$pred_prob_nodding_best <- as.numeric(predict(best_model_object, newdata = pred, type = "response"))
  pred$resid_prop_nodding_best <- pred$prop_nodding_species - pred$pred_prob_nodding_best

  write_csv(pred, file.path(PUB_FIG_DIR, "predicted_probability_best_model_recreated.csv"))
}

pred_col <- detect_pred_col(pred)
if (is.na(pred_col)) stop("Could not detect prediction column.")

if (!"resid_prop_nodding_best" %in% names(pred) && "prop_nodding_species" %in% names(pred)) {
  pred$resid_prop_nodding_best <- pred$prop_nodding_species - pred[[pred_col]]
}

p_pred <- ggplot(pred, aes(x = grid_lon, y = grid_lat, color = .data[[pred_col]], size = n_total_species)) +
  geom_point(alpha = 0.9) +
  coord_equal() +
  scale_color_viridis_c(limits = c(0, 1), option = "C") +
  theme_paper(11) +
  theme(panel.grid.major = element_line(colour = "grey90", linewidth = 0.25)) +
  paper_labs(
    x = "Longitude",
    y = "Latitude",
    color = "Predicted\nproportion",
    size = "Total\nspecies",
    title = "Predicted proportion of nodding-flowered species"
  )

save_plot(p_pred, "Fig6a_predicted_probability_map.png", width = 7.2, height = 6.4)

if ("resid_prop_nodding_best" %in% names(pred)) {
  p_resid <- ggplot(pred, aes(x = grid_lon, y = grid_lat, color = resid_prop_nodding_best, size = n_total_species)) +
    geom_point(alpha = 0.9) +
    coord_equal() +
    scale_color_gradient2(midpoint = 0) +
    theme_paper(11) +
    theme(panel.grid.major = element_line(colour = "grey90", linewidth = 0.25)) +
    paper_labs(
      x = "Longitude",
      y = "Latitude",
      color = "Observed -\npredicted",
      size = "Total\nspecies",
      title = "Residuals of nodding-species proportion"
    )

  save_plot(p_resid, "Fig6b_residual_map.png", width = 7.2, height = 6.4)
}

if ("prop_nodding_species" %in% names(pred)) {
  p_obs_pred <- ggplot(pred, aes(x = .data[[pred_col]], y = prop_nodding_species)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.35) +
    geom_point(aes(size = n_total_species), alpha = 0.75) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    theme_paper(11) +
    paper_labs(
      x = "Predicted proportion",
      y = "Observed proportion",
      size = "Total\nspecies",
      title = "Observed vs predicted proportion"
    )

  save_plot(p_obs_pred, "Fig6c_observed_vs_predicted.png", width = 6.2, height = 5.8)
}

if (exists("p_resid") && exists("p_obs_pred")) {
  p_combined <- (p_pred + p_resid) / p_obs_pred
  save_plot(p_combined, "Fig6_combined_prediction_residual_panel.png", width = 11, height = 10)
}

if ("standardized_excess_nodding" %in% names(pred)) {
  p_excess <- ggplot(pred, aes(x = grid_lon, y = grid_lat, color = standardized_excess_nodding, size = n_total_species)) +
    geom_point(alpha = 0.9) +
    coord_equal() +
    scale_color_gradient2(midpoint = 0) +
    theme_paper(11) +
    theme(panel.grid.major = element_line(colour = "grey90", linewidth = 0.25)) +
    paper_labs(
      x = "Longitude",
      y = "Latitude",
      color = "Standardized\nexcess",
      size = "Total\nspecies",
      title = "Hotspots and coldspots of nodding-flowered species"
    )

  save_plot(p_excess, "Fig7_standardized_excess_nodding_map.png", width = 7.2, height = 6.4)
}

## =========================================================
## 8. Output index and compact summary
## =========================================================

fig_index <- tibble::tibble(
  figure = c(
    "Fig1_global_model_comparison_top20",
    "Fig2_best_model_each_scheme",
    "Fig3_best_model_coefficients",
    "Fig4_pollinator_order_PCA_loadings",
    "Fig5_environmental_PCA_top_loadings",
    "Fig6a_predicted_probability_map",
    "Fig6b_residual_map",
    "Fig6c_observed_vs_predicted",
    "Fig6_combined_prediction_residual_panel",
    "Fig7_standardized_excess_nodding_map"
  ),
  png = file.path(PUB_FIG_DIR, paste0(figure, ".png")),
  pdf = file.path(PUB_FIG_DIR, paste0(figure, ".pdf")),
  exists_png = file.exists(png),
  exists_pdf = file.exists(pdf)
)

write_csv(fig_index, file.path(PUB_FIG_DIR, "00_publication_figure_index.csv"))

summary_tbl <- tibble::tibble(
  best_scheme = best_scheme,
  best_model = best_model,
  best_method = BEST_METHOD,
  best_AIC = best_row$AIC[1],
  global_delta_AIC = best_row$global_delta_AIC[1],
  global_akaike_weight = best_row$global_akaike_weight[1],
  n_grid = best_row$n_grid[1],
  AUC_presence = best_row$AUC_presence[1],
  logloss = best_row$logloss[1],
  brier = best_row$brier[1]
)

write_csv(summary_tbl, file.path(PUB_FIG_DIR, "00_best_model_summary.csv"))

message("\n🎉 DONE — publication-ready figures exported")
message("Best scheme: ", best_scheme)
message("Best model : ", best_model)
message("Output dir : ", PUB_FIG_DIR)
message("Figure index:")
print(fig_index, n = Inf, width = Inf)
