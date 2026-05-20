############################################################
## Cirsium head orientation grid analysis
## Pre-VIF environmental GRID PCA
## + merged Lepidoptera broad pollinator SSDMs
##
## Key design:
##   - Environmental PCA is rebuilt from PRE-VIF candidate
##     environmental variables at GRID level, not species level.
##   - Pollinator groups:
##       bee
##       lepidoptera = butterfly + hawkmoth
##       fly
##   - Response:
##       cbind(n_nodding_species, n_upward_species)
##     = regional/grid-level proportion of nodding Cirsium species.
##
## Inputs expected:
##   chelsa_pollinator_enmeval_rebuild_no_swe/
##     species_enmeval_run_summary.csv
##     env_candidates/chelsa_candidates_no_swe_japan.tif
##     cirsium_occurrences_with_guild_sdm_speciesM_no_swe.csv
##
## Output:
##   chelsa_pollinator_enmeval_rebuild_no_swe/
##     broad_group_bee_lepidoptera_fly_preVIF_gridEnvPCA_analysis/
############################################################

## =========================================================
## 0. Packages
## =========================================================

pkgs <- c(
  "terra", "dplyr", "readr", "stringr", "tibble", "tidyr",
  "ggplot2", "pROC", "forcats", "glmmTMB"
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
  library(forcats)
  library(glmmTMB)
})

set.seed(42)

## =========================================================
## 1. Paths and settings
## =========================================================

ROOT_DIR <- "chelsa_pollinator_enmeval_rebuild_no_swe"

SPECIES_SUMMARY_FILE <- file.path(ROOT_DIR, "species_enmeval_run_summary.csv")

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

CIRSIUM_OCC_FILE <- file.path(
  ROOT_DIR,
  "cirsium_occurrences_with_guild_sdm_speciesM_no_swe.csv"
)

OUT_DIR <- file.path(
  ROOT_DIR,
  "broad_group_bee_lepidoptera_fly_preVIF_gridEnvPCA_analysis"
)

STACK_DIR <- file.path(OUT_DIR, "broad_group_stacks")
FIG_DIR <- file.path(OUT_DIR, "figures")
SSDM_PNG_DIR <- file.path(FIG_DIR, "used_ssdm_png")
MAP_DIR <- file.path(FIG_DIR, "maps")
PCA_FIG_DIR <- file.path(FIG_DIR, "pca_diagnostics")
COEF_DIR <- file.path(FIG_DIR, "coefficients")
RESP_DIR <- file.path(FIG_DIR, "response_curves")
TMP_DIR <- file.path(OUT_DIR, "terra_tmp")

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(STACK_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(SSDM_PNG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(MAP_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(PCA_FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(COEF_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(RESP_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TMP_DIR, showWarnings = FALSE, recursive = TRUE)

terra::terraOptions(memfrac = 0.45, tempdir = TMP_DIR)

GRID_SIZE_DEG <- 0.5
MIN_TOTAL_SPECIES_PER_GRID <- 2
N_ENV_PC_USE <- 3
COR_FILTER_THRESHOLD <- 0.85

FORCE_REBUILD_STACKS <- FALSE

## =========================================================
## 2. Helper functions
## =========================================================

`%||%` <- function(a, b) {
  if (!is.null(a) && length(a) > 0 && !is.na(a) && nzchar(a)) a else b
}

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
  if (is.na(rhs) || rhs == "" || rhs == "1") {
    return(stats::as.formula("cbind(n_nodding_species, n_upward_species) ~ 1"))
  }
  stats::as.formula(paste("cbind(n_nodding_species, n_upward_species) ~", rhs))
}

get_rhs <- function(vars) {
  vars <- vars[!is.na(vars) & vars != ""]
  if (length(vars) == 0) return("1")
  paste(vars, collapse = " + ")
}

coef_from_glm <- function(m, model_name, method_name) {
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
    dplyr::mutate(model = model_name, method = method_name, .before = 1)
}

coef_from_glmmTMB <- function(m, model_name, method_name) {
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
    dplyr::mutate(model = model_name, method = method_name, .before = 1)
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
    error = error
  )
}

empty_coef <- function(model, method, error = NA_character_) {
  tibble::tibble(
    model = model, method = method, term = NA_character_,
    estimate = NA_real_, std.error = NA_real_,
    statistic = NA_real_, p.value = NA_real_, error = error
  )
}

safe_arrange_aic <- function(df) {
  min_aic <- suppressWarnings(min(df$AIC, na.rm = TRUE))
  if (!is.finite(min_aic)) {
    df$delta_AIC <- NA_real_
    df$akaike_weight <- NA_real_
  } else {
    df$delta_AIC <- df$AIC - min_aic
    ww <- exp(-0.5 * df$delta_AIC)
    df$akaike_weight <- ww / sum(ww, na.rm = TRUE)
  }
  df |> dplyr::arrange(.data$AIC)
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
        model = name, method = method,
        AIC = stats::AIC(m),
        BIC = stats::BIC(m),
        n_grid = stats::nobs(m),
        logLik = as.numeric(stats::logLik(m)),
        df = attr(stats::logLik(m), "df"),
        AUC_presence = safe_auc(y_presence, pred),
        logloss = logloss_binomial(data$n_nodding_species, data$n_total_species, pred),
        brier = brier_prop(data$prop_nodding_species, pred),
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
        model = name, method = method,
        AIC = stats::AIC(m),
        BIC = stats::BIC(m),
        n_grid = nrow(data),
        logLik = as.numeric(stats::logLik(m)),
        df = attr(stats::logLik(m), "df"),
        AUC_presence = safe_auc(y_presence, pred),
        logloss = logloss_binomial(data$n_nodding_species, data$n_total_species, pred),
        brier = brier_prop(data$prop_nodding_species, pred),
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

plot_model_compare <- function(df, title, out_file) {
  if (!"AIC" %in% names(df) || all(is.na(df$AIC))) return(NULL)
  plot_df <- df |>
    dplyr::filter(is.finite(.data$AIC)) |>
    dplyr::arrange(.data$AIC)
  if (nrow(plot_df) == 0) return(NULL)

  p <- plot_df |>
    dplyr::mutate(model = factor(model, levels = rev(model))) |>
    ggplot2::ggplot(ggplot2::aes(x = delta_AIC, y = model)) +
    ggplot2::geom_col() +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::labs(x = "Delta AIC", y = NULL, title = title)

  ggplot2::ggsave(out_file, p, width = 9, height = 6, dpi = 300)
  p
}

plot_ssdm_png <- function(raster_file, out_file, title = NULL, max_cells = 800000) {
  if (is.na(raster_file) || !file.exists(raster_file)) return(FALSE)

  r <- terra::rast(raster_file)
  if (terra::nlyr(r) > 1) r <- r[[1]]

  ncell_r <- terra::ncell(r)
  if (is.finite(ncell_r) && ncell_r > max_cells) {
    fact <- ceiling(sqrt(ncell_r / max_cells))
    r_plot <- terra::aggregate(r, fact = fact, fun = mean, na.rm = TRUE)
  } else {
    r_plot <- r
  }

  value_col <- names(r_plot)[1]
  clean_title <- title %||% value_col

  grDevices::png(out_file, width = 2400, height = 1700, res = 300)
  op <- graphics::par(mar = c(4.2, 4.2, 3.0, 5.2), bg = "white")
  ok <- tryCatch({
    terra::plot(
      r_plot, main = clean_title,
      xlab = "Longitude", ylab = "Latitude",
      axes = TRUE, plg = list(title = value_col)
    )
    TRUE
  }, error = function(e) {
    message("Failed to export SSDM PNG: ", raster_file, " | ", conditionMessage(e))
    FALSE
  })
  graphics::par(op)
  grDevices::dev.off()
  ok
}

## Greedy correlation filter.
## priority = earlier variables are kept preferentially.
select_by_correlation <- function(data, vars, priority, threshold = 0.85) {
  vars <- vars[vars %in% names(data)]
  priority <- unique(c(priority[priority %in% vars], setdiff(vars, priority)))

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

## =========================================================
## 3. Broad pollinator classification
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
  "Papilio", "Polygonia", "Satyrium", "Speyeria",
  "Cynthia", "Vanessa"
)

FLY_GENERA <- c(
  "Syrphus", "Eristalis", "Melanostoma", "Volucella",
  "Dasysyrphus", "Didea", "Episyrphus", "Eupeodes",
  "Ferdinandea", "Leucozona", "Meliscaeva", "Merodon",
  "Platycheirus", "Rhingia", "Sphaerophoria"
)

HAWKMOTH_GENERA <- c("Macroglossum")

## IMPORTANT:
## Hawkmoth is merged into Lepidoptera here.
BROAD_GROUPS <- list(
  bee = BEE_GENERA,
  lepidoptera = unique(c(BUTTERFLY_GENERA, HAWKMOTH_GENERA)),
  fly = FLY_GENERA
)

classify_broad_group <- function(genus) {
  dplyr::case_when(
    genus %in% BEE_GENERA ~ "bee",
    genus %in% unique(c(BUTTERFLY_GENERA, HAWKMOTH_GENERA)) ~ "lepidoptera",
    genus %in% FLY_GENERA ~ "fly",
    TRUE ~ "unknown"
  )
}

## =========================================================
## 4. Load environmental candidate stack
## =========================================================

stop_if_missing(SPECIES_SUMMARY_FILE)
stop_if_missing(CIRSIUM_OCC_FILE)

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
  stop("swe is still in candidate environment stack. Use no-SWE candidate output.")
}

message("Pre-VIF candidate environmental layers:")
print(names(env_candidate))

## =========================================================
## 5. Build merged broad-group SSDM stacks
## =========================================================

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

readr::write_csv(
  classification_summary,
  file.path(OUT_DIR, "00_pollinator_broad_group_classification_lepidoptera_merged.csv")
)

print(classification_summary, n = Inf, width = Inf)

unknowns <- classification_summary |> dplyr::filter(broad_group == "unknown")
if (nrow(unknowns) > 0) {
  warning("Some genera are unknown. Check 00_pollinator_broad_group_classification_lepidoptera_merged.csv")
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

readr::write_csv(
  stack_summary,
  file.path(OUT_DIR, "01_broad_group_ssdm_stack_summary.csv")
)

## Export SSDM PNGs
png_records <- list()
png_counter <- 1

for (i in seq_len(nrow(stack_summary))) {
  rr <- stack_summary[i, ]

  raster_items <- tibble::tibble(
    broad_group = rr$broad_group,
    raster_type = c("sum_cloglog", "mean_cloglog", "binary_richness_10p"),
    raster_file = c(rr$sum_cloglog_file, rr$mean_cloglog_file, rr$binary_richness_file)
  ) |>
    dplyr::filter(!is.na(raster_file), file.exists(raster_file))

  for (j in seq_len(nrow(raster_items))) {
    item <- raster_items[j, ]

    out_png <- file.path(
      SSDM_PNG_DIR,
      paste0("used_ssdm_", safe_name(item$broad_group), "_", item$raster_type, ".png")
    )

    pretty_group <- dplyr::recode(
      as.character(item$broad_group),
      bee = "Bee",
      lepidoptera = "Lepidoptera",
      fly = "Fly",
      .default = as.character(item$broad_group)
    )

    pretty_type <- dplyr::recode(
      as.character(item$raster_type),
      sum_cloglog = "summed suitability",
      mean_cloglog = "mean suitability",
      binary_richness_10p = "binary richness",
      .default = as.character(item$raster_type)
    )

    ok <- plot_ssdm_png(
      raster_file = item$raster_file,
      out_file = out_png,
      title = paste0(pretty_group, " SSDM: ", pretty_type)
    )

    png_records[[png_counter]] <- tibble::tibble(
      broad_group = item$broad_group,
      raster_type = item$raster_type,
      raster_file = item$raster_file,
      png_file = out_png,
      exported = ok
    )
    png_counter <- png_counter + 1
  }
}

png_tbl <- dplyr::bind_rows(png_records)
readr::write_csv(png_tbl, file.path(OUT_DIR, "01b_used_ssdm_png_files.csv"))

## =========================================================
## 6. Load Cirsium occurrences and extract env + SSDM
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

xy <- cbind(occ0$decimalLongitude, occ0$decimalLatitude)
colnames(xy) <- c("lon", "lat")

env_vals <- safe_extract_df(env_candidate, xy)
names(env_vals) <- paste0("gridPCAenv_", names(env_candidate))

all_stack_files <- c(
  stack_summary$sum_cloglog_file,
  stack_summary$mean_cloglog_file,
  stack_summary$binary_richness_file
)
all_stack_files <- all_stack_files[!is.na(all_stack_files) & file.exists(all_stack_files)]

if (length(all_stack_files) == 0) stop("No broad SSDM rasters found.")

broad_r <- terra::rast(all_stack_files)
names(broad_r) <- safe_name(tools::file_path_sans_ext(basename(all_stack_files)))

broad_vals <- safe_extract_df(broad_r, xy)
names(broad_vals) <- names(broad_r)

occ_all <- dplyr::bind_cols(occ0, env_vals, broad_vals, .name_repair = "unique")

readr::write_csv(
  occ_all,
  file.path(OUT_DIR, "02_cirsium_occurrences_with_preVIF_env_and_broad_SSDM.csv")
)

## =========================================================
## 7. Build grid table
## =========================================================

occ_grid <- occ_all |>
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

env_cols <- names(occ_grid)[stringr::str_detect(names(occ_grid), "^gridPCAenv_")]
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

readr::write_csv(
  grid_dat,
  file.path(OUT_DIR, "03_grid_raw_counts_preVIF_env_broad_SSDM.csv")
)

## =========================================================
## 8. Rebuild environmental PCA at GRID level from pre-VIF env variables
## =========================================================

env_grid <- grid_dat |>
  dplyr::select(dplyr::all_of(env_cols)) |>
  dplyr::select(where(is.numeric))

good_env_cols <- names(env_grid)[
  vapply(env_grid, function(x) {
    sx <- stats::sd(x, na.rm = TRUE)
    sum(is.finite(x)) >= 5 && is.finite(sx) && sx > 0
  }, logical(1))
]

env_grid <- env_grid[, good_env_cols, drop = FALSE]

if (ncol(env_grid) < 2) {
  stop("Too few pre-VIF environmental variables for grid PCA.")
}

env_grid_imp <- env_grid
for (nm in names(env_grid_imp)) {
  x <- env_grid_imp[[nm]]
  x[!is.finite(x)] <- stats::median(x, na.rm = TRUE)
  env_grid_imp[[nm]] <- x
}

env_pca <- stats::prcomp(env_grid_imp, center = TRUE, scale. = TRUE)
saveRDS(env_pca, file.path(OUT_DIR, "04_env_gridPCA_preVIF_model.rds"))

env_scores <- as.data.frame(env_pca$x) |>
  tibble::as_tibble()

pc_keep <- intersect(paste0("PC", seq_len(N_ENV_PC_USE)), names(env_scores))
env_scores_keep <- env_scores[, pc_keep, drop = FALSE]
names(env_scores_keep) <- paste0("env_", names(env_scores_keep))

grid_dat <- dplyr::bind_cols(grid_dat, env_scores_keep)

for (nm in names(env_scores_keep)) {
  grid_dat[[paste0("z_", nm)]] <- safe_scale(grid_dat[[nm]])
}

env_loadings <- as.data.frame(env_pca$rotation) |>
  tibble::rownames_to_column("variable") |>
  tibble::as_tibble()

env_var_tbl <- tibble::tibble(
  PC = paste0("PC", seq_along(env_pca$sdev)),
  variance = env_pca$sdev^2,
  prop_variance = variance / sum(variance),
  cum_variance = cumsum(prop_variance)
)

env_top_loadings <- env_loadings |>
  tidyr::pivot_longer(cols = dplyr::starts_with("PC"), names_to = "PC", values_to = "loading") |>
  dplyr::filter(PC %in% paste0("PC", 1:6)) |>
  dplyr::group_by(PC) |>
  dplyr::slice_max(order_by = abs(loading), n = 12, with_ties = FALSE) |>
  dplyr::arrange(PC, dplyr::desc(abs(loading))) |>
  dplyr::ungroup()

readr::write_csv(env_loadings, file.path(OUT_DIR, "04_env_gridPCA_preVIF_loadings.csv"))
readr::write_csv(env_var_tbl, file.path(OUT_DIR, "04_env_gridPCA_preVIF_variance_explained.csv"))
readr::write_csv(env_top_loadings, file.path(OUT_DIR, "04b_env_gridPCA_preVIF_top_loadings_PC1_PC6.csv"))

p_var <- env_var_tbl |>
  dplyr::slice_head(n = 10) |>
  ggplot2::ggplot(ggplot2::aes(x = PC, y = prop_variance)) +
  ggplot2::geom_col() +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::labs(
    x = "Environmental PC",
    y = "Proportion of variance",
    title = "Pre-VIF environmental grid PCA"
  )

ggplot2::ggsave(
  file.path(PCA_FIG_DIR, "env_gridPCA_preVIF_variance_explained.png"),
  p_var, width = 7, height = 5, dpi = 300
)

p_load <- env_top_loadings |>
  dplyr::filter(PC %in% paste0("PC", 1:3)) |>
  dplyr::mutate(variable = stringr::str_remove(variable, "^gridPCAenv_")) |>
  ggplot2::ggplot(ggplot2::aes(x = loading, y = forcats::fct_reorder(variable, loading))) +
  ggplot2::geom_vline(xintercept = 0, linetype = 2) +
  ggplot2::geom_col() +
  ggplot2::facet_wrap(~ PC, scales = "free_y") +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::labs(
    x = "PCA loading",
    y = NULL,
    title = "Top loadings: pre-VIF environmental grid PCA"
  )

ggplot2::ggsave(
  file.path(PCA_FIG_DIR, "env_gridPCA_preVIF_PC1_PC3_top_loadings.png"),
  p_load, width = 12, height = 8, dpi = 300
)

## =========================================================
## 9. Scale pollinator variables and build analysis grid
## =========================================================

## Use only sum_cloglog for main analysis.
POLL_VARS <- c(
  "broad_bee_sum_cloglog",
  "broad_lepidoptera_sum_cloglog",
  "broad_fly_sum_cloglog"
)
POLL_VARS <- POLL_VARS[POLL_VARS %in% names(grid_dat)]

for (nm in POLL_VARS) {
  grid_dat[[nm]][is.na(grid_dat[[nm]])] <- 0
  grid_dat[[paste0("z_", nm)]] <- safe_scale(grid_dat[[nm]])
}

analysis_grid <- grid_dat |>
  dplyr::filter(n_total_species >= MIN_TOTAL_SPECIES_PER_GRID) |>
  tidyr::drop_na(
    n_nodding_species, n_upward_species, n_total_species,
    dplyr::all_of(paste0("z_", names(env_scores_keep)))
  ) |>
  dplyr::mutate(
    prop_nodding_species = n_nodding_species / n_total_species,
    y_presence = ifelse(n_nodding_species > 0, 1, 0),
    z_lon = safe_scale(grid_lon),
    z_lat = safe_scale(grid_lat)
  )

readr::write_csv(
  analysis_grid,
  file.path(OUT_DIR, "05_analysis_grid_preVIF_gridEnvPCA_lepidoptera_merged.csv")
)

message("\nAnalysis grid summary:")
print(analysis_grid |> dplyr::summarise(
  n_grid = dplyr::n(),
  mean_total_species = mean(n_total_species),
  total_nodding_species_counts = sum(n_nodding_species),
  total_upward_species_counts = sum(n_upward_species),
  mean_prop_nodding = mean(prop_nodding_species)
))

## =========================================================
## 10. Correlation diagnostics and filtered predictor sets
## =========================================================

ENV_VARS <- paste0("z_", names(env_scores_keep))
ENV_VARS <- ENV_VARS[ENV_VARS %in% names(analysis_grid)]
POLL_Z_VARS <- paste0("z_", POLL_VARS)
POLL_Z_VARS <- POLL_Z_VARS[POLL_Z_VARS %in% names(analysis_grid)]

DIAG_VARS <- unique(c(
  ENV_VARS,
  POLL_Z_VARS,
  "z_lon", "z_lat",
  "prop_nodding_species"
))

cor_mat <- stats::cor(
  analysis_grid[, DIAG_VARS],
  use = "pairwise.complete.obs"
)

cor_pairs <- as.data.frame(as.table(cor_mat)) |>
  tibble::as_tibble() |>
  dplyr::rename(var1 = Var1, var2 = Var2, correlation = Freq) |>
  dplyr::filter(as.character(var1) < as.character(var2)) |>
  dplyr::arrange(dplyr::desc(abs(correlation)))

readr::write_csv(cor_pairs, file.path(OUT_DIR, "07_predictor_response_correlation_pairs.csv"))

cor_predictor_pairs <- cor_pairs |>
  dplyr::filter(var1 != "prop_nodding_species", var2 != "prop_nodding_species")

readr::write_csv(
  cor_predictor_pairs,
  file.path(OUT_DIR, "07_predictor_only_correlation_pairs.csv")
)

message("\nStrong correlations among predictors/response:")
print(cor_pairs |> dplyr::filter(abs(correlation) >= 0.60), n = Inf, width = Inf)

## Heatmap
cor_df <- as.data.frame(as.table(cor_mat)) |>
  tibble::as_tibble() |>
  dplyr::rename(var1 = Var1, var2 = Var2, correlation = Freq)

p_cor <- ggplot2::ggplot(cor_df, ggplot2::aes(var1, var2, fill = correlation)) +
  ggplot2::geom_tile() +
  ggplot2::geom_text(ggplot2::aes(label = round(correlation, 2)), size = 2.8) +
  ggplot2::scale_fill_gradient2(limits = c(-1, 1), midpoint = 0) +
  ggplot2::coord_equal() +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
  ggplot2::labs(title = "Predictor correlation: pre-VIF grid env PCA + broad SSDM", x = NULL, y = NULL)

ggplot2::ggsave(
  file.path(FIG_DIR, "predictor_correlation_heatmap.png"),
  p_cor, width = 10, height = 9, dpi = 300
)

## Correlation-filtered predictor sets
## Option A: keep environmental PCs preferentially.
priority_keep_env <- c(
  POLL_Z_VARS,
  "z_env_PC1", "z_env_PC2", "z_env_PC3",
  "z_lat", "z_lon"
)

## Option B: keep latitude as spatial proxy preferentially when it conflicts with env PCs.
priority_keep_lat <- c(
  POLL_Z_VARS,
  "z_lat",
  "z_env_PC1", "z_env_PC3", "z_env_PC2",
  "z_lon"
)

all_predictors_for_filter <- unique(c(ENV_VARS, POLL_Z_VARS, "z_lon", "z_lat"))

sel_env <- select_by_correlation(
  analysis_grid, all_predictors_for_filter,
  priority = priority_keep_env,
  threshold = COR_FILTER_THRESHOLD
)

sel_lat <- select_by_correlation(
  analysis_grid, all_predictors_for_filter,
  priority = priority_keep_lat,
  threshold = COR_FILTER_THRESHOLD
)

filter_decisions <- dplyr::bind_rows(
  sel_env$dropped |> dplyr::mutate(filter_set = "corrfiltered_keep_envPC", .before = 1),
  sel_lat$dropped |> dplyr::mutate(filter_set = "corrfiltered_keep_lat", .before = 1)
)

filter_retained <- dplyr::bind_rows(
  tibble::tibble(filter_set = "corrfiltered_keep_envPC", variable = sel_env$selected),
  tibble::tibble(filter_set = "corrfiltered_keep_lat", variable = sel_lat$selected)
)

readr::write_csv(filter_decisions, file.path(OUT_DIR, "07b_correlation_filter_decisions.csv"))
readr::write_csv(filter_retained, file.path(OUT_DIR, "07c_correlation_filter_retained_predictors.csv"))

## =========================================================
## 11. Model specifications
## =========================================================

rhs_env <- get_rhs(ENV_VARS)
rhs_poll <- get_rhs(POLL_Z_VARS)
rhs_env_poll <- get_rhs(c(ENV_VARS, POLL_Z_VARS))
rhs_lat <- "z_lat"
rhs_lonlat <- "z_lon + z_lat"
rhs_space_poly <- "z_lon + z_lat + I(z_lon^2) + I(z_lat^2) + z_lon:z_lat"

model_specs <- list(
  null = "1",
  lat_only = rhs_lat,
  lonlat = rhs_lonlat,
  space_poly = rhs_space_poly,

  env = rhs_env,
  poll_lepidoptera_merged = rhs_poll,
  env_plus_poll = rhs_env_poll,

  env_plus_lat = get_rhs(c(ENV_VARS, "z_lat")),
  poll_plus_lat = get_rhs(c(POLL_Z_VARS, "z_lat")),
  env_plus_poll_plus_lat = get_rhs(c(ENV_VARS, POLL_Z_VARS, "z_lat")),

  env_plus_space_poly = paste(rhs_env, rhs_space_poly, sep = " + "),
  poll_plus_space_poly = paste(rhs_poll, rhs_space_poly, sep = " + "),
  env_plus_poll_plus_space_poly = paste(rhs_env_poll, rhs_space_poly, sep = " + "),

  corrfiltered_keep_envPC = get_rhs(sel_env$selected),
  corrfiltered_keep_lat = get_rhs(sel_lat$selected)
)

model_specs <- model_specs[
  !vapply(model_specs, function(x) is.na(x) || x == "", logical(1))
]

readr::write_csv(
  tibble::tibble(model = names(model_specs), rhs = unlist(model_specs)),
  file.path(OUT_DIR, "08_model_specifications.csv")
)

## =========================================================
## 12. Fit GLM and beta-binomial
## =========================================================

glm_results <- list()
bb_results <- list()

for (nm in names(model_specs)) {
  f <- make_formula(model_specs[[nm]])
  message("Fitting: ", nm)
  glm_results[[nm]] <- fit_glm_binom(f, analysis_grid, nm)
  bb_results[[nm]] <- fit_glmmTMB_betabinom(f, analysis_grid, nm)
}

glm_summary <- dplyr::bind_rows(lapply(glm_results, `[[`, "summary")) |> safe_arrange_aic()
bb_summary <- dplyr::bind_rows(lapply(bb_results, `[[`, "summary")) |> safe_arrange_aic()

glm_coef <- dplyr::bind_rows(lapply(glm_results, `[[`, "coef"))
bb_coef <- dplyr::bind_rows(lapply(bb_results, `[[`, "coef"))

readr::write_csv(glm_summary, file.path(OUT_DIR, "09_binomial_GLM_model_comparison.csv"))
readr::write_csv(glm_coef, file.path(OUT_DIR, "09_binomial_GLM_coefficients.csv"))
readr::write_csv(bb_summary, file.path(OUT_DIR, "10_beta_binomial_model_comparison.csv"))
readr::write_csv(bb_coef, file.path(OUT_DIR, "10_beta_binomial_coefficients.csv"))

message("\nBinomial GLM model comparison:")
print(glm_summary, n = Inf, width = Inf)

message("\nBeta-binomial model comparison:")
print(bb_summary, n = Inf, width = Inf)

plot_model_compare(
  glm_summary,
  "Binomial GLM model comparison",
  file.path(FIG_DIR, "binomial_GLM_delta_AIC.png")
)

plot_model_compare(
  bb_summary,
  "Beta-binomial model comparison",
  file.path(FIG_DIR, "beta_binomial_delta_AIC.png")
)

## =========================================================
## 13. Best model coefficients and predictions
## =========================================================

best_method <- "beta_binomial_glmmTMB"
best_table <- bb_summary

if (all(is.na(best_table$AIC))) {
  best_method <- "binomial_GLM"
  best_table <- glm_summary
}

best_table_nonNA <- best_table |> dplyr::filter(is.finite(.data$AIC))

if (nrow(best_table_nonNA) > 0) {
  best_name <- best_table_nonNA$model[1]
} else {
  best_method <- "binomial_GLM"
  best_table_nonNA <- glm_summary |> dplyr::filter(is.finite(.data$AIC))
  best_name <- best_table_nonNA$model[1]
}

best_model_object <- NULL
best_coef <- NULL

if (best_method == "beta_binomial_glmmTMB") {
  best_model_object <- bb_results[[best_name]]$model
  best_coef <- bb_coef |> dplyr::filter(model == best_name)
} else {
  best_model_object <- glm_results[[best_name]]$model
  best_coef <- glm_coef |> dplyr::filter(model == best_name)
}

readr::write_csv(best_coef, file.path(OUT_DIR, "11_best_model_coefficients.csv"))

message("\nBest model:")
print(best_table_nonNA[1, ], width = Inf)

message("\nBest coefficients:")
print(best_coef, n = Inf, width = Inf)

## Coefficient plot
coef_plot <- best_coef |>
  dplyr::filter(term != "(Intercept)") |>
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

  ggplot2::ggsave(
    file.path(COEF_DIR, "best_model_coefficients.png"),
    p_coef, width = 10, height = 6, dpi = 300
  )
}

## Predictions
if (!is.null(best_model_object)) {
  saveRDS(
    best_model_object,
    file.path(OUT_DIR, paste0("11_best_model_object_", best_method, "_", best_name, ".rds"))
  )

  analysis_grid$pred_prob_nodding_best <- tryCatch(
    as.numeric(stats::predict(best_model_object, newdata = analysis_grid, type = "response")),
    error = function(e) {
      warning("Could not predict best model probability: ", conditionMessage(e))
      rep(NA_real_, nrow(analysis_grid))
    }
  )

  analysis_grid$pred_prob_upward_best <- 1 - analysis_grid$pred_prob_nodding_best
  analysis_grid$resid_prop_nodding_best <- analysis_grid$prop_nodding_species - analysis_grid$pred_prob_nodding_best
  analysis_grid$best_model_name <- best_name
  analysis_grid$best_model_method <- best_method

  readr::write_csv(
    analysis_grid,
    file.path(OUT_DIR, "12_predicted_probability_nodding_grid_best_model.csv")
  )
}

## Maps
p_obs <- ggplot2::ggplot(
  analysis_grid,
  ggplot2::aes(grid_lon, grid_lat, color = prop_nodding_species, size = n_total_species)
) +
  ggplot2::geom_point(alpha = 0.85) +
  ggplot2::coord_equal() +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::scale_color_viridis_c(limits = c(0, 1)) +
  ggplot2::labs(
    x = "Longitude", y = "Latitude",
    color = "Observed\nproportion",
    size = "Total\nspecies",
    title = "Observed proportion of nodding Cirsium species"
  )

ggplot2::ggsave(file.path(MAP_DIR, "map_observed_prop_nodding.png"), p_obs, width = 7, height = 7, dpi = 300)

if ("pred_prob_nodding_best" %in% names(analysis_grid)) {
  p_pred <- ggplot2::ggplot(
    analysis_grid,
    ggplot2::aes(grid_lon, grid_lat, color = pred_prob_nodding_best, size = n_total_species)
  ) +
    ggplot2::geom_point(alpha = 0.9) +
    ggplot2::coord_equal() +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::scale_color_viridis_c(limits = c(0, 1)) +
    ggplot2::labs(
      x = "Longitude", y = "Latitude",
      color = "Predicted\nP(nodding)",
      size = "Total\nspecies",
      title = "Predicted proportion of nodding Cirsium species"
    )

  ggplot2::ggsave(file.path(MAP_DIR, "map_predicted_probability_nodding_best_model.png"), p_pred, width = 7, height = 7, dpi = 300)

  p_resid <- ggplot2::ggplot(
    analysis_grid,
    ggplot2::aes(grid_lon, grid_lat, color = resid_prop_nodding_best, size = n_total_species)
  ) +
    ggplot2::geom_point(alpha = 0.9) +
    ggplot2::coord_equal() +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::scale_color_gradient2(midpoint = 0) +
    ggplot2::labs(
      x = "Longitude", y = "Latitude",
      color = "Observed -\npredicted",
      size = "Total\nspecies",
      title = "Residual proportion of nodding species"
    )

  ggplot2::ggsave(file.path(MAP_DIR, "map_residual_probability_nodding_best_model.png"), p_resid, width = 7, height = 7, dpi = 300)
}

## Response curves for best model variables
if (!is.null(best_model_object) && nrow(coef_plot) > 0) {
  response_terms <- coef_plot |> dplyr::pull(term)
  response_terms <- response_terms[response_terms %in% names(analysis_grid)]

  for (v in response_terms) {
    grid_new <- analysis_grid[rep(1, 100), , drop = FALSE]

    for (nm in names(grid_new)) {
      if (is.numeric(grid_new[[nm]])) {
        grid_new[[nm]] <- stats::median(analysis_grid[[nm]], na.rm = TRUE)
      }
    }

    grid_new[[v]] <- seq(
      min(analysis_grid[[v]], na.rm = TRUE),
      max(analysis_grid[[v]], na.rm = TRUE),
      length.out = 100
    )

    grid_new$pred <- tryCatch(
      as.numeric(stats::predict(best_model_object, newdata = grid_new, type = "response")),
      error = function(e) rep(NA_real_, nrow(grid_new))
    )

    p_resp <- ggplot2::ggplot(grid_new, ggplot2::aes(x = .data[[v]], y = pred)) +
      ggplot2::geom_line(linewidth = 1) +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::labs(
        x = v,
        y = "Predicted P(nodding)",
        title = paste0("Response curve: ", v)
      )

    ggplot2::ggsave(
      file.path(RESP_DIR, paste0("response_curve_", safe_name(v), ".png")),
      p_resp, width = 6, height = 4, dpi = 300
    )
  }
}

## =========================================================
## 14. Done
## =========================================================

message("\n🎉 DONE — Lepidoptera-merged broad SSDM + pre-VIF GRID environmental PCA analysis")
message("Output dir: ", OUT_DIR)
message("Key outputs:")
message("  04_env_gridPCA_preVIF_variance_explained.csv")
message("  04_env_gridPCA_preVIF_loadings.csv")
message("  04b_env_gridPCA_preVIF_top_loadings_PC1_PC6.csv")
message("  05_analysis_grid_preVIF_gridEnvPCA_lepidoptera_merged.csv")
message("  07_predictor_response_correlation_pairs.csv")
message("  07b_correlation_filter_decisions.csv")
message("  07c_correlation_filter_retained_predictors.csv")
message("  09_binomial_GLM_model_comparison.csv")
message("  10_beta_binomial_model_comparison.csv")
message("  11_best_model_coefficients.csv")
message("  12_predicted_probability_nodding_grid_best_model.csv")
message("Figures:")
message("  figures/used_ssdm_png/*.png")
message("  figures/pca_diagnostics/*.png")
message("  figures/beta_binomial_delta_AIC.png")
message("  figures/coefficients/best_model_coefficients.png")
message("  figures/maps/*.png")
message("  figures/response_curves/*.png")
