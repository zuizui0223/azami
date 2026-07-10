#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ape)
  library(jsonlite)
  library(rotl)
  library(V.PhyloMaker2)
})

args <- commandArgs(trailingOnly = TRUE)
value_after <- function(flag, default = NULL) {
  index <- match(flag, args)
  if (is.na(index)) return(default)
  if (index == length(args)) stop(sprintf("Missing value after %s", flag))
  args[[index + 1]]
}

species_path <- value_after("--species")
out_dir <- value_after("--out-dir")
n_random <- as.integer(value_after("--n-random", "50"))
seed <- as.integer(value_after("--seed", "20260710"))
if (is.null(species_path) || is.null(out_dir)) {
  stop("Required: --species and --out-dir")
}
if (!is.finite(n_random) || n_random < 1L) stop("--n-random must be positive")
set.seed(seed)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

species <- read.csv(species_path, stringsAsFactors = FALSE, check.names = FALSE)
if (!"taxon_name" %in% names(species)) stop("Species table lacks taxon_name")
taxa <- sort(unique(trimws(species$taxon_name)))
taxa <- taxa[nzchar(taxa)]
if (length(taxa) < 30L) stop("Too few taxa for historical sensitivity")
if (anyDuplicated(taxa)) stop("Taxon names are not unique")

# ---------------------------------------------------------------------------
# Open Tree of Life: name-resolution and synthetic-topology audit.
# This tree is not assigned biological branch lengths and is not the primary
# PGLS tree. It is an independent topology/provenance sensitivity resource.
# ---------------------------------------------------------------------------
open_tree_error <- ""
open_matches <- tryCatch(
  rotl::tnrs_match_names(
    names = taxa,
    context_name = "Asteraceae",
    do_approximate_matching = TRUE
  ),
  error = function(error) {
    open_tree_error <<- sprintf("%s: %s", class(error)[1], conditionMessage(error))
    NULL
  }
)

open_summary <- list(
  n_input_taxa = length(taxa),
  n_rows = 0L,
  n_resolved = 0L,
  n_exact_score_one = 0L,
  n_approximate = 0L,
  n_synonym = 0L,
  n_ambiguous_multiple_matches = 0L,
  n_unique_ott_ids = 0L,
  induced_tree_tips = 0L,
  error = open_tree_error
)

if (!is.null(open_matches)) {
  write.csv(
    open_matches,
    file.path(out_dir, "opentree_name_matches.csv"),
    row.names = FALSE,
    na = ""
  )
  resolved <- open_matches[!is.na(open_matches$ott_id), , drop = FALSE]
  open_summary$n_rows <- nrow(open_matches)
  open_summary$n_resolved <- nrow(resolved)
  open_summary$n_exact_score_one <- sum(
    !is.na(resolved$score) & resolved$score >= 0.999999 &
      !resolved$approximate_match
  )
  open_summary$n_approximate <- sum(resolved$approximate_match, na.rm = TRUE)
  open_summary$n_synonym <- sum(resolved$is_synonym, na.rm = TRUE)
  open_summary$n_ambiguous_multiple_matches <- sum(
    !is.na(resolved$number_matches) & resolved$number_matches > 1,
    na.rm = TRUE
  )
  unique_ott <- unique(resolved$ott_id)
  open_summary$n_unique_ott_ids <- length(unique_ott)
  if (length(unique_ott) >= 2L) {
    open_tree <- tryCatch(
      rotl::tol_induced_subtree(ott_ids = unique_ott),
      error = function(error) {
        open_summary$error <<- sprintf(
          "tol_induced_subtree: %s: %s",
          class(error)[1], conditionMessage(error)
        )
        NULL
      }
    )
    if (!is.null(open_tree)) {
      open_summary$induced_tree_tips <- length(open_tree$tip.label)
      write.tree(open_tree, file.path(out_dir, "opentree_synthetic_induced_subtree.tre"))
      open_tree_tip_table <- data.frame(
        tip_label = open_tree$tip.label,
        stringsAsFactors = FALSE
      )
      write.csv(
        open_tree_tip_table,
        file.path(out_dir, "opentree_induced_tip_labels.csv"),
        row.names = FALSE
      )
    }
  }
}

open_about <- tryCatch(rotl::tol_about(), error = function(error) NULL)
if (!is.null(open_about)) {
  write_json(
    open_about,
    file.path(out_dir, "opentree_synthesis_about.json"),
    pretty = TRUE,
    auto_unbox = TRUE,
    null = "null"
  )
}

# ---------------------------------------------------------------------------
# GBOTB/V.PhyloMaker2: dated backbone plus explicit grafting scenarios.
# Direct backbone representation is audited before any PGLS is run.
# IMPORTANT: output.tree=FALSE returns only the pruned user-species hypotheses.
# output.tree=TRUE would retain one full ~73,000-tip backbone per replicate and
# is neither necessary nor computationally appropriate here.
# ---------------------------------------------------------------------------
sp_list <- data.frame(
  species = taxa,
  genus = sub(" .*", "", taxa),
  family = rep("Asteraceae", length(taxa)),
  stringsAsFactors = FALSE
)
canonical_tip <- gsub(" ", "_", taxa, fixed = TRUE)
in_backbone <- canonical_tip %in% GBOTB.extended.LCVP$tip.label
backbone_audit <- data.frame(
  taxon_name = taxa,
  canonical_tip = canonical_tip,
  genus = sp_list$genus,
  family = sp_list$family,
  directly_present_in_GBOTB_extended_LCVP = in_backbone,
  historical_information_tier = ifelse(
    in_backbone,
    "direct_backbone_tip",
    "grafted_within_genus_hypothesis"
  ),
  stringsAsFactors = FALSE
)
write.csv(
  backbone_audit,
  file.path(out_dir, "gbotb_lcvp_backbone_audit.csv"),
  row.names = FALSE
)

maker <- V.PhyloMaker2::phylo.maker(
  sp.list = sp_list,
  tree = V.PhyloMaker2::GBOTB.extended.LCVP,
  nodes = V.PhyloMaker2::nodes.info.1.LCVP,
  output.tree = FALSE,
  scenarios = c("S1", "S2", "S3"),
  r = n_random
)

scenario1 <- maker$scenario.1
scenario2 <- maker$scenario.2
scenario3 <- maker$scenario.3
if (is.null(scenario1) || is.null(scenario2) || is.null(scenario3)) {
  stop("V.PhyloMaker2 did not return all requested pruned scenarios")
}
if (!inherits(scenario2, "multiPhylo")) class(scenario2) <- "multiPhylo"
expected <- sort(canonical_tip)
for (label in c("scenario1", "scenario3")) {
  tree <- get(label)
  if (!identical(sort(tree$tip.label), expected)) {
    stop(sprintf("%s tip labels do not match all input taxa", label))
  }
}
if (length(scenario2) != n_random) stop("Unexpected number of randomized trees")
if (any(vapply(
  scenario2,
  function(tree) !identical(sort(tree$tip.label), expected),
  logical(1)
))) stop("At least one randomized tree does not match all taxa")

write.tree(scenario1, file.path(out_dir, "gbotb_lcvp_scenario1.tre"))
write.tree(scenario3, file.path(out_dir, "gbotb_lcvp_scenario3.tre"))
write.tree(scenario2, file.path(out_dir, "gbotb_lcvp_scenario2_randomized.trees"))
if (!is.null(maker$species.list)) {
  write.csv(
    maker$species.list,
    file.path(out_dir, "vphylomaker_species_list_output.csv"),
    row.names = FALSE,
    na = ""
  )
}

edge_audit <- function(tree, scenario, replicate = NA_integer_) {
  positive <- tree$edge.length[is.finite(tree$edge.length) & tree$edge.length > 0]
  data.frame(
    scenario = scenario,
    replicate = replicate,
    n_tips = length(tree$tip.label),
    n_nodes = tree$Nnode,
    ultrametric = isTRUE(ape::is.ultrametric(tree)),
    binary = isTRUE(ape::is.binary.tree(tree)),
    n_zero_or_negative_edges = sum(
      !is.finite(tree$edge.length) | tree$edge.length <= 0
    ),
    min_positive_edge = if (length(positive)) min(positive) else NA_real_,
    max_edge = if (length(positive)) max(positive) else NA_real_,
    stringsAsFactors = FALSE
  )
}
edge_rows <- rbind(
  edge_audit(scenario1, "S1"),
  edge_audit(scenario3, "S3"),
  do.call(rbind, lapply(seq_along(scenario2), function(index) {
    edge_audit(scenario2[[index]], "S2_random", index)
  }))
)
write.csv(
  edge_rows,
  file.path(out_dir, "historical_tree_edge_audit.csv"),
  row.names = FALSE
)

report <- list(
  created_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  seed = seed,
  n_input_taxa = length(taxa),
  opentree = open_summary,
  gbotb_lcvp = list(
    n_direct_backbone_tips = sum(in_backbone),
    direct_backbone_fraction = mean(in_backbone),
    n_grafted_within_genus = sum(!in_backbone),
    n_randomized_scenario2_trees = length(scenario2),
    scenario1_tips = length(scenario1$tip.label),
    scenario3_tips = length(scenario3$tip.label),
    all_trees_have_nonnegative_branch_lengths = all(
      edge_rows$n_zero_or_negative_edges == 0
    )
  ),
  interpretation = list(
    opentree_role = paste(
      "Independent synthetic topology and taxonomic-name audit; branch lengths",
      "are not treated as a dated Cirsium species history."
    ),
    gbotb_role = paste(
      "Dated vascular-plant backbone sensitivity. Directly represented tips",
      "and genus-grafted hypotheses are reported separately."
    ),
    biological_boundary = paste(
      "Neither resource resolves reticulate evolution, chloroplast capture or",
      "allopolyploid parentage. PGLS results are historical-constraint",
      "sensitivities rather than proof that one tree is the true species tree."
    )
  )
)
write_json(
  report,
  file.path(out_dir, "historical_tree_audit_report.json"),
  pretty = TRUE,
  auto_unbox = TRUE,
  null = "null"
)
cat(toJSON(report, pretty = TRUE, auto_unbox = TRUE), "\n")
