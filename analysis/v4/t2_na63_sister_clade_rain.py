"""T2 descriptive sister-clade rainfall contrasts on the author-published NA63 nuclear tree.

The four pure-nodding tip clades are NOT necessarily four independently
established gains. Modern BIO12 niche differences are NOT ancestral rain.
This complements, and cannot overrule, the existing NA63 phylogenetic Mk null.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from analysis.v4.t2_pgls_placement_screen import (
    AUTHOR_NUCLEAR_TREE_GIT_BLOB, AUTHOR_NUCLEAR_TREE_PATH,
    SOURCE_DOI, SOURCE_COMMIT, SOURCE_REPO,
    download_author_states, read_author_nuclear_tree, parsimony,
)
from analysis.v4.t2_phylogenetic_climate_null import (
    ENV_PATH, ENV_GIT_BLOB, fetch_author, read_published_climate,
)

EXPECTED_N = 63
EXPECTED_DOWN = 5
EXPECTED_PURE_DOWN_CLADES = 4


def sister_clade_contrasts(tree, states, climate):
    """All-nodding maximal tip clades vs their immediate all-erect sisters."""
    tip_names = [x.name for x in tree.get_terminals()]
    if (len(tip_names) != EXPECTED_N or len(set(tip_names)) != EXPECTED_N
            or set(tip_names) != set(states) or set(tip_names) != set(climate)
            or sum(states.values()) != EXPECTED_DOWN):
        raise ValueError("Published 63-tip tree/state/climate source contract changed")

    parents = {c: p for p in tree.find_clades() for c in p.clades}
    descendants = {}
    for n in tree.find_clades(order="postorder"):
        descendants[n] = ([n.name] if n.is_terminal()
                          else [x for c in n.clades for x in descendants[c]])
    all_down = {n: all(states[x] == 1 for x in xs) for n, xs in descendants.items()}
    groups = [n for n in tree.find_clades() if all_down[n]
              and (n not in parents or not all_down[parents[n]])]
    if len(groups) != EXPECTED_PURE_DOWN_CLADES:
        raise ValueError("Four author-rooted pure-nodding tip clades not found")

    def mean(taxa, variable):
        return sum(climate[t][variable] for t in taxa) / len(taxa)

    pairs = []
    for node in groups:
        parent = parents.get(node)
        if parent is None or len(parent.clades) != 2:
            raise ValueError("No auditable rooted immediate sister clade")
        sister = parent.clades[0] if parent.clades[1] is node else parent.clades[1]
        down = sorted(descendants[node])
        up = sorted(descendants[sister])
        if any(states[t] != 0 for t in up):
            raise ValueError("Immediate sister includes nodding; cannot report pure comparison")
        d12, u12 = mean(down, "BIO12"), mean(up, "BIO12")
        pairs.append({
            "nodding_taxa": down,
            "sister_erect_taxa": up,
            "n_nodding_taxa": len(down),
            "n_sister_erect_taxa": len(up),
            "nodding_BIO12_mean_mm": d12,
            "sister_erect_BIO12_mean_mm": u12,
            "nodding_minus_sister_erect_BIO12_mm": d12-u12,
            "nodding_BIO1_original_source_unit_mean": mean(down, "BIO1"),
            "sister_erect_BIO1_original_source_unit_mean": mean(up, "BIO1"),
        })
    pairs.sort(key=lambda r: r["nodding_taxa"][0])
    pars = parsimony(tree, states)
    n_wet = sum(r["nodding_minus_sister_erect_BIO12_mm"] > 0 for r in pairs)
    n_dry = sum(r["nodding_minus_sister_erect_BIO12_mm"] < 0 for r in pairs)
    if (pars["root_erect_min_steps"] != 4 or
            pars["root_nodding_min_steps"] != 6 or
            n_wet != 2 or n_dry != 2):
        raise ValueError("Original published data comparison drifted; do not promote claims")
    return {
        "scientific_status": "DESCRIPTIVE_SISTER_CLADE_SIGNS_NOT_ADAPTIVE_EVOLUTION",
        "n_author_tips": len(tip_names),
        "n_nodding_author_scored": EXPECTED_DOWN,
        "n_maximal_pure_nodding_tip_clades": len(groups),
        "root_erect_minimum_total_changes": pars["root_erect_min_steps"],
        "root_nodding_minimum_total_changes": pars["root_nodding_min_steps"],
        "n_nodding_clades_wetter_than_immediate_erect_sister": n_wet,
        "n_nodding_clades_drier_than_immediate_erect_sister": n_dry,
        "comparisons": pairs,
        "repeated_independent_gain_count_established": False,
        "historical_rain_at_origins_estimated": False,
        "rain_dependent_evolutionary_rate_test_performed": False,
        "same_244_taxa_as_global_S1": False,
        "covers_both_S2_species": False,
        "interpretation": (
            "Present-day author species-level BIO12 of maximal pure-nodding tip clades "
            "versus immediate all-erect sister clades. Only 2/4 contrasts wetter. "
            "Rooted maximum-parsimony minimum total changes are conditional, not gain counts. "
            "The existing NA63 BIO12/BIO1/elevation Mk null remains primary."
        ),
    }


def run_original():
    states, trait_sha256 = download_author_states()
    tree_raw, tree_sha256 = fetch_author(
        AUTHOR_NUCLEAR_TREE_PATH, AUTHOR_NUCLEAR_TREE_GIT_BLOB)
    env_raw, env_sha256 = fetch_author(ENV_PATH, ENV_GIT_BLOB)
    tree = read_author_nuclear_tree(tree_raw, states)
    climate = read_published_climate(env_raw)
    output = sister_clade_contrasts(tree, states, climate)
    output["provenance"] = {
        "published_DOI": SOURCE_DOI,
        "github_repo": SOURCE_REPO,
        "github_commit": SOURCE_COMMIT,
        "original_trait_sha256": trait_sha256,
        "original_nuclear_tree_sha256": tree_sha256,
        "original_environment_sha256": env_sha256,
        "source_paths": {
            "traits": "traits/traits.csv",
            "tree": AUTHOR_NUCLEAR_TREE_PATH,
            "environment": ENV_PATH,
        },
        "original_git_blob_ids": {
            "traits": "fa127f0587b676e6e004ea121424d26069586bcb",
            "tree": AUTHOR_NUCLEAR_TREE_GIT_BLOB,
            "environment": ENV_GIT_BLOB,
        }
    }
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    result = run_original()
    dest = Path(args.out)
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(
        json.dumps(result, indent=2, ensure_ascii=False, sort_keys=True) + "\n",
        encoding="utf-8")
    print("T2_NA63_SISTER_CLADE_RAIN=" + json.dumps({
        "root_erect_minimum_total_changes": result["root_erect_minimum_total_changes"],
        "pure_nodding_clades": result["n_maximal_pure_nodding_tip_clades"],
        "wetter": result["n_nodding_clades_wetter_than_immediate_erect_sister"],
        "drier": result["n_nodding_clades_drier_than_immediate_erect_sister"]
    }, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
