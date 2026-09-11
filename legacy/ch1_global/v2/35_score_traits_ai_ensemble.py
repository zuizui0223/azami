#!/usr/bin/env python3
"""Generate fully automated, model-derived Cirsium trait measurements.

This pipeline does not use human annotations. It combines two zero-shot CLIP
backbones, prompt ensembles, and two deterministic image variants (original and
horizontal flip). The resulting states are model-derived measurements with
explicit agreement and stability metadata, not manually validated observations.
"""

from __future__ import annotations

import argparse
import gc
import hashlib
import json
import os
import random
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from PIL import Image, ImageOps


DEFAULT_MODELS = ["openai/clip-vit-base-patch32", "openai/clip-vit-base-patch16"]
IMAGE_VARIANTS = ("original", "hflip")
NON_SCOREABLE_STATE = "unassessable"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Create fully automated ensemble trait measurements from a detected-head packet.")
    parser.add_argument("--packet-manifest", required=True)
    parser.add_argument("--packet-root", required=True)
    parser.add_argument("--ontology", required=True)
    parser.add_argument("--prompt-spec", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--model-id", action="append", default=[])
    parser.add_argument("--device", default="cpu")
    parser.add_argument("--batch-size", type=int, default=24)
    parser.add_argument("--threads", type=int, default=0)
    parser.add_argument("--seed", type=int, default=20260705)
    parser.add_argument("--strict-margin-percentile", type=float, default=0.50)
    parser.add_argument("--majority-vote-fraction", type=float, default=0.75)
    return parser.parse_args()


def text(value: Any) -> str:
    return "" if pd.isna(value) else str(value).strip()


def as_bool(value: Any) -> bool:
    return text(value).lower() in {"true", "1", "yes"}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def json_sha256(value: Any) -> str:
    canonical = json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


def normalize(tensor: Any):
    return tensor / tensor.norm(dim=-1, keepdim=True).clamp_min(1e-12)


def packet_path(value: Any, root: Path) -> Path:
    relative = Path(text(value))
    if relative.is_absolute() or ".." in relative.parts:
        raise ValueError(f"Unsafe packet image path: {relative}")
    return root / relative


def read_variant(path: Path, variant: str) -> Image.Image:
    with Image.open(path) as opened:
        image = opened.convert("RGB")
    if variant == "original":
        return image
    if variant == "hflip":
        return ImageOps.mirror(image)
    raise ValueError(f"Unknown image variant: {variant}")


def load_trait_specs(ontology: pd.DataFrame, prompt_spec: dict[str, Any]) -> list[dict[str, Any]]:
    required = {"trait_id", "allowed_states", "image_annotatable", "layer", "allow_multiple"}
    if missing := required.difference(ontology.columns):
        raise ValueError(f"Ontology missing columns: {sorted(missing)}")
    configs = prompt_spec.get("traits")
    if not isinstance(configs, dict):
        raise ValueError("Prompt specification must include a traits object")
    result: list[dict[str, Any]] = []
    for row in ontology.to_dict("records"):
        trait_id = text(row["trait_id"])
        if not as_bool(row["image_annotatable"]) or text(row["layer"]).lower() == "detector":
            continue
        config = configs.get(trait_id)
        if not isinstance(config, dict):
            raise ValueError(f"Image trait {trait_id!r} is missing a prompt specification")
        allowed_states = [state for state in text(row["allowed_states"]).split("|") if state]
        mode = text(config.get("mode", "score")) or "score"
        if mode == "not_scored":
            result.append({
                "trait_id": trait_id,
                "mode": "not_scored",
                "allowed_states": allowed_states,
                "reason": text(config.get("reason", "Not automatically estimated.")),
            })
            continue
        view = text(config.get("view"))
        if view not in {"head", "context"}:
            raise ValueError(f"Trait {trait_id!r} requires view=head or view=context")
        state_prompts = config.get("state_prompts")
        if not isinstance(state_prompts, dict) or not state_prompts:
            raise ValueError(f"Trait {trait_id!r} has no valid state prompt specification")
        scoreable_states = [state for state in allowed_states if state != NON_SCOREABLE_STATE]
        if set(state_prompts) != set(scoreable_states):
            raise ValueError(
                f"Trait {trait_id!r} prompt states do not match ontology; "
                f"expected={sorted(scoreable_states)} actual={sorted(state_prompts)}"
            )
        for state, prompts in state_prompts.items():
            if not isinstance(prompts, list) or not prompts or not all(isinstance(prompt, str) and prompt.strip() for prompt in prompts):
                raise ValueError(f"Trait {trait_id!r}, state {state!r} needs non-empty prompts")
        result.append({
            "trait_id": trait_id,
            "mode": "score",
            "view": view,
            "state_prompts": state_prompts,
            "allowed_states": allowed_states,
            "allow_multiple": as_bool(row["allow_multiple"]),
        })
    if not result:
        raise ValueError("No image-annotatable traits found")
    return result


def encode_images(
    paths: list[Path],
    processor: Any,
    model: Any,
    device: str,
    batch_size: int,
    variant: str,
) -> np.ndarray:
    import torch

    encoded: list[np.ndarray] = []
    for start in range(0, len(paths), batch_size):
        batch_paths = paths[start:start + batch_size]
        images = [read_variant(path, variant) for path in batch_paths]
        inputs = processor(images=images, return_tensors="pt")
        with torch.inference_mode():
            features = normalize(model.get_image_features(pixel_values=inputs["pixel_values"].to(device)))
        encoded.append(features.cpu().numpy().astype(np.float32))
        print(f"[INFO] {variant}: {min(start + len(batch_paths), len(paths))}/{len(paths)}")
    return np.concatenate(encoded, axis=0)


def encode_state_prompts(state_prompts: dict[str, list[str]], processor: Any, model: Any, device: str) -> tuple[list[str], np.ndarray]:
    import torch

    states = list(state_prompts)
    state_vectors: list[np.ndarray] = []
    for state in states:
        inputs = processor(text=state_prompts[state], padding=True, truncation=True, return_tensors="pt")
        inputs = {key: value.to(device) for key, value in inputs.items()}
        with torch.inference_mode():
            features = normalize(model.get_text_features(**inputs))
        mean_vector = normalize(features.mean(dim=0, keepdim=True)).squeeze(0)
        state_vectors.append(mean_vector.cpu().numpy().astype(np.float32))
    return states, np.stack(state_vectors, axis=0)


def build_variant_rows(
    unit_ids: list[str],
    trait: dict[str, Any],
    image_embeddings: dict[str, np.ndarray],
    processor: Any,
    model: Any,
    device: str,
    model_id: str,
    model_revision: str,
) -> list[dict[str, Any]]:
    states, state_vectors = encode_state_prompts(trait["state_prompts"], processor, model, device)
    rows: list[dict[str, Any]] = []
    for variant, features in image_embeddings.items():
        similarities = features @ state_vectors.T
        order = np.argsort(-similarities, axis=1)
        top = order[:, 0]
        runner = order[:, 1] if len(states) > 1 else top
        for index, unit_id in enumerate(unit_ids):
            rows.append({
                "annotation_unit_id": unit_id,
                "trait_id": trait["trait_id"],
                "image_view": trait["view"],
                "model_id": model_id,
                "model_revision": model_revision,
                "image_variant": variant,
                "variant_candidate_state": states[int(top[index])],
                "variant_runner_up_state": states[int(runner[index])],
                "variant_top_similarity": float(similarities[index, top[index]]),
                "variant_runner_up_similarity": float(similarities[index, runner[index]]),
                "variant_similarity_margin": float(similarities[index, top[index]] - similarities[index, runner[index]]),
            })
    return rows


def aggregate_group(group: pd.DataFrame) -> dict[str, Any]:
    counts = Counter(group["variant_candidate_state"].tolist())
    max_votes = max(counts.values())
    tied = sorted(state for state, count in counts.items() if count == max_votes)
    if len(tied) == 1:
        chosen = tied[0]
    else:
        score = group.groupby("variant_candidate_state")["variant_top_similarity"].mean().to_dict()
        chosen = sorted(tied, key=lambda state: (-float(score[state]), state))[0]
    winners = group.loc[group["variant_candidate_state"].eq(chosen)]
    model_votes = group.groupby("model_id")["variant_candidate_state"].agg(lambda values: "|".join(sorted(set(values)))).to_dict()
    return {
        "annotation_unit_id": text(group.iloc[0]["annotation_unit_id"]),
        "trait_id": text(group.iloc[0]["trait_id"]),
        "ai_candidate_state": chosen,
        "n_variant_votes": int(len(group)),
        "n_votes_for_candidate": int(max_votes),
        "vote_fraction": float(max_votes / len(group)),
        "n_distinct_variant_states": int(len(counts)),
        "variant_vote_counts_json": json.dumps(dict(sorted(counts.items())), ensure_ascii=False, sort_keys=True),
        "model_candidate_states_json": json.dumps(dict(sorted(model_votes.items())), ensure_ascii=False, sort_keys=True),
        "median_similarity_margin": float(group["variant_similarity_margin"].median()),
        "minimum_similarity_margin": float(group["variant_similarity_margin"].min()),
        "mean_top_similarity": float(group["variant_top_similarity"].mean()),
        "mean_winner_top_similarity": float(winners["variant_top_similarity"].mean()),
    }


def status_from_consensus(vote_fraction: float, margin_percentile: float, strict_threshold: float, majority_threshold: float) -> str:
    if vote_fraction == 1.0 and margin_percentile >= strict_threshold:
        return "strict_ensemble_consensus_higher_margin"
    if vote_fraction == 1.0:
        return "strict_ensemble_consensus_lower_margin"
    if vote_fraction >= majority_threshold:
        return "majority_ensemble_candidate"
    return "model_variant_disagreement"


def safe_pivot(long: pd.DataFrame, values: list[str]) -> pd.DataFrame:
    output = long[["annotation_unit_id"]].drop_duplicates().sort_values("annotation_unit_id").reset_index(drop=True)
    for value in values:
        subset = long.pivot(index="annotation_unit_id", columns="trait_id", values=value)
        subset.columns = [f"ai__{trait_id}__{value}" for trait_id in subset.columns]
        output = output.merge(subset.reset_index(), on="annotation_unit_id", how="left", validate="one_to_one")
    return output


def main() -> None:
    args = parse_args()
    if args.batch_size < 1:
        raise ValueError("--batch-size must be >= 1")
    if args.threads < 0:
        raise ValueError("--threads must be >= 0")
    if not 0 <= args.strict_margin_percentile <= 1:
        raise ValueError("--strict-margin-percentile must be in [0, 1]")
    if not 0 < args.majority_vote_fraction <= 1:
        raise ValueError("--majority-vote-fraction must be in (0, 1]")
    model_ids = args.model_id or list(DEFAULT_MODELS)
    if len(model_ids) < 2:
        raise ValueError("At least two --model-id values are required for an ensemble")
    if len(set(model_ids)) != len(model_ids):
        raise ValueError("Model IDs must be unique")

    packet_manifest_path = Path(args.packet_manifest)
    packet_root = Path(args.packet_root)
    ontology_path = Path(args.ontology)
    prompt_path = Path(args.prompt_spec)
    if not packet_manifest_path.is_file() or not ontology_path.is_file() or not prompt_path.is_file() or not packet_root.is_dir():
        raise FileNotFoundError("Packet manifest/root, ontology, or prompt spec is missing")

    packet = pd.read_csv(packet_manifest_path, dtype=str, keep_default_na=False)
    required_packet = {"annotation_unit_id", "source_image", "crop_path", "context_crop_path"}
    if missing := required_packet.difference(packet.columns):
        raise ValueError(f"Packet manifest missing columns: {sorted(missing)}")
    if packet.empty or packet["annotation_unit_id"].duplicated().any():
        raise ValueError("Packet manifest must have nonempty unique annotation_unit_id values")
    ontology = pd.read_csv(ontology_path, dtype=str, keep_default_na=False)
    prompt_spec = json.loads(prompt_path.read_text(encoding="utf-8"))
    traits = load_trait_specs(ontology, prompt_spec)
    scored_traits = [trait for trait in traits if trait["mode"] == "score"]
    not_scored_traits = [trait for trait in traits if trait["mode"] == "not_scored"]

    try:
        import torch
        import transformers
        from transformers import CLIPModel, CLIPProcessor
    except ImportError as error:
        raise SystemExit("Install torch and transformers before running the AI ensemble.") from error

    random.seed(args.seed)
    np.random.seed(args.seed)
    torch.manual_seed(args.seed)
    if args.threads:
        torch.set_num_threads(args.threads)
    if args.device != "cpu" and not torch.cuda.is_available():
        raise ValueError(f"Requested device={args.device!r} but CUDA is unavailable")
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False

    view_paths: dict[str, list[Path]] = {}
    for view, field in (("head", "crop_path"), ("context", "context_crop_path")):
        if any(trait.get("view") == view for trait in scored_traits):
            paths = [packet_path(value, packet_root) for value in packet[field]]
            missing_images = [str(path) for path in paths if not path.is_file()]
            if missing_images:
                raise FileNotFoundError(f"Missing {view} images: {missing_images[:5]}")
            view_paths[view] = paths

    unit_ids = packet["annotation_unit_id"].tolist()
    variant_rows: list[dict[str, Any]] = []
    model_locks: list[dict[str, str]] = []
    for model_id in model_ids:
        print(f"[INFO] loading model {model_id}")
        processor = CLIPProcessor.from_pretrained(model_id)
        model = CLIPModel.from_pretrained(model_id).to(args.device)
        model.eval()
        model_revision = text(getattr(getattr(model, "config", None), "_commit_hash", "")) or "unresolved_by_transformers"
        model_locks.append({"model_id": model_id, "resolved_revision": model_revision})
        embeddings: dict[str, dict[str, np.ndarray]] = {}
        for view, paths in view_paths.items():
            embeddings[view] = {
                variant: encode_images(paths, processor, model, args.device, args.batch_size, variant)
                for variant in IMAGE_VARIANTS
            }
        for trait in scored_traits:
            variant_rows.extend(build_variant_rows(
                unit_ids=unit_ids,
                trait=trait,
                image_embeddings=embeddings[trait["view"]],
                processor=processor,
                model=model,
                device=args.device,
                model_id=model_id,
                model_revision=model_revision,
            ))
        del embeddings
        del model
        del processor
        gc.collect()
        if args.device != "cpu":
            torch.cuda.empty_cache()

    variants = pd.DataFrame(variant_rows)
    expected_variant_rows = len(packet) * len(scored_traits) * len(model_ids) * len(IMAGE_VARIANTS)
    if len(variants) != expected_variant_rows:
        raise ValueError(f"Unexpected variant-prediction row count: got={len(variants)} expected={expected_variant_rows}")
    aggregated = pd.DataFrame([
        aggregate_group(group)
        for _, group in variants.groupby(["annotation_unit_id", "trait_id"], sort=True)
    ])
    expected_scored_rows = len(packet) * len(scored_traits)
    if len(aggregated) != expected_scored_rows:
        raise ValueError(f"Unexpected ensemble measurement row count: got={len(aggregated)} expected={expected_scored_rows}")
    aggregated["median_margin_percentile_within_trait"] = aggregated.groupby("trait_id")["median_similarity_margin"].rank(method="average", pct=True)
    aggregated["ai_measurement_status"] = [
        status_from_consensus(vote_fraction, percentile, args.strict_margin_percentile, args.majority_vote_fraction)
        for vote_fraction, percentile in zip(aggregated["vote_fraction"], aggregated["median_margin_percentile_within_trait"])
    ]
    aggregated["analysis_state_ai_all"] = aggregated["ai_candidate_state"]
    aggregated["analysis_state_ai_conservative"] = aggregated["ai_candidate_state"].where(
        aggregated["ai_measurement_status"].eq("strict_ensemble_consensus_higher_margin"),
        "",
    )
    aggregated["measurement_basis"] = "fully_automated_zero_shot_clip_ensemble"
    aggregated["model_count"] = len(model_ids)
    aggregated["image_variant_count"] = len(IMAGE_VARIANTS)

    for trait in not_scored_traits:
        additional = pd.DataFrame({
            "annotation_unit_id": unit_ids,
            "trait_id": trait["trait_id"],
            "ai_candidate_state": "",
            "n_variant_votes": 0,
            "n_votes_for_candidate": 0,
            "vote_fraction": np.nan,
            "n_distinct_variant_states": 0,
            "variant_vote_counts_json": "{}",
            "model_candidate_states_json": "{}",
            "median_similarity_margin": np.nan,
            "minimum_similarity_margin": np.nan,
            "mean_top_similarity": np.nan,
            "mean_winner_top_similarity": np.nan,
            "median_margin_percentile_within_trait": np.nan,
            "ai_measurement_status": "not_auto_estimated_targeted_macro_required",
            "analysis_state_ai_all": "",
            "analysis_state_ai_conservative": "",
            "measurement_basis": "automated_exclusion_generic_observer_photos_insufficient",
            "model_count": len(model_ids),
            "image_variant_count": len(IMAGE_VARIANTS),
            "not_auto_estimated_reason": trait["reason"],
        })
        aggregated = pd.concat([aggregated, additional], ignore_index=True, sort=False)

    aggregated["not_auto_estimated_reason"] = aggregated.get("not_auto_estimated_reason", "").fillna("")
    aggregated = aggregated.sort_values(["annotation_unit_id", "trait_id"]).reset_index(drop=True)
    long_columns = [
        "annotation_unit_id", "trait_id", "analysis_state_ai_all", "analysis_state_ai_conservative", "ai_candidate_state",
        "ai_measurement_status", "measurement_basis", "vote_fraction", "n_variant_votes", "n_votes_for_candidate",
        "n_distinct_variant_states", "median_similarity_margin", "minimum_similarity_margin",
        "median_margin_percentile_within_trait", "mean_top_similarity", "mean_winner_top_similarity",
        "variant_vote_counts_json", "model_candidate_states_json", "model_count", "image_variant_count", "not_auto_estimated_reason",
    ]
    measurements = aggregated[long_columns].copy()
    wide = safe_pivot(measurements, [
        "analysis_state_ai_all", "analysis_state_ai_conservative", "ai_measurement_status", "vote_fraction",
        "median_margin_percentile_within_trait", "measurement_basis",
    ])
    status_summary = (
        measurements.groupby(["trait_id", "ai_measurement_status"], dropna=False)
        .size()
        .rename("n_measurements")
        .reset_index()
        .sort_values(["trait_id", "ai_measurement_status"])
    )

    output = Path(args.out_dir)
    output.mkdir(parents=True, exist_ok=True)
    variants.to_csv(output / "ai_ensemble_variant_predictions.csv", index=False, encoding="utf-8-sig")
    measurements.to_csv(output / "ai_ensemble_trait_measurements_long.csv", index=False, encoding="utf-8-sig")
    wide.to_csv(output / "ai_ensemble_trait_measurements_wide.csv", index=False, encoding="utf-8-sig")
    status_summary.to_csv(output / "ai_ensemble_status_summary.csv", index=False, encoding="utf-8-sig")
    report = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "semantic_status": "Fully automated model-derived trait measurements. No human annotations were used. These are not manually validated observations or calibrated biological ground truth.",
        "n_annotation_units": int(len(packet)),
        "n_image_traits_scored": int(len(scored_traits)),
        "n_image_traits_not_auto_estimated": int(len(not_scored_traits)),
        "n_measurement_rows": int(len(measurements)),
        "n_variant_prediction_rows": int(len(variants)),
        "model_locks": model_locks,
        "image_variants": list(IMAGE_VARIANTS),
        "seed": args.seed,
        "device": args.device,
        "batch_size": args.batch_size,
        "torch_version": torch.__version__,
        "transformers_version": transformers.__version__,
        "packet_manifest_sha256": sha256_file(packet_manifest_path),
        "ontology_sha256": sha256_file(ontology_path),
        "prompt_spec_sha256": json_sha256(prompt_spec),
        "strict_margin_percentile": args.strict_margin_percentile,
        "majority_vote_fraction": args.majority_vote_fraction,
        "status_counts_by_trait": {
            trait_id: group.set_index("ai_measurement_status")["n_measurements"].astype(int).to_dict()
            for trait_id, group in status_summary.groupby("trait_id", sort=True)
        },
        "interpretation": {
            "analysis_state_ai_all": "Always the ensemble winner for scored traits; carry ai_measurement_status into any downstream analysis.",
            "analysis_state_ai_conservative": "Only present for strict agreement across all model/augmentation variants with above-median within-trait margin.",
            "external_mucilage_visible": "Automatically retained as not_auto_estimated_targeted_macro_required because generic observer photographs are not an adequate measurement channel.",
        },
    }
    (output / "ai_ensemble_measurement_provenance.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({
        "n_annotation_units": report["n_annotation_units"],
        "n_measurement_rows": report["n_measurement_rows"],
        "n_variant_prediction_rows": report["n_variant_prediction_rows"],
        "model_locks": model_locks,
        "prompt_spec_sha256": report["prompt_spec_sha256"],
    }, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
