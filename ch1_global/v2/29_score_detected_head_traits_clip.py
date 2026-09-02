#!/usr/bin/env python3
"""Generate zero-shot CLIP trait proposals for a blinded detected-head packet.

Outputs are provisional model suggestions only. They are intentionally written to
files separate from annotation responses: model candidates never overwrite
human labels, and they must not be interpreted as validated trait observations.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from PIL import Image

NONSCORABLE_ALLOWED = {"unassessable"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Score detected-head trait candidates with zero-shot CLIP.")
    parser.add_argument("--packet-manifest", required=True)
    parser.add_argument("--packet-root", required=True)
    parser.add_argument("--ontology", required=True)
    parser.add_argument("--prompt-spec", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--model-id", default="openai/clip-vit-base-patch32")
    parser.add_argument("--device", default="cpu")
    parser.add_argument("--batch-size", type=int, default=16)
    parser.add_argument("--n-priority-units", type=int, default=250)
    return parser.parse_args()


def text(value: Any) -> str:
    return "" if pd.isna(value) else str(value).strip()


def json_hash(payload: Any) -> str:
    canonical = json.dumps(payload, ensure_ascii=False, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


def path_from_packet(value: Any, packet_root: Path) -> Path:
    path = Path(text(value))
    return path if path.is_absolute() else packet_root / path


def read_image(path: Path) -> Image.Image:
    with Image.open(path) as image:
        return image.convert("RGB")


def normalize(array: Any):
    return array / array.norm(dim=-1, keepdim=True).clamp_min(1e-12)


def trait_spec_or_raise(ontology: pd.DataFrame, prompt_spec: dict[str, Any]) -> list[dict[str, Any]]:
    required = {"trait_id", "allowed_states", "image_annotatable", "layer"}
    missing = required.difference(ontology.columns)
    if missing:
        raise ValueError(f"Ontology missing columns: {sorted(missing)}")
    traits_config = prompt_spec.get("traits")
    if not isinstance(traits_config, dict):
        raise ValueError("Prompt spec must contain an object named traits")
    selected: list[dict[str, Any]] = []
    for row in ontology.to_dict("records"):
        trait_id = text(row.get("trait_id"))
        image_annotatable = text(row.get("image_annotatable")).lower() in {"true", "1", "yes"}
        if not image_annotatable or text(row.get("layer")).lower() == "detector":
            continue
        config = traits_config.get(trait_id)
        if config is None:
            raise ValueError(f"No prompt specification for image-annotatable trait {trait_id!r}")
        states = [state for state in text(row.get("allowed_states")).split("|") if state]
        mode = text(config.get("mode", "score")) or "score"
        if mode == "not_scored":
            selected.append({"trait_id": trait_id, "mode": mode, "states": states, "config": config})
            continue
        view = text(config.get("view"))
        if view not in {"head", "context"}:
            raise ValueError(f"Trait {trait_id!r} must specify view=head or context")
        state_prompts = config.get("state_prompts")
        if not isinstance(state_prompts, dict) or not state_prompts:
            raise ValueError(f"Trait {trait_id!r} requires a non-empty state_prompts object")
        scoreable_states = [state for state in states if state not in NONSCORABLE_ALLOWED]
        missing_states = set(scoreable_states).difference(state_prompts)
        extra_states = set(state_prompts).difference(scoreable_states)
        if missing_states or extra_states:
            raise ValueError(f"Trait {trait_id!r} prompt-state mismatch; missing={sorted(missing_states)}, extra={sorted(extra_states)}")
        for state, prompts in state_prompts.items():
            if not isinstance(prompts, list) or not prompts or not all(isinstance(prompt, str) and prompt.strip() for prompt in prompts):
                raise ValueError(f"Trait {trait_id!r}, state {state!r} needs a non-empty prompt list")
        selected.append({
            "trait_id": trait_id,
            "mode": mode,
            "states": scoreable_states,
            "config": config,
            "view": view,
            "allow_multiple": text(row.get("allow_multiple")).lower() in {"true", "1", "yes"},
        })
    if not selected:
        raise ValueError("No eligible ontology traits")
    return selected


def encode_images(paths: list[Path], processor: Any, model: Any, device: str, batch_size: int) -> np.ndarray:
    import torch

    vectors: list[np.ndarray] = []
    for start in range(0, len(paths), batch_size):
        images = [read_image(path) for path in paths[start:start + batch_size]]
        inputs = processor(images=images, return_tensors="pt")
        with torch.inference_mode():
            features = normalize(model.get_image_features(pixel_values=inputs["pixel_values"].to(device)))
        vectors.append(features.cpu().numpy().astype(np.float32))
        print(f"[INFO] encoded images {min(start + len(images), len(paths))}/{len(paths)}")
    return np.concatenate(vectors, axis=0)


def encode_state_text(state_prompts: dict[str, list[str]], processor: Any, model: Any, device: str) -> tuple[list[str], np.ndarray]:
    import torch

    states = list(state_prompts)
    state_vectors: list[np.ndarray] = []
    for state in states:
        inputs = processor(text=state_prompts[state], return_tensors="pt", padding=True, truncation=True)
        inputs = {key: value.to(device) for key, value in inputs.items()}
        with torch.inference_mode():
            features = normalize(model.get_text_features(**inputs))
        averaged = normalize(features.mean(dim=0, keepdim=True)).squeeze(0)
        state_vectors.append(averaged.cpu().numpy().astype(np.float32))
    return states, np.stack(state_vectors, axis=0)


def margin_percentile(values: np.ndarray) -> np.ndarray:
    return pd.Series(values).rank(method="average", pct=True).to_numpy(dtype=float)


def confidence_tier(percentile: float) -> str:
    if percentile <= 1 / 3:
        return "low_margin"
    if percentile <= 2 / 3:
        return "mid_margin"
    return "high_margin"


def priority_weight(config: dict[str, Any]) -> float:
    try:
        value = float(config.get("review_weight", 1.0))
    except (TypeError, ValueError) as error:
        raise ValueError(f"Invalid review_weight: {config.get('review_weight')!r}") from error
    if value <= 0:
        raise ValueError("review_weight must be positive")
    return value


def main() -> None:
    args = parse_args()
    if args.batch_size < 1 or args.n_priority_units < 1:
        raise ValueError("--batch-size and --n-priority-units must be >= 1")
    packet_root = Path(args.packet_root)
    manifest_path = Path(args.packet_manifest)
    if not packet_root.is_dir() or not manifest_path.is_file():
        raise FileNotFoundError("Packet root or manifest is missing")
    packet = pd.read_csv(manifest_path, dtype=str, keep_default_na=False)
    required_packet = {"annotation_unit_id", "source_image", "crop_path", "context_crop_path"}
    if missing := required_packet.difference(packet.columns):
        raise ValueError(f"Packet manifest missing columns: {sorted(missing)}")
    if packet.empty or packet["annotation_unit_id"].duplicated().any():
        raise ValueError("Packet must have unique, nonempty annotation_unit_id values")
    ontology = pd.read_csv(args.ontology, dtype=str, keep_default_na=False)
    prompt_spec = json.loads(Path(args.prompt_spec).read_text(encoding="utf-8"))
    traits = trait_spec_or_raise(ontology, prompt_spec)

    try:
        import torch
        import transformers
        from transformers import CLIPModel, CLIPProcessor
    except ImportError as error:
        raise SystemExit("Install torch and transformers before CLIP trait scoring.") from error
    if args.device != "cpu" and not torch.cuda.is_available():
        raise ValueError(f"Requested device={args.device!r} but CUDA is unavailable")
    processor = CLIPProcessor.from_pretrained(args.model_id)
    model = CLIPModel.from_pretrained(args.model_id).to(args.device)
    model.eval()

    view_cache: dict[str, np.ndarray] = {}
    for view, field in (("head", "crop_path"), ("context", "context_crop_path")):
        if any(trait.get("view") == view for trait in traits if trait["mode"] != "not_scored"):
            paths = [path_from_packet(value, packet_root) for value in packet[field]]
            missing_paths = [str(path) for path in paths if not path.is_file()]
            if missing_paths:
                raise FileNotFoundError(f"Packet has missing {view} images: {missing_paths[:5]}")
            view_cache[view] = encode_images(paths, processor, model, args.device, args.batch_size)

    proposal_rows: list[dict[str, Any]] = []
    trait_summaries: dict[str, Any] = {}
    for trait in traits:
        trait_id, config = trait["trait_id"], trait["config"]
        if trait["mode"] == "not_scored":
            for unit_id in packet["annotation_unit_id"]:
                proposal_rows.append({
                    "annotation_unit_id": unit_id, "trait_id": trait_id,
                    "proposal_status": "not_scored_targeted_macro_required",
                    "ai_candidate_state": "", "runner_up_state": "", "top_similarity": "", "runner_up_similarity": "",
                    "similarity_margin": "", "margin_percentile_within_trait": "", "confidence_tier": "", "image_view": "",
                    "model_id": args.model_id, "prompt_spec_version": text(prompt_spec.get("prompt_spec_version")),
                    "allow_multiple_trait_state": "false", "review_priority": "targeted_macro_required",
                    "reason": text(config.get("reason", "This trait is not auto-scored from generic observer photographs.")),
                })
            trait_summaries[trait_id] = {"mode": "not_scored", "n_rows": int(len(packet))}
            continue

        states, state_features = encode_state_text(config["state_prompts"], processor, model, args.device)
        similarities = view_cache[trait["view"]] @ state_features.T
        rank = np.argsort(-similarities, axis=1)
        top_index = rank[:, 0]
        second_index = rank[:, 1] if len(states) > 1 else rank[:, 0]
        top_scores = similarities[np.arange(len(packet)), top_index]
        second_scores = similarities[np.arange(len(packet)), second_index]
        margins = top_scores - second_scores
        percentiles = margin_percentile(margins)
        weight = priority_weight(config)
        base_priority = text(config.get("base_review_priority", "medium"))
        for index, unit_id in enumerate(packet["annotation_unit_id"]):
            percentile = float(percentiles[index])
            tier = confidence_tier(percentile)
            review_priority = "high" if base_priority == "high" or tier == "low_margin" else "medium" if base_priority == "medium" or tier == "mid_margin" else "low"
            proposal_rows.append({
                "annotation_unit_id": unit_id, "trait_id": trait_id,
                "proposal_status": "zero_shot_candidate_not_validated",
                "ai_candidate_state": states[int(top_index[index])], "runner_up_state": states[int(second_index[index])],
                "top_similarity": float(top_scores[index]), "runner_up_similarity": float(second_scores[index]),
                "similarity_margin": float(margins[index]), "margin_percentile_within_trait": percentile,
                "confidence_tier": tier, "image_view": trait["view"], "model_id": args.model_id,
                "prompt_spec_version": text(prompt_spec.get("prompt_spec_version")),
                "allow_multiple_trait_state": str(bool(trait.get("allow_multiple", False))).lower(),
                "review_priority": review_priority,
                "reason": "Zero-shot CLIP state comparison. This is a model proposal, not an annotation or validated trait observation.",
                "_priority_score": (1.0 - percentile) * weight,
            })
        trait_summaries[trait_id] = {
            "mode": "zero_shot_clip", "image_view": trait["view"], "states": states, "n_rows": int(len(packet)),
            "similarity_margin_quantiles": {str(q): float(np.quantile(margins, q)) for q in (0, .1, .25, .5, .75, .9, 1)},
            "review_weight": weight,
        }

    proposals = pd.DataFrame(proposal_rows)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    public_columns = [column for column in proposals.columns if column != "_priority_score"]
    proposals[public_columns].to_csv(out_dir / "zero_shot_trait_proposals_long.csv", index=False, encoding="utf-8-sig")
    scorable = proposals.loc[proposals["proposal_status"].eq("zero_shot_candidate_not_validated")].copy()
    if scorable.empty:
        raise ValueError("No traits were scored")
    unit_scores = scorable.groupby("annotation_unit_id", as_index=False).agg(
        highest_uncertainty_score=("_priority_score", "max"),
        n_high_priority_traits=("review_priority", lambda values: int(sum(value == "high" for value in values))),
        n_low_margin_traits=("confidence_tier", lambda values: int(sum(value == "low_margin" for value in values))),
    )
    unit_scores["stable_tiebreak"] = unit_scores["annotation_unit_id"].map(lambda value: hashlib.sha256(f"20260704:{value}".encode("utf-8")).hexdigest())
    unit_scores = unit_scores.sort_values(["highest_uncertainty_score", "n_high_priority_traits", "n_low_margin_traits", "stable_tiebreak"], ascending=[False, False, False, True])
    unit_scores["triage_rank"] = range(1, len(unit_scores) + 1)
    unit_scores.head(min(args.n_priority_units, len(unit_scores))).to_csv(out_dir / "zero_shot_trait_priority_review_subset.csv", index=False, encoding="utf-8-sig")

    wide = packet[["annotation_unit_id", "annotation_batch", "source_image", "crop_path", "context_crop_path"]].copy()
    for trait_id, part in scorable.groupby("trait_id", sort=True):
        fields = part[["annotation_unit_id", "ai_candidate_state", "runner_up_state", "similarity_margin", "margin_percentile_within_trait", "confidence_tier", "review_priority"]].rename(columns={
            "ai_candidate_state": f"ai__{trait_id}__candidate_state",
            "runner_up_state": f"ai__{trait_id}__runner_up_state",
            "similarity_margin": f"ai__{trait_id}__similarity_margin",
            "margin_percentile_within_trait": f"ai__{trait_id}__margin_percentile",
            "confidence_tier": f"ai__{trait_id}__confidence_tier",
            "review_priority": f"ai__{trait_id}__review_priority",
        })
        wide = wide.merge(fields, on="annotation_unit_id", how="left", validate="one_to_one")
    wide = wide.merge(unit_scores.drop(columns=["stable_tiebreak"]), on="annotation_unit_id", how="left", validate="one_to_one")
    wide.to_csv(out_dir / "zero_shot_trait_proposals_wide.csv", index=False, encoding="utf-8-sig")

    report = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "model_id": args.model_id,
        "model_commit_hash_if_available": text(getattr(getattr(model, "config", None), "_commit_hash", "")),
        "transformers_version": transformers.__version__, "torch_version": torch.__version__,
        "device": args.device, "batch_size": args.batch_size, "n_annotation_units": int(len(packet)),
        "n_trait_rows_total": int(len(proposals)), "n_trait_rows_scored": int(len(scorable)),
        "n_priority_review_units": int(min(args.n_priority_units, len(unit_scores))),
        "prompt_spec_version": text(prompt_spec.get("prompt_spec_version")), "prompt_spec_sha256": json_hash(prompt_spec),
        "trait_summaries": trait_summaries,
        "semantic_status": "zero-shot model candidates only. They are not annotations, not detector validation, and not inputs for biological inference until independently checked.",
        "known_limits": [
            "CLIP similarities are relative comparisons within the supplied state prompts, not calibrated probabilities.",
            "Morphological microtraits and external mucilage are especially unsuitable for automatic acceptance from generic observer photographs.",
            "Traits with allow_multiple=true are represented by a primary candidate plus runner-up, not a final multi-state observation."
        ]
    }
    (out_dir / "zero_shot_trait_scoring_provenance.json").write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({key: report[key] for key in ("model_id", "n_annotation_units", "n_trait_rows_scored", "n_priority_review_units", "prompt_spec_sha256")}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
