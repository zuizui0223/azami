from __future__ import annotations

import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
REGISTRY_PATH = Path(__file__).with_name("workflow_surface_registry.json")
WORKFLOW_DIR = REPO_ROOT / ".github" / "workflows"


def main() -> int:
    registry = json.loads(REGISTRY_PATH.read_text(encoding="utf-8"))
    categories = registry["categories"]

    registered: list[str] = []
    for category, names in categories.items():
        if not isinstance(names, list):
            raise SystemExit(f"category {category!r} is not a list")
        registered.extend(names)

    duplicates = sorted({name for name in registered if registered.count(name) > 1})
    if duplicates:
        raise SystemExit(f"workflow registry contains duplicates: {duplicates}")

    existing = sorted(
        p.name
        for p in WORKFLOW_DIR.iterdir()
        if p.is_file()
        and p.name.startswith("ch1-v3-")
        and p.suffix in {".yml", ".yaml"}
    )
    registered_set = set(registered)
    existing_set = set(existing)

    unregistered = sorted(existing_set - registered_set)
    missing = sorted(registered_set - existing_set)

    if unregistered or missing:
        if unregistered:
            print("UNREGISTERED ch1-v3 workflows:")
            for name in unregistered:
                print(f"  - {name}")
        if missing:
            print("REGISTERED but missing workflows:")
            for name in missing:
                print(f"  - {name}")
        return 1

    policy = registry.get("policy", {})
    if policy.get("new_scientific_analysis_allowed") is not False:
        raise SystemExit("workflow-only phase must keep new_scientific_analysis_allowed=false")
    if policy.get("unregistered_ch1_v3_workflow_allowed") is not False:
        raise SystemExit("workflow-only phase must reject unregistered ch1-v3 workflows")

    print(f"PASS: {len(existing)} Chapter 1 v3 workflows are registered")
    for category, names in categories.items():
        print(f"  {category}: {len(names)}")
    print("Scientific analysis scope remains frozen; this guard checks workflow surface only.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
