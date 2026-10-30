"""Read-only inventory of GitHub branch and PR references.

This script intentionally has NO branch deletion, ref update, PR closure, merge,
or release publication operation. A candidate means 'needs manual review'.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
from collections import defaultdict
from pathlib import Path
from urllib.request import Request, urlopen

REPO_PATTERN = re.compile(r"^[A-Za-z0-9_.-]+/[A-Za-z0-9_.-]+$")


def get_all(repo: str, endpoint: str, token: str = "") -> list[dict]:
    """GET all pages, stopping only at an actual short page."""
    rows: list[dict] = []
    for page in range(1, 51):
        separator = "&" if "?" in endpoint else "?"
        url = f"https://api.github.com/repos/{repo}/{endpoint}{separator}per_page=100&page={page}"
        headers = {
            "Accept": "application/vnd.github+json",
            "User-Agent": "azami-read-only-ref-audit",
            "X-GitHub-Api-Version": "2022-11-28",
        }
        if token:
            headers["Authorization"] = f"Bearer {token}"
        req = Request(url, headers=headers, method="GET")
        with urlopen(req, timeout=30) as response:
            payload = json.load(response)
        if not isinstance(payload, list):
            raise ValueError(f"Expected a list from {endpoint} page {page}")
        rows.extend(payload)
        if len(payload) < 100:
            return rows
    raise RuntimeError(f"Pagination limit hit for {endpoint}; incomplete inventory")


def build_report(
    branches: list[dict],
    pulls: list[dict],
    tags: list[dict],
    releases: list[dict],
    default_branch: str,
) -> dict:
    """Classify references without asserting that any can safely be deleted."""
    name_sha: dict[str, str] = {}
    grouped: dict[str, list[str]] = defaultdict(list)
    for branch in branches:
        name = branch["name"]
        sha = branch["commit"]["sha"]
        if name in name_sha:
            raise ValueError(f"Duplicate branch in paginated results: {name}")
        name_sha[name] = sha
        grouped[sha].append(name)

    pr_heads = {p["head"]["ref"] for p in pulls}
    open_pulls = [p for p in pulls if p["state"] == "open"]
    active_pr_refs = {
        ref
        for p in open_pulls
        for ref in (p["head"]["ref"], p["base"]["ref"])
    }
    protected = active_pr_refs | {default_branch}
    groups = []
    for sha, names in sorted(grouped.items(), key=lambda x: (-len(x[1]), x[0])):
        if len(names) < 2:
            continue
        group_names = sorted(names)
        no_pr_head = sorted(set(names) - pr_heads)
        groups.append({
            "sha": sha,
            "count": len(names),
            "branch_names": group_names,
            "historical_pr_head_match_count": len(names) - len(no_pr_head),
            "no_matching_pr_head_names": no_pr_head,
            "open_pr_or_default_refs": sorted(set(names) & protected),
            "assessment": "REVIEW_ONLY_NEVER_AUTODELETE",
        })

    pr_graph = [
        {
            "number": p["number"],
            "title": p.get("title", ""),
            "head": p["head"]["ref"],
            "base": p["base"]["ref"],
            "head_sha": p["head"]["sha"],
            "head_branch_exists": p["head"]["ref"] in name_sha,
            "base_branch_exists": p["base"]["ref"] in name_sha,
        }
        for p in open_pulls
    ]
    pr_graph.sort(key=lambda p: p["number"], reverse=True)
    tag_rows = [
        {"ref": t["ref"], "target_sha": t["object"]["sha"], "target_type": t["object"].get("type", "unknown")}
        for t in tags
    ]
    tag_rows.sort(key=lambda x: x["ref"])
    return {
        "schema_version": 1,
        "scope": "public GitHub refs; read-only; no deletion recommended",
        "branch_count": len(name_sha),
        "pull_request_count": len(pulls),
        "open_pull_request_count": len(open_pulls),
        "default_branch": default_branch,
        "shared_tip_groups": groups,
        "open_pr_dependencies": pr_graph,
        "tag_refs": tag_rows,
        "release_count": len(releases),
        "draft_release_count": sum(bool(r.get("draft")) for r in releases),
        "notes": [
            "A matching historical PR head name does not mean the original PR head SHA equals the current branch tip.",
            "A branch lacking a PR head-name match can still be referenced by scripts, external links, or deployments.",
            "Annotated tag objects may need dereferencing before commit-level equivalence claims.",
            "Never delete a branch, close a PR, or change a ref based solely on this report.",
        ],
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", default="zuizui0223/azami")
    parser.add_argument("--default-branch", default="main")
    parser.add_argument("--input-json", type=Path, help="Offline fixture with branches, pulls, tags, releases")
    parser.add_argument("--output-json", type=Path, help="Local output only; defaults to stdout")
    args = parser.parse_args(argv)
    if not REPO_PATTERN.fullmatch(args.repo):
        parser.error("Invalid owner/repository syntax")

    if args.input_json:
        inputs = json.loads(args.input_json.read_text(encoding="utf-8"))
    else:
        token = os.environ.get("GITHUB_TOKEN", "")
        inputs = {
            "branches": get_all(args.repo, "branches", token),
            "pulls": get_all(args.repo, "pulls?state=all", token),
            "tags": get_all(args.repo, "git/refs/tags", token),
            "releases": get_all(args.repo, "releases", token),
        }
    report = build_report(
        inputs["branches"], inputs["pulls"], inputs.get("tags", []),
        inputs.get("releases", []), args.default_branch,
    )
    payload = json.dumps(report, ensure_ascii=False, indent=2) + "\n"
    if args.output_json:
        args.output_json.parent.mkdir(parents=True, exist_ok=True)
        args.output_json.write_text(payload, encoding="utf-8")
    else:
        sys.stdout.write(payload)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
