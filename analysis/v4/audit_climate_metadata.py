"""Retrieve checksum-locked cache metadata; never expose observation records."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import zipfile

from analysis.v4.native_665115_pipeline import CLIMATE

ARTIFACTS = dict(zip(CLIMATE, (
    12243016948, 12241949860, 12242703088, 12242508971,
    12242917323, 12242842598, 12243355603, 12243306531, 12242657215)))


def execute(output):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    result = {}
    for variable, artifact in ARTIFACTS.items():
        column, expected = CLIMATE[variable]
        archive = output.parent / ("metadata-audit-" + variable + ".zip")
        with archive.open("wb") as handle:
            subprocess.run(["gh", "api",
                            f"repos/zuizui0223/azami/actions/artifacts/{artifact}/zip"],
                           stdout=handle, check=True)
        with archive.open("rb") as handle:
            digest = hashlib.file_digest(handle, "sha256").hexdigest()
        if digest != expected:
            raise ValueError("Cache checksum mismatch: " + variable)
        with zipfile.ZipFile(archive) as z:
            names = [n for n in z.namelist() if n.endswith("grid.json")]
            if len(names) != 1:
                raise ValueError("Ambiguous metadata: " + variable)
            metadata = json.loads(z.read(names[0]))
        result[column] = {k:v for k,v in metadata.items()
                          if k not in ("tiles", "cache_sha256")}
        result[column]["cache_artifact_id"] = artifact
        result[column]["cache_archive_sha256"] = digest
        archive.unlink()
        print(column, result[column], flush=True)
    (output / "all_environment_source_metadata.json").write_text(
        json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    execute(sys.argv[1])
