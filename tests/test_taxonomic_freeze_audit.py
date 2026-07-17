import importlib.util
from pathlib import Path


MODULE_PATH = Path(__file__).parents[1] / "analysis" / "audit_taxonomic_freeze.py"
spec = importlib.util.spec_from_file_location("taxonomic_freeze", MODULE_PATH)
module = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(module)


def decision(input_name: str, accepted_name: str, status: str = "accepted") -> dict[str, str]:
    return {
        "input_name": input_name,
        "accepted_name": accepted_name,
        "decision": status,
        "authority_source": "POWO",
        "authority_record": "record",
        "checked_date": "2026-07-17",
        "notes": "",
    }


def test_valid_decisions_are_indexed_once():
    indexed, errors = module.validate_decisions(
        [decision("Cirsium a", "Cirsium a"), decision("Cnicus b", "Cirsium b", "synonym")]
    )
    assert errors == []
    assert indexed["Cnicus b"]["accepted_name"] == "Cirsium b"


def test_duplicate_input_name_fails():
    _, errors = module.validate_decisions(
        [decision("Cirsium a", "Cirsium a"), decision("Cirsium a", "Cirsium b")]
    )
    assert any("exactly one required" in error for error in errors)


def test_provisional_requires_an_accepted_name_and_is_allowed_for_review_stage():
    indexed, errors = module.validate_decisions(
        [decision("Cirsium uncertain", "Cirsium uncertain", "provisional")]
    )
    assert errors == []
    assert indexed["Cirsium uncertain"]["decision"] == "provisional"


def test_missing_authority_is_rejected():
    row = decision("Cirsium a", "Cirsium a")
    row["authority_source"] = ""
    _, errors = module.validate_decisions([row])
    assert any("authority_source is blank" in error for error in errors)
