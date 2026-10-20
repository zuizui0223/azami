"""Display increments for decoded CHELSA v2.1 layers; never change model fitting."""
import csv

# Units: original raster affine decoding plus CHELSA BIOCLIM+ documentation.
# https://essd.copernicus.org/articles/14/5573/2022/
# BIO4 remains explicitly an index; do not label its decoded values as degrees C.
UNITS = {
    "chelsa_bio01": ("deg C", 1.),
    "chelsa_bio04": ("BIO4 index", 100.),
    "chelsa_bio12": ("mm/year", 100.),
    "chelsa_bio15": ("CV percentage points", 10.),
    "chelsa_rsds_mean": ("MJ/m2/day", 1.),
    "chelsa_vpd_mean": ("Pa", 100.),
    "chelsa_sfcwind_mean": ("m/s", 1.),
    "chelsa_gsp": ("mm/growing season", 100.),
    "chelsa_npp": ("g dry matter/m2/year", 100.),
}


def interval_per_increment(effect, sd, increment):
    if sd <= 0 or increment <= 0:
        raise ValueError("Positive scales required")
    factor = increment / sd
    return (effect["raw_units_per_environment_sd"] * factor,
            [v * factor for v in effect["raw_ci95"]])


def write_report(result, output):
    if result["input_receipt"].get("environment_decoding_status") != \
            "physical_units_gsp_invalid_encoding_masked_v1":
        raise ValueError("Physical-unit reporting requires verified decoded inputs")
    rows = []
    lines = ["# 環境の実単位あたりの形質変化", "",
             "他の環境条件を固定した線形の変化量。形質は元の測定単位。"
             "区間は点別95%で、モジュールの共同検定とは区別する。",
             "種間と種内には同じ環境増分を使用。BIO4は指数単位として明示する。", "",
             "|位置づけ|形質成分|尺度|環境|増分|変化|95%区間|",
             "|---|---|---|---|---|---:|---|"]
    for module, model in result["models"].items():
        for endpoint, component in model.get("components", {}).items():
            for term, effect in component.get("effects", {}).items():
                variable, level = term.rsplit("__", 1)
                unit, increment = UNITS[variable]
                beta, interval = interval_per_increment(effect, result["scales"][variable]["sd"], increment)
                role = "SI" if module == "bract_projection" else "Main"
                rows.append(dict(module=module, endpoint=endpoint, level=level,
                                 environment=variable, increment=increment, unit=unit,
                                 change=beta, lower=interval[0], upper=interval[1]))
                lines.append(f"|{role}|{endpoint}|{'種間' if level=='between' else '種内'}|"
                             f"{variable}|{increment:g} {unit}|{beta:+.4g}|"
                             f"[{interval[0]:+.4g}, {interval[1]:+.4g}]|")
    (output/"environment_effects_physical_units_ja.md").write_text("\n".join(lines)+"\n", encoding="utf-8")
    if rows:
        with (output/"environment_effects_physical_units.csv").open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
