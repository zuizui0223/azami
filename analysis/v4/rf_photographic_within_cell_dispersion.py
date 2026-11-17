"""Original-S1 photo-trait dispersion across 1/5-degree taxon-cells.

Exploratory diagnostic, not a camera measurement-error estimate, causal test,
Random Forest fit, or replacement for GEB's frozen S1/S2/H3c estimates.
Source S1 is decoded/checksum-authenticated inside a protected Actions job;
only coarse anonymous summaries can leave that job. No environmental source
complete-case gates are applied, so this is NOT the RF28/soil34 exact cohort.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from analysis.v4.native_integration import FEATURES

EXPECTED_N = 313350
GRAINS = (1, 5)
GROUP_KEYS = ("taxon", "lat_cell", "lon_cell")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(8 << 20), b''):
            h.update(block)
    return h.hexdigest()


def _normalise(v: np.ndarray) -> np.ndarray:
    v = np.asarray(v, dtype=float)
    norm = np.linalg.norm(v, axis=-1, keepdims=True)
    return np.divide(v, norm, out=np.full_like(v, np.nan), where=norm > 1e-8)


def _digest_distribution(values):
    vals = np.asarray(values, dtype=float)
    vals = vals[np.isfinite(vals)]
    if not len(vals):
        return {"n": 0, "q05": None, "q50": None, "q95": None}
    q = np.quantile(vals, [0.05, 0.5, 0.95])
    return {"n": int(len(vals)), "q05": float(q[0]),
            "q50": float(q[1]), "q95": float(q[2])}


def audit_construct(frame: pd.DataFrame, construct: str, *, degrees: int,
                    min_per_cell: int = 3, min_cells: int = 8,
                    min_photos_taxon: int = 40) -> dict:
    """Equal-cell weighted per-taxon dispersion; no private IDs in return.

    Within-cell RMS: mean over retained cells of the per-cell photo residual
    squared against the cell median (joint circular-chord loss for hue).
    Among-cell RMS: mean squared cell-median deviation from the within-taxon
    median of cell medians (joint circular-chord loss for hue). The ratio is
    not an ICC, a noise fraction, or an estimate of camera bias.
    """
    if construct not in FEATURES or degrees not in GRAINS or min_per_cell < 2:
        raise ValueError('INVALID_REGISTERED_DISPERSION_DESIGN')
    vals = list(FEATURES[construct])
    required = {"accepted_key", "latitude", "longitude", *vals,
                *(v + "__eligible" for v in vals)}
    if not required.issubset(frame):
        raise ValueError('DISPERSION_SOURCE_TRAIT_OR_QC_COLUMNS_MISSING')
    n = len(frame)
    lat = pd.to_numeric(frame.latitude, errors='coerce').to_numpy(float)
    lon = pd.to_numeric(frame.longitude, errors='coerce').to_numpy(float)
    y = frame[vals].apply(pd.to_numeric, errors='coerce').to_numpy(float)
    finite = (np.isfinite(y).all(axis=1) & np.isfinite(lat) & np.isfinite(lon)
              & (np.abs(lat) <= 90) & (np.abs(lon) <= 180)
              & frame.accepted_key.notna().to_numpy())
    for col in vals:
        finite &= frame[col + '__eligible'].eq(True).to_numpy(bool)
    if construct == 'hue':
        finite &= np.linalg.norm(y, axis=1) >= 0.20
    qualified_photo_count = int(finite.sum())
    if not qualified_photo_count:
        return {"construct": construct, "degrees": degrees,
                "status": "no_eligible_photographs", "n_original_rows_seen": n}
    good = pd.DataFrame({"taxon": frame.loc[finite,'accepted_key'].astype(str).to_numpy(),
                         "lat_cell": np.floor((lat[finite]+90)/degrees).astype(int),
                         "lon_cell": np.floor((lon[finite]+180)/degrees).astype(int)})
    for i,col in enumerate(vals):
        good[col] = y[finite,i]
    if construct == 'hue':
        unit = _normalise(good[vals].to_numpy(float))
        for i,col in enumerate(vals):
            good[col] = unit[:,i]
    grp = good.groupby(list(GROUP_KEYS), sort=False, observed=True)
    counts = grp.size().rename('n_photos')
    medians = grp[vals].median()
    cells = medians.join(counts).reset_index()
    cells = cells.loc[cells.n_photos >= min_per_cell].copy()
    if construct == 'hue':
        norm = np.linalg.norm(cells[vals].to_numpy(float),axis=1)
        cells = cells.loc[norm >= .20].copy()
    tax = cells.groupby('taxon').agg(n_cells=('n_photos','size'), n_photos=('n_photos','sum'))
    keep = tax.index[(tax.n_cells >= min_cells) & (tax.n_photos >= min_photos_taxon)]
    cells = cells.loc[cells.taxon.isin(keep)].copy()
    if cells.empty:
        return {"construct": construct, "degrees": degrees,
                "status": "insufficient_repeated_taxon_cells", "n_original_rows_seen": n,
                "n_eligible_photo_QC_before_cell_gate": qualified_photo_count,
                "n_qualified_taxa":0,"n_qualified_cells":0}
    joined = good.merge(cells[list(GROUP_KEYS)+vals], on=list(GROUP_KEYS),
                        how='inner', validate='many_to_one', suffixes=('','_median'))
    center = cells[vals].to_numpy(float)
    actual = joined[vals].to_numpy(float)
    centroid = joined[[v+'_median' for v in vals]].to_numpy(float)
    if construct == 'hue':
        center = _normalise(center)
        actual = _normalise(actual)
        centroid = _normalise(centroid)
    squared = np.sum((actual-centroid)**2,axis=1)
    residual = joined[list(GROUP_KEYS)].copy()
    residual['within_sq'] = squared
    per_cell = residual.groupby(list(GROUP_KEYS),observed=True)['within_sq'].mean().reset_index()
    per_taxon_within = per_cell.groupby('taxon')['within_sq'].mean()
    if construct == 'hue':
        for i,col in enumerate(vals):
            cells[col] = center[:,i]
    tax_center = cells.groupby('taxon')[vals].median()
    if construct == 'hue':
        tax_center.loc[:,:] = _normalise(tax_center.to_numpy(float))
    baseline = cells[['taxon']].join(tax_center, on='taxon',rsuffix='_taxmedian')
    ref = baseline[vals].to_numpy(float)
    if construct == 'hue':
        ref = _normalise(ref)
    baseline_sq = np.sum((center-ref)**2,axis=1)
    cells['between_sq'] = baseline_sq
    per_taxon_between = cells.groupby('taxon')['between_sq'].mean()
    within_rmse = np.sqrt(per_taxon_within.loc[keep].to_numpy(float))
    between_rmse = np.sqrt(per_taxon_between.loc[keep].to_numpy(float))
    ratios = np.divide(between_rmse, within_rmse,
                       out=np.full_like(between_rmse,np.nan),where=within_rmse>1e-12)
    finite_ratio = np.isfinite(ratios)
    return {
        "construct": construct, "degrees": int(degrees), "status":"estimated_QC_only",
        "n_original_rows_seen":int(n),
        "n_eligible_photo_QC_before_cell_gate":qualified_photo_count,
        "n_photo_records_in_supported_cells":int(cells.n_photos.sum()),
        "n_qualified_taxa":int(len(keep)),
        "n_qualified_cells":int(len(cells)),
        "photos_per_supported_taxon_cell":_digest_distribution(cells.n_photos.to_numpy(int)),
        "eligible_cells_per_taxon":_digest_distribution(cells.groupby('taxon').size().to_numpy(int)),
        "taxon_equal_within_cell_RMSE":_digest_distribution(within_rmse),
        "taxon_equal_among_cell_median_RMSE":_digest_distribution(between_rmse),
        "among_over_within_RMSE_per_taxon":_digest_distribution(ratios),
        "n_taxa_among_spread_below_within_spread":int((ratios[finite_ratio]<1).sum()),
        "n_finite_taxon_ratios":int(finite_ratio.sum()),
        "distance_unit":("circular_chord_on_unit_sin_cos" if construct=='hue'
                         else "original_photo_trait_numeric_units"),
        "not_photographic_measurement_error_or_reliability_ICC":True,
        "same_28_variable_complete_case_population_as_RF":False,
        "individual_taxon_ID_coordinates_photo_and_cell_rows_not_exported":True,
    }


def run(native:Path, receipt_path:Path, output:Path, *,expected_n:int=EXPECTED_N,
        constructs:tuple[str,...]|None=None, grains:tuple[int,...]=GRAINS):
    receipt=json.loads(Path(receipt_path).read_text(encoding='utf-8'))
    if (sha256(native)!=receipt.get('prepared_table_sha256') or
            receipt.get('native_observations')!=expected_n or
            receipt.get('environment_decoding_status')!=
            'physical_units_gsp_invalid_encoding_masked_v1'):
        raise ValueError('ORIGINAL_S1_DECODE_AND_SHA_NOT_AUTHENTICATED')
    chosen=tuple(FEATURES) if constructs is None else tuple(constructs)
    if not chosen or len(chosen)!=len(set(chosen)) or not set(chosen)<=set(FEATURES):
        raise ValueError('INVALID_CONSTRUCT_SELECTION')
    if not grains or len(grains)!=len(set(grains)) or not set(grains)<=set(GRAINS):
        raise ValueError('INVALID_GEO_GRAINS')
    required=['obs_id','accepted_key','latitude','longitude']
    for c in chosen:
        for trait in FEATURES[c]:
            required.extend([trait, trait+'__eligible'])
    raw=pd.read_csv(native, usecols=list(dict.fromkeys(required)),
                    dtype={'obs_id':str},low_memory=False)
    if len(raw)!=expected_n or not raw.obs_id.is_unique or raw.obs_id.isna().any():
        raise ValueError('SOURCE_S1_COUNT_OR_OBS_ID_DRIFT')
    raw=raw.drop(columns=['obs_id'])
    results=[audit_construct(raw,c,degrees=d) for c in chosen for d in grains]
    out={"schema":"azami_v4_original_s1_photo_within_cell_dispersion_QC_v1",
         "status":"source_authenticated_exploratory_diagnostic_not_model_fit",
         "source_original_s1_observations":int(expected_n),
         "source_original_decoded_sha256":sha256(native),
         "photo_QC_cohort_not_RF_28_or_soil_34_complete_case_cohort":True,
         "design":{"constructs":list(chosen),"degrees":list(grains),
                   "min_source_photos_per_cell":3,
                   "min_cells_per_taxon":8,"min_photos_per_taxon":40},
         "interpretation":"Within-cell photo dispersion mixes biological, temporal, and photographic variation; not measurement-error identification or RF predictability",
         "not_a_new_Holm_p_value_family_or_frozen_GEB_result":True,
         "private_row_taxon_observer_site_coordinates_not_exported":True,
         "results":results}
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(out,indent=2,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf-8')
    return out


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--native',type=Path,required=True)
    p.add_argument('--receipt',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    r=run(a.native,a.receipt,a.output)
    print(json.dumps({'status':r['status'],
                      'source_s1':r['source_original_s1_observations'],
                      'summary_rows':len(r['results']),
                      'constructs':r['design']['constructs']}),flush=True)

if __name__=='__main__':main()