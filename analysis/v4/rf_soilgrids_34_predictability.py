"""Source-real, budgeted SoilGrids 0–5 cm Q0.5 add-on to an archived 28-variable RF.

This is a posthoc, prospective-within-scope, *sparse soil-pixel* exploratory
analysis: 1-degree eligible photo-trait taxon-cells are hash-capped at 14 per
taxon and three original photo positions per retained cell. Soil sampling is
outcome-blind within qualified cells; no photo values, IDs, coordinates, taxon
names or per-location soil pixels are exported. Requires private authenticated
S1 and the original/checked DEM-derived terrain field; does NOT modify the
existing GEB article, positive Holm results or previous RF analysis.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import time
from pathlib import Path
import numpy as np
import pandas as pd

from analysis.v4.native_integration import FEATURES
from analysis.v4.rf_archived_atmosphere import verify_original_1981_2010_airfields
from analysis.v4.rf_extended_environment import BIO9, SOIL6, TERRAIN5, verify_extended_source
from analysis.v4.rf_multitrait_predictability import (
    _valid_measurement, aggregate_endpoint, score_models, sha256, verified_source,
)
from analysis.v4.rf_soilgrids_private_pilot import VRT_PINNED_SHA256, PHYSICAL_BOUNDS
from analysis.v4.rf_soilgrids_webdav_probe import (
    CONVERSIONS, _download_vrt, official_url, parsed_vrt_metadata,
)
from analysis.v4.rf_worldclim_terrain import GROUPS, MODEL_SET

SEED = 'azami-s1-soilgrids-1deg-sparse-rf-20261031-v1'
MODEL28 = tuple(MODEL_SET['chelsa_vpd_wind23_plus_worldclim_terrain28'])
MODEL34 = MODEL28 + tuple(SOIL6)
MODELS = {'terrain28_soilmatched': MODEL28, 'terrain28_plus_soil34': MODEL34}
GROUPED_IMPORTANCE = {**{k: v for k, v in GROUPS.items() if k != 'soil'},
                      'soil_surface_isric_q050': tuple(SOIL6)}
CELL_KEYS = ['taxon', 'lat_cell', 'lon_cell']
MAX_CELLS_PER_TAXON = 14
PHOTO_POINTS_PER_CELL = 3
MIN_COMPLETE_SOIL_POINTS_PER_CELL = 2
MAX_ORIGINAL_SOURCE_POINTS = 3000
ALLOWED_CONSTRUCTS = ('angle',)  # first original-source stage, frozen before fit
assert len(MODEL28) == 28 and len(MODEL34) == 34 and len(set(MODEL34)) == 34


def stable_hash(value: str) -> str:
    return hashlib.blake2b((SEED + '|' + str(value)).encode(), digest_size=12).hexdigest()


def selected_source_positions(frame: pd.DataFrame, construct: str = 'angle', *,
                              max_cells_per_taxon: int = MAX_CELLS_PER_TAXON):
    """No outcome-based ranking; first fix 1-degree eligible cells, then hash-subset.

    Caller must supply the same 28 verified predictor values for both models.
    The cell photo medians are estimated over ALL 28-eligible original photos;
    only soil uses a fixed three-photo location sample per qualifying cell.
    """
    if construct not in ALLOWED_CONSTRUCTS or max_cells_per_taxon != 14:
        raise ValueError('UNDECLARED_SPARSITY_OR_CONSTRUCT')
    cells = aggregate_endpoint(frame, construct, predictors=MODEL28, degrees=1,
                               min_per_cell=3, min_cells=8, min_photos=40,
                               min_taxa=25)
    if cells.taxon.nunique() < 25 or cells.region.nunique() < 3:
        raise ValueError('SOIL_SOURCE_1DEG_TAXON_COHORT_TOO_SMALL')
    cells = cells.copy()
    cells['__rank'] = [stable_hash(f'{t}|{r}|{c}') for t, r, c in
                       zip(cells.taxon, cells.lat_cell, cells.lon_cell)]
    subset = (cells.sort_values(['taxon', '__rank'], kind='stable')
                   .groupby('taxon', sort=False).head(max_cells_per_taxon)
                   .drop(columns=['__rank']).reset_index(drop=True))
    # Duplicate aggregate_endpoint's endpoint-specific photo eligibility. No
    # colour/angle values enter the hash ranking, only a fixed QC mask.
    names = list(FEATURES[construct])
    z = frame[[*MODEL28, 'latitude', 'longitude', *names]].apply(
        pd.to_numeric, errors='coerce')
    valid = np.isfinite(z.to_numpy(dtype=float)).all(axis=1)
    valid &= z.latitude.between(-90, 90).to_numpy(bool)
    valid &= z.longitude.between(-180, 180).to_numpy(bool)
    valid &= frame.accepted_key.notna().to_numpy(bool)
    valid &= _valid_measurement(frame, names)
    photo = frame.loc[valid, ['obs_id', 'accepted_key', 'latitude', 'longitude']].copy()
    photo['taxon'] = photo.accepted_key.astype(str)
    photo['lat_cell'] = np.floor(photo.latitude.add(90).div(1)).astype(int)
    photo['lon_cell'] = np.floor(photo.longitude.add(180).div(1)).astype(int)
    photo = photo.merge(subset[CELL_KEYS], on=CELL_KEYS,
                        how='inner', validate='many_to_one')
    photo['__rank'] = photo.obs_id.map(stable_hash)
    chosen = (photo.sort_values(CELL_KEYS + ['__rank'], kind='stable')
                   .groupby(CELL_KEYS, sort=False).head(PHOTO_POINTS_PER_CELL))
    samples_per_cell = chosen.groupby(CELL_KEYS).size()
    if (len(samples_per_cell) != len(subset) or
            not samples_per_cell.eq(PHOTO_POINTS_PER_CELL).all() or
            chosen.obs_id.duplicated().any()):
        raise ValueError('SOIL_PHOTO_SITE_SELECTION_NOT_THREE_PER_QUALIFIED_CELL')
    sites = chosen[['obs_id', 'latitude', 'longitude']].drop_duplicates('obs_id')
    if len(sites) > MAX_ORIGINAL_SOURCE_POINTS:
        raise ValueError('SOILGRIDS_PREDECLARED_SOURCE_POINT_BUDGET_EXCEEDED')
    return subset, chosen[['obs_id', *CELL_KEYS]].copy(), sites.reset_index(drop=True), {
        'n_baseline_eligible_taxa': int(cells.taxon.nunique()),
        'n_baseline_eligible_cells': int(len(cells)),
        'n_outcome_blind_selected_taxa': int(subset.taxon.nunique()),
        'n_selected_taxon_cells': int(len(subset)),
        'n_selected_source_photo_positions': int(len(sites)),
        'photo_sites_per_selected_cell': PHOTO_POINTS_PER_CELL,
        'max_selected_cells_per_taxon': MAX_CELLS_PER_TAXON,
        'cell_selection_uses_outcome_magnitudes': False,
        'source_photo_identifiers_taxon_or_coordinates_exported': False,
    }


def soil_reader(url: str, latitude, longitude, *, max_seconds=750):
    """Source VRT /vsicurl reads native original tiles; abort a too-slow layer."""
    import rasterio
    from rasterio.warp import transform
    from rasterio.transform import rowcol
    start = time.monotonic()
    virtual = '/vsicurl?max_retry=1&retry_delay=1&list_dir=no&url=' + url
    with rasterio.Env(GDAL_DISABLE_READDIR_ON_OPEN='EMPTY_DIR',
                      GDAL_HTTP_MAX_RETRY='1', GDAL_HTTP_RETRY_DELAY='1',
                      GDAL_HTTP_TIMEOUT='15', GDAL_HTTP_CONNECTTIMEOUT='9',
                      GDAL_CACHEMAX=64):
        with rasterio.open(virtual) as src:
            if src.count != 1 or src.crs is None:
                raise ValueError('SOIL_SOURCE_RASTER_BAND_OR_CRS_INVALID')
            xs, ys = transform('EPSG:4326', src.crs,
                               np.asarray(longitude).tolist(),
                               np.asarray(latitude).tolist())
            # Dataset.index is a scalar convenience wrapper; rowcol accepts
            # the vector of projected point coordinates from the photo cohort.
            rr, cc = rowcol(src.transform, xs, ys)
            order = np.lexsort((np.asarray(cc)//256, np.asarray(rr)//256))
            arr = np.full(len(xs), np.nan, dtype=float)
            points = np.column_stack([xs, ys])
            for i, (idx, sampled) in enumerate(zip(order, src.sample(points[order], masked=True, indexes=1))):
                if i % 16 == 0 and time.monotonic() - start > max_seconds:
                    raise TimeoutError('SOILGRIDS_BUDGETED_REMOTE_TILE_EXTRACTION_TIMEOUT')
                v = sampled[0]
                if not np.ma.is_masked(v) and np.isfinite(float(v)):
                    arr[idx] = float(v)
            return arr


def sample_six_source_sites(sites: pd.DataFrame, *, fetcher=None, reader=None):
    """Six provenance-pinned layers at the exact SAME subset of 1° cell sites."""
    fetcher = _download_vrt if fetcher is None else fetcher
    reader = soil_reader if reader is None else reader
    if len(sites) < 25 or len(sites) > MAX_ORIGINAL_SOURCE_POINTS:
        raise ValueError('INVALID_PREDECLARED_SOIL_SAMPLE_SIZE')
    output = sites[['obs_id']].copy()
    receipt = {}
    for name in SOIL6:
        t0 = time.monotonic()
        code, divisor, unit = CONVERSIONS[name]
        url = official_url(name)
        payload, _ = fetcher(url)
        digest = hashlib.sha256(payload).hexdigest()
        if digest != VRT_PINNED_SHA256[name]:
            raise ValueError('SOILGRIDS_SIX_LAYER_VRT_SHA_MISMATCH_' + code)
        meta = parsed_vrt_metadata(payload, code)
        raw = np.asarray(reader(url, sites.latitude.to_numpy(float),
                                   sites.longitude.to_numpy(float)), dtype=float)
        if raw.shape != (len(sites),) or np.isinf(raw).any():
            raise ValueError('SOILGRIDS_SOURCE_PIXEL_VECTOR_INVALID_' + code)
        transformed = raw / divisor
        valid = np.isfinite(transformed)
        low, high = PHYSICAL_BOUNDS[name]
        if (valid & ((transformed < low) | (transformed > high))).any():
            raise ValueError('SOILGRIDS_SOURCE_PIXEL_PHYSICAL_BOUND_FAILURE_' + code)
        output[name] = transformed
        receipt[name] = {'vrt_sha256':digest,
                         'official_vrt_url':url,
                         'source_tile_bytes_sha_not_verified':True,
                         'depth':'0-5cm','statistic':'Q0.5','unit':unit,
                         'finite_original_site_count':int(valid.sum()),
                         'finite_fraction':float(valid.mean()),
                         'raster_tile_references':meta['raster_tile_references'],
                         'source_pixel_extraction_seconds':round(time.monotonic()-t0,2)}
        print(json.dumps({'property':code,'finite':int(valid.sum()),
                          'status':'genuine_soilgrids_pixel_extracted'}),flush=True)
    complete = output[list(SOIL6)].notna().all(axis=1)
    aggregate={'n_distinct_original_photo_locations_requested':int(len(sites)),
               'n_six_layer_complete_original_photo_locations':int(complete.sum()),
               'six_layer_complete_fraction':float(complete.mean()),
               'n_imputed_missing_soil_pixels':0,
               'source_soilgrids_original_tiles_not_byte_pinned':True,
               'variables':receipt}
    return output, aggregate


def attach_soil_to_cells(selected:pd.DataFrame, top_photos:pd.DataFrame,
                         soil_sites:pd.DataFrame, *, min_taxa=25):
    """Two or more of the THREE preselected photo sites must have all 6 soils."""
    if top_photos.obs_id.duplicated().any() or soil_sites.obs_id.duplicated().any():
        raise ValueError('SOURCE_SOIL_PHOTO_ID_GRAIN_DRIFT')
    p = top_photos.merge(soil_sites, on='obs_id', how='left', validate='one_to_one', indicator=True)
    if not p['_merge'].eq('both').all():
        raise ValueError('SOIL_SITE_IDENTITY_MISSING')
    p['__joint'] = np.isfinite(p[list(SOIL6)].to_numpy(dtype=float)).all(axis=1)
    usable = p.loc[p['__joint']].copy()
    grouped=usable.groupby(CELL_KEYS, sort=True)
    soil=grouped[list(SOIL6)].median().reset_index()
    soil['__ncomplete']=grouped.size().to_numpy(int)
    soil = soil.loc[soil['__ncomplete']>=MIN_COMPLETE_SOIL_POINTS_PER_CELL].copy()
    joined = selected.merge(soil, on=CELL_KEYS, how='inner', validate='one_to_one')
    support = joined.groupby('taxon').agg(n_cells=('photos','size'), n_photos=('photos','sum'))
    names=support.index[(support.n_cells>=8)&(support.n_photos>=40)]
    eligible=joined.loc[joined.taxon.isin(names)].reset_index(drop=True)
    audit={'n_qualified_cells_before_soil':int(len(selected)),
           'n_cells_at_least_2_joint_soil_samples':int(len(soil)),
           'n_cells_after_taxon_requalification':int(len(eligible)),
           'n_taxa_after_soil_requalification':int(eligible.taxon.nunique()),
           'n_original_source_photos_under_retained_cells':int(eligible.photos.sum()),
           'minimum_complete_six_soil_points_per_cell':MIN_COMPLETE_SOIL_POINTS_PER_CELL,
           'exclusion_uses_soil_missingness_not_image_trait_value':True}
    if eligible.taxon.nunique()<min_taxa or eligible.region.nunique()!=3:
        return eligible,{**audit,'status':'insufficient_taxa_or_sectors_after_soil'},False
    return eligible,{**audit,'status':'fit_eligible'},True


def verified_private_28(native,receipt,seven,seven_sha,nine,nine_sha,nine_manifest,
                        terrain,terrain_manifest,five_degree,metadata):
    ref=json.loads(Path(five_degree).read_text())
    if (ref.get('status')!='original_source_terrain_derived_and_evaluated' or
            ref.get('schema')!='azami_original_s1_chelsa_atmosphere23_vs_worldclim_terrain28_exploratory_v1' or
            ref.get('terrain_provenance',{}).get('mismatch_n')!=0 or
            not ref.get('frozen_original_GEB_S1_S2_H3c_unchanged')):
        raise ValueError('PREVIOUS_WORLDCLIM_DEM_V4_28_COEFFICIENT_SOURCE_NOT_VERIFIED')
    if (sha256(terrain)!=ref['terrain_provenance']['terrain_private_table_sha256'] or
            sha256(terrain_manifest)!=ref['terrain_provenance']['terrain_private_manifest_sha256']):
        raise ValueError('SOIL_DEM_SOURCE_DIGEST_MISMATCH')
    frame,source=verified_source(native,receipt,seven,seven_sha)
    if source['source_s1_sha256']!=ref['source']['source_s1_sha256']:
        raise ValueError('SOIL_S1_SOURCE_PREVIOUS_DEM_RUN_NOT_SAME')
    frame,_=verify_extended_source(frame,nine,nine_sha,nine_manifest,required_columns=BIO9)
    verify_original_1981_2010_airfields(frame,receipt,metadata)
    frame,_=verify_extended_source(frame,terrain,sha256(terrain),terrain_manifest,
                                   required_columns=TERRAIN5)
    return frame,source,ref['terrain_provenance']


def run(native,receipt,seven,seven_sha,nine,nine_sha,nine_manifest,terrain,
        terrain_manifest,five_degree,metadata,output, *,trees=80,repeats=2,
        fetcher=None,reader=None):
    if trees<8 or repeats<1:
        raise ValueError('SOIL_SOURCE_RF_FIXED_HYPERPARAMS_INVALID')
    frame,source,dem=verified_private_28(native,receipt,seven,seven_sha,nine,nine_sha,
                                         nine_manifest,terrain,terrain_manifest,
                                         five_degree,metadata)
    selected,photo,sites,selection=selected_source_positions(frame)
    soil,soil_receipt=sample_six_source_sites(sites,fetcher=fetcher,reader=reader)
    cells,support,eligible=attach_soil_to_cells(selected,photo,soil)
    results=[]
    if eligible:
        for scheme in ('unseen_taxa_balanced','unseen_longitude_sectors'):
            if scheme=='unseen_longitude_sectors' and cells.groupby('region').taxon.nunique().min()<5:
                results.append({'generalisation':scheme,'status':'insufficient_five_taxa_per_heldout_sector'})
                continue
            result=score_models(cells,'angle',scheme,trees=trees,repeats=repeats,
                                predictor_sets=MODELS,predictor_groups=GROUPED_IMPORTANCE)
            result['status']='source_real_sparse_soil_pilot_posthoc'
            results.append(result)
    report={
        'schema':'azami_v4_source_real_28v34_sparse_1degree_soilgrids_q050_v1',
        'status':'source_real_rf_completed_posthoc_sparse_soil' if results and eligible else 'source_real_soil_sparsity_insufficient_rf',
        'source':{'s1_original_n':source['s1_observations'],
                  's1_source_sha256':source['source_s1_sha256'],
                  'verified_DEM_original_elevation_mismatch_n':dem['mismatch_n']},
        'selection':selection,'soil_source':soil_receipt,'soil_cell_support':support,
        'model_predictors':{k:list(v) for k,v in MODELS.items()},
        'construct':'angle','geographic_grain_degree':1,
        'splits':['unseen_taxa_balanced','unseen_longitude_sectors'],
        'trees':trees,'permutation_repeats':repeats,
        'small_fixed_photo_site_sample_not_full_1degree_soil_coverage':True,
        'soilgrids_is_modelled_soil_not_measured_soil':True,
        'original_frozen_GEB_Holm_H3c_and_prior_RF_unchanged':True,
        'confidential_obs_ids_taxon_names_coordinates_model_weights_and_predictions_not_exported':True,
        'results':results,
    }
    output=Path(output);output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(report,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for item in ('native','receipt','archived-seven','archived-seven-sha',
                 'bio-nine','bio-nine-sha','nine-manifest','terrain-table',
                 'terrain-manifest','five-degree','metadata','output'):
        p.add_argument('--'+item,type=str if item.endswith('sha') else Path,required=True)
    p.add_argument('--trees',type=int,default=80)
    p.add_argument('--repeats',type=int,default=2)
    a=p.parse_args()
    out=run(a.native,a.receipt,a.archived_seven,a.archived_seven_sha,
            a.bio_nine,a.bio_nine_sha,a.nine_manifest,a.terrain_table,
            a.terrain_manifest,a.five_degree,a.metadata,a.output,
            trees=a.trees,repeats=a.repeats)
    print(json.dumps({'status':out['status'],
                      'original_s1':out['source']['s1_original_n'],
                      'tested_soil_cell_count':out['soil_cell_support']['n_qualified_cells_before_soil'],
                      'n_taxa_fit':out['soil_cell_support']['n_taxa_after_soil_requalification'],
                      'models_run':len(out['results'])}),flush=True)

if __name__=='__main__':
    main()