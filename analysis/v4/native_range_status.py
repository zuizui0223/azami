"""Native-range-first cohort rules; WCVP source semantics, not generic GBIF defaults.

The blank-establishment native convention applies ONLY to an existing WCVP
TDWG distribution row. A missing region record is UNLISTED, never native.
No trait, detector, environmental result or sampling abundance enters eligibility.
"""
from __future__ import annotations


def norm(value):
    return ' '.join(str(value or '').strip().split())


def distribution_status(record):
    if not norm(record.get('locationId')).startswith('TDWG:'):
        return 'uncertain'
    occurrence = norm(record.get('occurrenceStatus')).upper()
    means = norm(record.get('establishmentMeans')).upper()
    if occurrence in {'DOUBTFUL', 'ABSENT', 'PRESENT_DOUBTFUL'}:
        return 'uncertain'
    if means in {'', 'NATIVE'}:
        return 'native'
    return 'introduced' if means == 'INTRODUCED' else 'uncertain'


def resolve_exact(name, records):
    exact = [r for r in records if norm(r.get('canonicalName') or r.get('scientificName')).casefold() == name.casefold()]
    targets = {}
    for r in exact:
        accepted = norm(r.get('taxonomicStatus')).upper() == 'ACCEPTED' and not r.get('synonym', False)
        key = r.get('key') if accepted else r.get('acceptedKey')
        if key is not None:
            targets.setdefault(int(key), []).append(r)
    out = {'input_name': name, 'resolution_status': 'unresolved_taxon', 'accepted_key': None,
           'accepted_name': '', 'n_exact_records': len(exact), 'n_accepted_targets': len(targets)}
    if len(targets) == 1:
        key, rows = next(iter(targets.items()))
        names = sorted({norm(r.get('accepted') or r.get('acceptedScientificName') or r.get('scientificName')) for r in rows})
        out.update(resolution_status='resolved_unique_accepted_key', accepted_key=key,
                   accepted_name=names[0] if names else '')
    return out


def bool_value(value):
    text = str(value).strip().lower()
    return True if text in {'true', '1'} else False if text in {'false', '0'} else None


def inclusion_reason(status, coordinate_usable, obscured, captive):
    if bool_value(coordinate_usable) is not True:
        return 'coordinates_unusable'
    if bool_value(obscured) is not False:
        return 'obscured'
    flag = bool_value(captive)
    if flag is not False:
        return 'captive' if flag is True else 'captive_unknown'
    return 'include' if status == 'native' else 'range_' + str(status)


def classify_roster(roster, source):
    """Classify every upstream observation; no phenotype-dependent filtering.

    source contains resolutions, distributions and pinned TDWG GeoJSON.
    Boundary, polygon overlap, invalid or missing coordinates remain unknown.
    """
    import numpy as np
    import pandas as pd
    from shapely import points
    from shapely.geometry import shape
    from shapely.strtree import STRtree
    if roster.obs_id.isna().any() or roster.obs_id.duplicated().any():
        raise ValueError('The source roster must contain unique observation IDs')
    out = roster.copy().reset_index(drop=True)
    resolved = {x['input_name']: x for x in source['resolutions']}
    unknown = sorted(set(out.taxon_name) - set(resolved))
    if unknown:
        raise ValueError('Source-name roster not completely queried: ' + repr(unknown))
    distribution = {}
    for key, rows in source['distributions'].items():
        for record in rows:
            code = norm(record.get('locationId')).removeprefix('TDWG:').upper()
            distribution.setdefault((int(key), code), set()).add(distribution_status(record))
    features = source['geometry']['features']
    geometries = [shape(x['geometry']) for x in features]
    codes = [x['properties']['LEVEL3_COD'].upper() for x in features]
    lat = pd.to_numeric(out.latitude, errors='coerce').to_numpy(float)
    lon = pd.to_numeric(out.longitude, errors='coerce').to_numpy(float)
    valid = np.isfinite(lat) & np.isfinite(lon) & (abs(lat) <= 90) & (abs(lon) <= 180)
    valid_indices = np.flatnonzero(valid)
    pair = STRtree(geometries).query(points(lon[valid], lat[valid]), predicate='within')
    matches = {}
    for point, poly in zip(*pair):
        matches.setdefault(int(valid_indices[point]), set()).add(codes[poly])
    out['tdwg_level3_code'] = [next(iter(matches[i])) if len(matches.get(i, set())) == 1 else '' for i in range(len(out))]
    out['resolution_status'] = out.taxon_name.map(lambda n: resolved[n]['resolution_status'])
    out['accepted_key'] = out.taxon_name.map(lambda n: resolved[n]['accepted_key']).astype('Int64')
    out['accepted_name'] = out.taxon_name.map(lambda n: resolved[n]['accepted_name'])
    statuses = []
    for name, code, good in zip(out.taxon_name, out.tdwg_level3_code, valid):
        r = resolved[name]
        if not good:
            status = 'invalid_coordinates'
        elif r['resolution_status'] != 'resolved_unique_accepted_key':
            status = 'unresolved_taxon'
        elif not code:
            status = 'unmapped_tdwg'
        else:
            choices = distribution.get((r['accepted_key'], code), set())
            status = next(iter(choices)) if len(choices) == 1 else 'unlisted' if not choices else 'ambiguous'
        statuses.append(status)
    out['native_range_status'] = statuses
    out['eligibility_reason'] = [inclusion_reason(*r) for r in zip(statuses, out.coordinate_usable_for_environment, out.obscured, out.captive)]
    out['native_primary'] = out.eligibility_reason.eq('include')
    return out
