"""Trait-blind crossed cases draws with full-source component closure."""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
import pandas as pd


def _labels(values, name):
    values = np.asarray(values)
    if values.ndim != 1 or pd.isna(values).any():
        raise ValueError(f'Missing or malformed {name}')
    values = values.astype(str)
    if np.any(values==''):
        raise ValueError(f'Empty {name}')
    return values


def _closure(levels, components):
    names, codes = np.unique(levels, return_inverse=True)
    parent = np.arange(len(names))
    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i
    first = {}
    for code, component in zip(codes, components):
        if component in first:
            a,b = find(code),find(first[component])
            parent[max(a,b)] = min(a,b)
        else:
            first[component] = code
    closed = np.array([find(code) for code in codes])
    return np.unique(closed, return_inverse=True)[1]


@dataclass(frozen=True)
class SourcePartition:
    observation_ids: np.ndarray
    taxa: np.ndarray
    components: np.ndarray
    taxon_blocks: np.ndarray
    spatial_blocks: np.ndarray
    grid_degrees: int


def source_partition(observation_ids, taxa, components, latitude, longitude, *, grid_degrees):
    """Build BEFORE endpoint filtering, using source identities/coordinates only."""
    ids, taxa, components = (_labels(v,n) for v,n in
                            [(observation_ids,'observation IDs'),(taxa,'taxa'),(components,'components')])
    lat,lon = np.asarray(latitude,float),np.asarray(longitude,float)
    if (grid_degrees not in (2,5) or len(ids)==0 or len(np.unique(ids))!=len(ids)
            or len(taxa)!=len(ids) or len(components)!=len(ids)
            or lat.shape!=ids.shape or lon.shape!=ids.shape
            or not np.isfinite(lat).all() or not np.isfinite(lon).all()
            or np.any(abs(lat)>90) or np.any(abs(lon)>180)):
        raise ValueError('Source identities, grid or coordinates invalid')
    lat_bin = np.floor(np.minimum(lat+90,np.nextafter(180.,-np.inf))/grid_degrees).astype(int)
    wrapped = (lon+180)%360
    wrapped[abs(lat)==90] = 180
    lon_bin = np.floor(wrapped/grid_degrees).astype(int)
    cells = lat_bin*(360//grid_degrees)+lon_bin
    arrays = [ids,taxa,components,_closure(taxa,components),_closure(cells,components)]
    for array in arrays:
        array.flags.writeable = False
    return SourcePartition(*arrays,int(grid_degrees))


def cohort_partition(source, positions):
    positions = np.asarray(positions)
    if (positions.ndim!=1 or len(positions)==0 or positions.dtype.kind not in 'iu'
            or len(np.unique(positions))!=len(positions) or np.any(positions<0)
            or np.any(positions>=len(source.observation_ids))):
        raise ValueError('Use distinct aligned positions in the full source')
    t,s = source.taxon_blocks[positions],source.spatial_blocks[positions]
    if len(np.unique(t))<2 or len(np.unique(s))<2:
        raise ValueError('A dependence factor collapsed; do not substitute row resampling')
    return positions,t,s


def crossed_draw(source, positions, *, seed, replicate):
    """Return indices into the supplied cohort and refit taxon-copy labels.

    Full-source closures are computed once; filtering cannot break a component
    bridge. No independent observation draw, trait value or fitted weight enters.
    """
    positions,t,s = cohort_partition(source,positions)
    if any(not isinstance(v,(int,np.integer)) or isinstance(v,bool) or v<0 for v in (seed,replicate)):
        raise ValueError('Nonnegative integer seed and replicate required')
    rng = np.random.default_rng(np.random.SeedSequence([seed,replicate,source.grid_degrees]))
    tl,sl = np.unique(t),np.unique(s)
    taxon_draws = rng.choice(tl,len(tl),replace=True)
    space_draws = rng.choice(sl,len(sl),replace=True)
    counts = {level:int((space_draws==level).sum()) for level in sl}
    _,taxon_codes = np.unique(source.taxa,return_inverse=True)
    parts,labels = [],[]
    for copy,level in enumerate(taxon_draws):
        indices = np.flatnonzero(t==level)
        repeats = np.array([counts[s[i]] for i in indices])
        selected = np.repeat(indices,repeats)
        parts.append(selected)
        # Encode original taxon as an integer, avoiding label concatenation aliases.
        labels.extend([f'{copy}:{taxon_codes[positions[i]]}' for i in selected])
    selected = np.concatenate(parts)
    if not len(selected):
        raise ValueError('Empty crossed draw; retain this failure without redraw')
    return selected,np.asarray(labels)
