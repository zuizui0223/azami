from pathlib import Path
from analysis.v5.native_pca_env_fullcoverage import prepare_original_wild_sources,apply_original_five_CHELSA
from analysis.v5.native_pca_env import background_mesh,_background_species_audit,_pca_env_2d,_probability_content
from analysis.v4.native_cached_climate import VARIABLES
import json,os
import numpy as np
import pandas as pd

def main():
    root=Path('work/submission'); private=root/'private'; out=root/'overlap';out.mkdir(exist_ok=True)
    wild,auth=prepare_original_wild_sources(private/'original-metadata.aesgcm',private/'authority350.zip',os.environ['AZAMI_METADATA_KEY'])
    land=background_mesh(auth)
    wild,land,_=apply_original_five_CHELSA(wild,land,source_meta=Path('analysis/v4/environment_source_metadata.json'),work=private/'rasters')
    frozen=json.loads(Path('analysis/v4/results/20261021_unified_abcde/S2_TWO_SPECIES_FULL_WCVP_NICHE_VERIFIED_20261021.json').read_text())
    for sp,sites in wild.dropna(subset=list(VARIABLES)).groupby('taxon_name'):
        env,_=_background_species_audit(auth,land.dropna(subset=list(VARIABLES)),sp,land)
        _,edges,a,b=_pca_env_2d(env[env._range=='native'],env[env._range=='introduced'],sites[sites.range_state=='native'],sites[sites.range_state=='introduced'],return_grid=True)
        native95=_probability_content(a); intro95=_probability_content(b)
        f=frozen['species'][sp]
        assert np.isclose(np.minimum(a,b).sum(),f['PCA_env_background_corrected_Schoener_D'])
        assert np.isclose(b[~native95].sum(),f['introduced_density_outside_native_95percent_probability_content'])
        np.savez_compressed(out/(sp.replace(' ','_')+'.npz'),x=edges[0],y=edges[1],native=a,introduced=b,native95=native95,introduced95=intro95,variance=f['PCA_env_2axis_variance_fraction'])
if __name__=='__main__': main()
