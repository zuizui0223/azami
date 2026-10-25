"""Verify six native/introduced floral-colour coefficients from frozen V4 science.

Only reads the immutable original aggregate (not private identifiers or photos);
does not fit models or conduct native-versus-introduced climate-niche comparisons.
"""
from __future__ import annotations
import argparse, csv, hashlib, json, math, zipfile
from pathlib import Path

ZIP_SHA="6e081dc081bb8ee2311a05ac5b7fcddd3f6309be2fa909e87d432789f4e0e205"
MEMBER="V4_S1_S2_four_vs_vif5_vif10_alltraits.json"
MEMBER_SHA="3738292e0a64da183dfe219c854c16cda199a7f9769f4928da8ab07657d60b70"
# Source-verified from 2026-10-24. Each tuple: taxon, range, visible trait,
# environmental covariate, n, standardized beta, raw P, BH P, Holm P.
EXPECTED=(
 ("Cirsium arvense","native","lightness","topo_elevation",
  22936,0.05626673194073065,0.141070055011146,0.35405817728287625,1.0),
 ("Cirsium arvense","introduced","lightness","topo_elevation",
  26171,-0.17373187423223405,1.0252616543029317e-05,
  0.0013123349175077526,0.0013123349175077526),
 ("Cirsium vulgare","native","chroma","chelsa_rsds_mean",
  17335,0.1368829067685122,0.0008325007823301402,0.017760016689709657,0.10239759622660725),
 ("Cirsium vulgare","introduced","chroma","chelsa_rsds_mean",
  27200,0.1345037874025455,0.00011470971522737719,
  0.006519919664304471,0.014568133833876904),
 ("Cirsium vulgare","native","lightness","chelsa_rsds_mean",
  17335,-0.04880796672586844,0.33030030767570273,0.5154122231858496,1.0),
 ("Cirsium vulgare","introduced","lightness","chelsa_rsds_mean",
  27200,-0.14280073287481068,0.00015281061713213603,
  0.006519919664304471,0.019254137758649138),
)
COORD={"lightness":"corolla_lab_lightness","chroma":"corolla_lab_chroma"}

def close(a,b,label):
    if not math.isclose(float(a),float(b),rel_tol=2e-7,abs_tol=1e-10):
        raise ValueError(f"{label}: source changed: {a} != {b}")

def run(artifact:Path, output:Path|None):
    data=artifact.read_bytes()
    if hashlib.sha256(data).hexdigest()!=ZIP_SHA:
        raise ValueError("Original V4 source aggregate ZIP SHA-256 mismatch")
    with zipfile.ZipFile(artifact) as z:
        raw=z.read(MEMBER)
    if hashlib.sha256(raw).hexdigest()!=MEMBER_SHA:
        raise ValueError("Original V4 full-model JSON SHA-256 mismatch")
    source=json.loads(raw)
    s2=source["multiplicity"]["S2/four_bio18"]
    if s2["n_tests"]!=128 or s2["n_Holm"]!=3:
        raise ValueError("Original S2 128/Holm-3 family changed")
    collected=[]
    for sp,st,trait,env,n,b,p,q,h in EXPECTED:
        source_key=sp+"/"+st+"/"+trait
        matching=[r for r in s2["rows"] if r["group"]==source_key and
                  r["coordinate"]==COORD[trait] and r["environment"]==env and
                  r["scale"]=="within_species_range"]
        if len(matching)!=1: raise ValueError("Missing/duplicate fixed S2 endpoint: "+source_key)
        row=matching[0]
        if row["n"]!=n:raise ValueError("Fixed sample size changed: "+source_key)
        for col,expected in (("beta_std",b),("p_original",p),("p_BH",q),("p_Holm",h)):
            close(row[col],expected,source_key+":"+col)
        ci=row["ci95_wald"]
        if len(ci)!=2 or not ci[0]<b<ci[1]:
            raise ValueError("Source confidence interval invalid")
        collected.append(dict(species=sp,range_status=st,trait=trait,exposure=env,n=n,
                              beta_std=b,ci_lower=ci[0],ci_upper=ci[1],
                              raw_p=p,bh_p=q,holm_p=h))
    allrows=s2["rows"]
    for st,n_expected in (("native",0),("introduced",3)):
        n_pass=sum(r["p_Holm"]<.05 for r in allrows if "/"+st+"/" in r["group"])
        if n_pass!=n_expected:raise ValueError("Frozen Holm distribution changed")
    if output is not None:
        output.parent.mkdir(parents=True,exist_ok=True)
        with output.open("w",encoding="utf-8",newline="") as fp:
            w=csv.DictWriter(fp,fieldnames=list(collected[0]))
            w.writeheader();w.writerows(collected)
    print(f"PASS: original fixed S2 family 128; {len(collected)} native/introduced colour rows; Holm 0 native / 3 introduced")

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--artifact",type=Path,required=True)
    p.add_argument("--output",type=Path)
    a=p.parse_args();run(a.artifact,a.output)

if __name__=="__main__":main()
