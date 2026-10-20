"""Two-dimensional anatomy measurements; no inferred defence or spine counts."""
import numpy as np

def involucre(mask,base,rim,phyllaries=()):
    """Explicit anatomical body silhouette and annotated free tips; 2D only."""
    if mask is None or base is None or rim is None:return {'status':'missing_anatomical_annotation'}
    yy,xx=np.nonzero(mask)
    if len(xx)<10:return {'status':'insufficient_mask_pixels'}
    axis=np.asarray(rim,float)-np.asarray(base,float);length=np.linalg.norm(axis)
    if length<=1e-8:return {'status':'degenerate_axis'}
    axis/=length;perp=np.array([-axis[1],axis[0]])
    xy=np.column_stack([xx,yy])-base;along=xy@axis;across=xy@perp
    height=float(np.quantile(along,.98)-np.quantile(along,.02))
    width=float(np.quantile(across,.98)-np.quantile(across,.02))
    if width<=0:return {'status':'degenerate_width'}
    return {'status':'measured_2d_anatomical_body','body_length_width_ratio':height/width,
            'body_width_pixels':width,'visible_annotated_phyllaries':len(phyllaries),
            'phyllary_length_relative_to_body_width':[float(np.linalg.norm(np.asarray(b,float)-a)/width) for a,b in phyllaries],
            'botanical_spine_count_estimated':False}
