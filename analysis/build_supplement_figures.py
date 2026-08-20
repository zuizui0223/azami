#!/usr/bin/env python3
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

plt.rcParams['svg.hashsalt'] = 'azami-supplement'

ENDPOINT_ORDER = [
    'Orientation','Lightness','Chroma','Hue sine','Hue cosine',
    'Aspect ratio','Circularity','Solidity','Width-profile CV'
]
LABELS = {
    'orientation_angle_degrees_median':'Orientation',
    'corolla_lab_lightness_median':'Lightness',
    'corolla_lab_chroma_median':'Chroma',
    'corolla_hue_sin_median':'Hue sine',
    'corolla_hue_cos_median':'Hue cosine',
    'shape_aspect_ratio_median':'Aspect ratio',
    'shape_circularity_median':'Circularity',
    'shape_solidity_median':'Solidity',
    'shape_width_cv_median':'Width-profile CV',
}
PRED_LABELS = {
    'chelsa_bio01':'BIO1','chelsa_bio04':'BIO4',
    'chelsa_bio12':'BIO12','chelsa_bio15':'BIO15'
}
SUPPORT_ENDPOINT = {
    'orientation_angle':'Orientation',
    'corolla_chroma':'Chroma',
    'hue_sin':'Hue sine',
    'hue_cos':'Hue cosine',
    'shape_aspect_ratio':'Aspect ratio',
}

def save(fig, out, stem):
    fig.tight_layout()
    fig.savefig(out / f'{stem}.svg', bbox_inches='tight', metadata={'Date': None})
    fig.savefig(out / f'{stem}.png', dpi=300, bbox_inches='tight', metadata={'Software':'azami'})
    fig.savefig(out / f'{stem}.pdf', bbox_inches='tight', metadata={'Creator':'azami','CreationDate':None,'ModDate':None})
    plt.close(fig)

def fig_s1(data, out):
    df = pd.read_csv(data / 'FigureS1_nested_variance.csv')
    df['label'] = pd.Categorical(df['label'], ENDPOINT_ORDER, ordered=True)
    df = df.sort_values('label')
    y = np.arange(len(df))
    fig, ax = plt.subplots(figsize=(8.4, 5.7))
    ax.errorbar(
        df['within_assigned_species_fraction'], y,
        xerr=[df['within_assigned_species_fraction'] - df['within_assigned_species_ci95_low'],
              df['within_assigned_species_ci95_high'] - df['within_assigned_species_fraction']],
        fmt='o', capsize=2, label='All assessable heads'
    )
    ax.scatter(df['one_head_per_photo_within_fraction'], y, marker='s', label='One head per photo')
    ax.errorbar(
        df['balanced_10_photo_within_median'], y,
        xerr=[df['balanced_10_photo_within_median'] - df['balanced_10_photo_ci95_low'],
              df['balanced_10_photo_ci95_high'] - df['balanced_10_photo_within_median']],
        fmt='^', capsize=2, label='10 photos per taxon'
    )
    ax.set_yticks(y, df['label'].astype(str))
    ax.set_xlim(0.45, 1.0)
    ax.set_xlabel('Fraction of visible variance below assigned-taxon mean')
    ax.set_title('Nested visible-variance sensitivities')
    ax.invert_yaxis()
    ax.legend(frameon=False, loc='upper left')
    save(fig, out, 'Figure_S1_nested_variance_sensitivities')

def expanded_primary(raw):
    rows = []
    for _, r in raw.iterrows():
        pred = PRED_LABELS[r['predictor']]
        if r['endpoint_type'] == 'linear':
            rows.append((LABELS[r['trait']], pred, float(r['beta_std_within'])))
        else:
            rows.append(('Hue sine', pred, float(r['hue_beta_sin'])))
            rows.append(('Hue cosine', pred, float(r['hue_beta_cos'])))
    return pd.DataFrame(rows, columns=['endpoint','predictor','beta']).drop_duplicates(['endpoint','predictor'])

def fig_s2(data, tables, out):
    raw = pd.read_csv(data / 'FigureS2_primary_coefficients.csv')
    df = expanded_primary(raw)
    support = pd.read_csv(tables / 'S03_primary_within_taxon_BH_supported_rows.csv')
    support['endpoint_label'] = support['endpoint'].map(SUPPORT_ENDPOINT)
    supported = set(zip(support['endpoint_label'], support['predictor']))
    preds = ['BIO1','BIO4','BIO12','BIO15']
    mat = df.pivot(index='endpoint', columns='predictor', values='beta').reindex(index=ENDPOINT_ORDER, columns=preds)
    vmax = np.nanmax(np.abs(mat.values))
    fig, ax = plt.subplots(figsize=(7.2, 6.3))
    im = ax.imshow(mat.values, vmin=-vmax, vmax=vmax, cmap='coolwarm', aspect='auto')
    ax.set_xticks(range(len(preds)), preds)
    ax.set_yticks(range(len(ENDPOINT_ORDER)), ENDPOINT_ORDER)
    for i, endpoint in enumerate(ENDPOINT_ORDER):
        for j, pred in enumerate(preds):
            v = mat.iloc[i, j]
            if np.isfinite(v):
                star = '*' if (endpoint, pred) in supported else ''
                ax.text(j, i, f'{v:+.3f}{star}', ha='center', va='center', fontsize=8,
                        fontweight='bold' if star else 'normal')
    fig.colorbar(im, ax=ax, label='Within-taxon standardized coefficient')
    ax.set_title('Primary within-taxon climate coefficients (36 component tests)')
    ax.text(0.0, -0.11, '* BH q < 0.05; hue sine/cosine require joint circular interpretation',
            transform=ax.transAxes, ha='left', va='top', fontsize=8)
    save(fig, out, 'Figure_S2_primary_coefficient_map')

def fig_s3(data, out):
    df = pd.read_csv(data / 'FigureS3_spde_model_group_summary.csv')
    df['endpoint'] = df['trait'].map(LABELS)
    groups = ['climate','climate_topography','climate_soil','full']
    group_labels = ['Climate','Climate + topography','Climate + soil','Full']
    mat = df.pivot(index='endpoint', columns='model_group', values='delta_waic_within_trait').reindex(index=ENDPOINT_ORDER, columns=groups)
    fig, ax = plt.subplots(figsize=(7.6, 6.2))
    im = ax.imshow(mat.values, aspect='auto', cmap='viridis')
    ax.set_xticks(range(len(groups)), group_labels, rotation=25, ha='right')
    ax.set_yticks(range(len(ENDPOINT_ORDER)), ENDPOINT_ORDER)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            v = mat.iloc[i, j]
            if np.isfinite(v):
                ax.text(j, i, f'{v:.1f}', ha='center', va='center', fontsize=8)
    fig.colorbar(im, ax=ax, label='ΔWAIC within endpoint (0 = best)')
    ax.set_title('Grouped SPDE-INLA model comparison')
    save(fig, out, 'Figure_S3_spde_delta_waic')

def ordered_spatial(data):
    df = pd.read_csv(data / 'FigureS45_spatial_robustness.csv')
    aliases = {
        'Corolla lightness':'Lightness', 'Corolla chroma':'Chroma',
        'Corolla hue sine':'Hue sine', 'Corolla hue cosine':'Hue cosine',
    }
    df['endpoint'] = df['endpoint'].replace(aliases)
    unknown = sorted(set(df['endpoint']) - set(ENDPOINT_ORDER))
    if unknown:
        raise ValueError(f'Unknown spatial-robustness endpoint labels: {unknown}')
    df['endpoint'] = pd.Categorical(df['endpoint'], ENDPOINT_ORDER, ordered=True)
    return df.sort_values('endpoint')

def fig_s4(data, out):
    df = ordered_spatial(data)
    y = np.arange(len(df))
    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    ax.scatter(df['residual_morans_I'], y, s=42)
    ax.axvline(0, linewidth=1, linestyle='--')
    ax.set_yticks(y, df['endpoint'].astype(str))
    ax.invert_yaxis()
    ax.set_xlabel("Residual Moran's I")
    ax.set_title('Residual spatial autocorrelation after frozen SPDE control')
    ax.text(0.0, -0.11, '999 permutations; no endpoint had P < 0.05',
            transform=ax.transAxes, ha='left', va='top', fontsize=8)
    save(fig, out, 'Figure_S4_residual_morans_I')

def fig_s5(data, out):
    df = ordered_spatial(data)
    y = np.arange(len(df))
    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    ax.scatter(df['minimum_leave_one_region_out_spearman_rho'], y, s=42)
    ax.set_yticks(y, df['endpoint'].astype(str))
    ax.invert_yaxis()
    ax.set_xlim(0.84, 1.0)
    ax.set_xlabel('Minimum leave-one-broad-region-out Spearman rho')
    ax.set_title('Taxon-rank stability under broad-region omission')
    save(fig, out, 'Figure_S5_leave_one_region_out_stability')

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--data-dir', type=Path, required=True)
    ap.add_argument('--tables-dir', type=Path)
    ap.add_argument('--out-dir', type=Path, required=True)
    a = ap.parse_args()
    tables = a.tables_dir or a.data_dir.parent / 'tables'
    a.out_dir.mkdir(parents=True, exist_ok=True)
    fig_s1(a.data_dir, a.out_dir)
    fig_s2(a.data_dir, tables, a.out_dir)
    fig_s3(a.data_dir, a.out_dir)
    fig_s4(a.data_dir, a.out_dir)
    fig_s5(a.data_dir, a.out_dir)
    expected = (
        list(a.out_dir.glob('Figure_S*.svg'))
        + list(a.out_dir.glob('Figure_S*.png'))
        + list(a.out_dir.glob('Figure_S*.pdf'))
    )
    if len(expected) != 15:
        raise SystemExit(f'Expected 15 figure files, got {len(expected)}')

if __name__ == '__main__':
    main()
