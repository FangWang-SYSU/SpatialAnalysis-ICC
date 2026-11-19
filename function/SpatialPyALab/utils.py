import numpy as np
import random
import os 
import pandas as pd
from typing import Optional, Tuple, Any

def compute_px_by_um(adata,
                     old_key ='spatial',
                     new_key ='spatial_manual_fine',
                     px_by_um_before = 26.6, 
                     n_pairs=50000, 
                     seed=42):
    np.random.seed(seed)
    indices = np.random.choice(adata.n_obs, size=n_pairs, replace=False)

    px_by_um_list = []

    for i in range(n_pairs):

        idx1 = indices[i]
        idx2 = (indices[(i + 1) % n_pairs])  

        p1_before = adata.obsm[old_key][idx1]
        p2_before = adata.obsm[old_key][idx2]
        p1_after = adata.obsm[new_key][idx1]
        p2_after = adata.obsm[new_key][idx2]

        d_px_before = np.linalg.norm(p2_before - p1_before)
        d_um = d_px_before / px_by_um_before
        d_px_after = np.linalg.norm(p2_after - p1_after)

        px_by_um_after = d_px_after / d_um
        px_by_um_list.append(px_by_um_after)

    px_by_um = np.mean(px_by_um_list)

    return px_by_um


def _mu_scale_from_fc_uns(adata, prefix: str = "fc", value_path=("puck_variables", "coord_by_um")) -> float:
    vals = []

    for k, v in (adata.uns or {}).items():
        if not isinstance(k, str) or not k.startswith(prefix):
            continue
        if v is None:
            continue
 
        candidates = v if isinstance(v, (list, tuple)) else [v]
        for item in candidates:
            if not isinstance(item, dict):
                continue
            cur = item
            ok = True
            for p in value_path:
                if isinstance(cur, dict) and p in cur:
                    cur = cur[p]
                else:
                    ok = False
                    break
            if not ok:
                continue
            try:
                x = float(cur)
                if np.isfinite(x):
                    vals.append(x)
            except Exception:
                pass

    return float(np.mean(vals))


def add_mu_scale(
    adata,
    px_by_um_csv_path: Optional[str] = None,
    old_key: str = "spatial",
    new_key: Optional[str] = None,
    n_pairs: int = 50000,
):
    if "mu_scale" not in adata.uns:
        adata.uns["mu_scale"] = {}

    mu_scale_before = None

    if px_by_um_csv_path and os.path.exists(px_by_um_csv_path):
        try:
            df = pd.read_csv(px_by_um_csv_path)
            if "px_by_um" not in df.columns:
                raise ValueError(f"CSV 'px_by_um' {px_by_um_csv_path}")
            mu_scale_before = float(df["px_by_um"].mean())
        except Exception as e:
            mu_scale_before = None
            csv_err = e
    else:
        csv_err = None

    if mu_scale_before is None:
        mu_scale_before = _mu_scale_from_fc_uns(adata, prefix="fc", value_path=("puck_variables", "coord_by_um"))

    mu_scale_before = round(mu_scale_before, 4)
    adata.uns["mu_scale"][old_key] = mu_scale_before

    if new_key:
        mu_scale_new = compute_px_by_um(
            adata,
            old_key=old_key,
            new_key=new_key,
            px_by_um_before=mu_scale_before,
            n_pairs=n_pairs,
        )
        mu_scale_new = round(float(mu_scale_new), 4)
        adata.uns["mu_scale"][new_key] = mu_scale_new

    return adata


import pandas as pd
def assign_label_to_hex(df_a, df_b, threshold=0.5,label='Label'):
    merged_df = pd.merge(df_a, df_b, left_on='barcode', right_on='Sequence', how='left')
    label_counts = merged_df.groupby('Hex_ID')[label].value_counts(normalize=True).unstack(fill_value=0)
    label_assignments = []
    for hex_id, counts in label_counts.iterrows():
        max_label = counts.idxmax()
        max_count = counts.max()
        if max_count >= threshold:
            label_assignments.append((hex_id, max_label))
        else:
            label_assignments.append((hex_id, 'conflict'))
    
    label_df = pd.DataFrame(label_assignments, columns=['Hex_ID', 'Assigned_Label'])
    merged_result = pd.merge(df_b, label_df, on='Hex_ID', how='left')
    
    return merged_result
