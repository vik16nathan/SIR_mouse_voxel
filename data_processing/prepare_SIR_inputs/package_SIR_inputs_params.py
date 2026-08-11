import numpy as np
import pandas as pd
import pickle
import os

base_dir = "../derivatives/SIR_inputs_200um/symmetric_source_target_masked/"
params_output_dir = "../derivatives/SIR_inputs_200um/symmetric_source_target_masked/params/"
downsampled_conn_dir = "../derivatives/downsampled_connectome/"
prefixes = ['str', 'hip', '']
percentiles = [20, 30, 40]

#####add new argument for whether log is true or false
log = True

###if log is true: add a 
def threshold_top_percent(mat, pct):
    """Keep top pct% of NONZERO weights (or min, if logged), preserve original values."""
    flat = mat.ravel()
    nonzero = flat[flat > np.min(flat)]

    if nonzero.size == 0:
        raise ValueError("Matrix has no nonzero entries.")

    thresh = np.percentile(nonzero, 100 - pct)
    out = mat.copy()
    out[out < thresh] = 0
    return out


def check_empty(mat):
    """Return masks of empty rows and columns."""
    row_empty_mask = (mat.sum(axis=1) == 0)
    col_empty_mask = (mat.sum(axis=0) == 0)
    return row_empty_mask, col_empty_mask


def enforce_no_empty_weighted_fast(thresh_mat, original_mat):
    """
    Fix empty rows/columns by restoring strongest original edge.
    Skips rows/cols that are entirely zero in original.
    """
    mat = thresh_mat.copy()

    row_empty_mask, col_empty_mask = check_empty(mat)
    rows_to_fix = np.where(row_empty_mask)[0]
    cols_to_fix = np.where(col_empty_mask)[0]

    # Precompute argmax once
    row_argmax = np.argmax(original_mat, axis=1)
    col_argmax = np.argmax(original_mat, axis=0)

    edges_added = 0

    # Fix empty rows
    for r in rows_to_fix:
        if np.all(original_mat[r, :] == 0):
            continue
        c = row_argmax[r]
        if mat[r, c] == 0:
            mat[r, c] = original_mat[r, c]
            edges_added += 1

    # Fix empty cols
    for c in cols_to_fix:
        if np.all(original_mat[:, c] == 0):
            continue
        r = col_argmax[c]
        if mat[r, c] == 0:
            mat[r, c] = original_mat[r, c]
            edges_added += 1

    return mat, len(rows_to_fix), len(cols_to_fix), edges_added


def make_filename(prefix, pct, retro=False, log=False):
    base = "source_target_indices_filt_overlap"

    if prefix != '':
        base += f"_{prefix}"

    base += f"_top{pct}"

    if retro:
        base += "_retro"
    
    if log:
        base += "_log"

    return base + ".pkl"


if __name__ == "__main__":

    os.makedirs(params_output_dir, exist_ok=True)

    for prefix in prefixes:

        # Load base data
        weights = pd.read_pickle(base_dir + "connectivity/" + prefix + "/voxel_model_200um_full_source_target_filt.pkl")
        eps = 1e-12
        if log == True:
            weights = np.log(weights + eps)
        distance = pd.read_pickle(base_dir + "connectivity/" + prefix + "/voxel_model_200um_full_pairwise_distance_source_target_filt.pkl")

        if prefix == '':
            source_target_pkl = pd.read_pickle(downsampled_conn_dir + "source_target_indices_filt_overlap.pkl")
        else:
            source_target_pkl = pd.read_pickle(downsampled_conn_dir + prefix + "_source_target_indices_filt_overlap.pkl")

        W = weights.values  # faster

        for pct in percentiles:

            print(f"\nProcessing prefix='{prefix}' | top {pct}%")

            # ---------- ANTERO ----------
            params = dict()
            params['distance'] = distance
            params['sources'] = source_target_pkl['source_indices']
            params['targets'] = source_target_pkl['target_indices']
            params['Snca'] = pd.read_pickle(base_dir + "GE_pkl/" + prefix + "/Snca_coronal_target_filt.pkl")

            W_thr = threshold_top_percent(W, pct)

            row_empty, col_empty = check_empty(W_thr)

            if row_empty.any() or col_empty.any():
                W_fixed, n_rows, n_cols, edges_added = enforce_no_empty_weighted_fast(W_thr, W)
                changed = not np.array_equal(W_thr, W_fixed)

                print(f"  ANTERO: FIXED | rows_fixed={n_rows}, cols_fixed={n_cols}, edges_added={edges_added}, changed={changed}")

                W_thr = W_fixed
            else:
                print("  ANTERO: UNCHANGED | no empty rows/cols")

            params['weights'] = pd.DataFrame(W_thr, index=weights.index, columns=weights.columns)

            out_name = make_filename(prefix, pct, retro=False, log=log)
            with open(params_output_dir + out_name, 'wb') as f:
                pickle.dump(params, f)

            # ---------- RETRO ----------
            n_src = len(params['sources'])
            conn_ipsi = W[:, :n_src]
            conn_contra = W[:, n_src:]
            conn_retro = np.hstack((conn_ipsi.T, conn_contra.T))

            W_thr_retro = threshold_top_percent(conn_retro, pct)

            row_empty, col_empty = check_empty(W_thr_retro)

            if row_empty.any() or col_empty.any():
                W_fixed_retro, n_rows, n_cols, edges_added = enforce_no_empty_weighted_fast(
                    W_thr_retro, conn_retro
                )
                changed = not np.array_equal(W_thr_retro, W_fixed_retro)

                print(f"  RETRO: FIXED | rows_fixed={n_rows}, cols_fixed={n_cols}, edges_added={edges_added}, changed={changed}")

                W_thr_retro = W_fixed_retro
            else:
                print("  RETRO: UNCHANGED | no empty rows/cols")

            params_retro = params.copy()
            params_retro['weights'] = pd.DataFrame(
                W_thr_retro,
                index=weights.index,
                columns=weights.columns
            )

            out_name_retro = make_filename(prefix, pct, retro=True, log=log)
            with open(params_output_dir + out_name_retro, 'wb') as f:
                pickle.dump(params_retro, f)