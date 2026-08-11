"""
Vikram Nathan, 03/25/2026

#Prerequisite: find_target_voxels_symmetric_to_source_voxels.py

#Goal: filter all SIR inputs to define values at common target indices (sources + reflected sources), and zero for all other voxels
#Output:
###(1) volumes for each SIR input
###(2) .pkl files with values in the SAME ORDERING as target indices

Note: can be done for any subsetted list of source/target voxels from the Allen (including STR/HIP only voxels)
"""

import os
import numpy as np
from pyminc.volumes.factory import *
import pyminc.volumes.factory as pyminc
import pandas as pd
import pickle
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
import sys
import subprocess
from pathlib import Path
from os import listdir
from os.path import isfile, join

###integer mask volume containing ALL possible target indices 
mask_vol_path="../derivatives/downsampled_connectome/target_indices_filt_overlap_source_200um.mnc"

###NOTE: can also run for full target indices, with no prefix (prefix="")
##path "../derivatives/downsampled_connectome/source_target_indices_filt_overlap.pkl"
target_indices_paths=["../derivatives/downsampled_connectome/str_source_target_indices_filt_overlap.pkl",
                      "../derivatives/downsampled_connectome/hip_source_target_indices_filt_overlap.pkl",
                       "../derivatives/downsampled_connectome/source_target_indices_filt_overlap.pkl"]

#target_indices_paths=["../derivatives/downsampled_connectome/source_target_indices_filt_overlap.pkl"]
prefixes=["str","hip", ""]
#prefixes=['']
conn_path="../derivatives/downsampled_connectome/voxel_model_200um_full.pkl"
conn_distance_path="../derivatives/downsampled_connectome/voxel_model_200um_full_pairwise_distance.pkl"

output_dir="../derivatives/SIR_inputs_200um/symmetric_source_target_masked/"
atrophy_dir="../derivatives/voxel_atrophy_maps/"
ge_dir="../derivatives/gene_expression_kNN/CCFv3_PIR/" 
#connectivity_dir= ###save for when Trillium is back

N_WORKERS = min(32, os.cpu_count())

def makedirs_ok(*paths):
    for p in paths:
        Path(p).mkdir(parents=True, exist_ok=True)

# mask_data is a numpy array (not a pyminc volume) so it can be pickled by ProcessPoolExecutor
def filter_mnc_file_source_target_indices(input_vol_path, mask_data, output_vol_dir, output_pkl_dir, target_indices):

    input_vol = pyminc.volumeFromFile(input_vol_path)
    basename=Path(input_vol_path).stem
    output_path=output_vol_dir+basename+"_target_filt.mnc"
    output_vol = pyminc.volumeLikeFile(input_vol_path, output_path)
    output_vol.data = np.zeros(output_vol.data.shape)
    output_vol.data[mask_data > 0] = input_vol.data[mask_data > 0]

    output_vol.writeFile()
    output_vol.closeVolume()

    ###todo: figure out how to subset so that .pkl file has same order as combined target indices
    # With a different volume (same shape as mask_vol):
    results = {
        i: input_vol.data[np.where(mask_data == i)]
        for i in target_indices
    }

    result_values = np.array([results[i] for i in target_indices]).T  # shape: (len(target_indices), n_voxels_per_label)
    result_df = pd.DataFrame(result_values, columns=target_indices)

    with open(output_pkl_dir + basename + "_target_filt.pkl", 'wb') as f:
        pickle.dump(result_df, f)

    input_vol.closeVolume()

def filter_downsampled_connectome(conn_path, output_conn_dir, source_indices, target_indices):

    basename=Path(conn_path).stem
    conn_df = pd.read_pickle(conn_path)
    conn_df_filt=conn_df.loc[source_indices, target_indices]
    print(np.array(conn_df_filt).shape)

    with open(output_conn_dir + basename + "_source_target_filt.pkl", 'wb') as f:
        pickle.dump(conn_df_filt, f)

def parallel_filter_files(file_list, mask_data, output_vol_dir, output_pkl_dir, target_indices, n_workers=N_WORKERS):
    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        futures = {
            executor.submit(filter_mnc_file_source_target_indices, f, mask_data, output_vol_dir, output_pkl_dir, target_indices): f
            for f in file_list
        }
        for future in as_completed(futures):
            fpath = futures[future]
            try:
                future.result()
            except Exception as exc:
                print(f"[ERROR] {Path(fpath).name}: {exc}")

if __name__ == "__main__":

    ###load mask once in main process as plain numpy array
    _mask_vol = pyminc.volumeFromFile(mask_vol_path, dtype='int')
    mask_data = np.asarray(_mask_vol.data, dtype=np.int32)
    _mask_vol.closeVolume()

    for i in range(len(prefixes)):
        prefix = prefixes[i]
        target_indices_path = target_indices_paths[i]
        atrophy_pir_files=[atrophy_dir+"CP/M83PBS_vs_M83HuPff_t_stats_linear_allen200_PIR_CCFv3.mnc",
                                    atrophy_dir+"CP/M83PBS_vs_M83MsPff_t_stats_linear_allen200_PIR_CCFv3.mnc",
                                    atrophy_dir + "DG/M83PBS_vs_M83HuPff_hipp_t_stats_linear_allen200_PIR_CCFv3.mnc"]
        ge_pir_files=[ge_dir + f for f in listdir(ge_dir) if isfile(join(ge_dir, f))]

        ###0. load mask volume/target indices
        source_target_indices=pd.read_pickle(target_indices_path)
        source_indices=source_target_indices['source_indices'].flatten()
        target_indices=source_target_indices['target_indices'].flatten()

        ###1. connectivity (TODO - different helper function)
        output_conn_dir=output_dir+"connectivity/"+prefix+"/"
        makedirs_ok(output_conn_dir)
        filter_downsampled_connectome(conn_path, output_conn_dir, source_indices, target_indices)
        print("Done filtering connectivity")

        ###2. pairwise distance
        output_conn_dir=output_dir+"connectivity/"+prefix+"/"
        makedirs_ok(output_conn_dir)
        filter_downsampled_connectome(conn_distance_path, output_conn_dir, source_indices, target_indices)
        print("Done filtering connectivity (Euclidean distance)")

        ###3. atrophy
        output_vol_dir=output_dir+"atrophy/"+prefix+"/"
        output_pkl_dir=output_dir+"atrophy_pkl/"+prefix+"/"
        makedirs_ok(output_vol_dir, output_pkl_dir)
        parallel_filter_files(atrophy_pir_files, mask_data, output_vol_dir, output_pkl_dir, target_indices)
        
        ###4. GE
        ###glob list of files in CCFv3 directory
        output_vol_dir=output_dir+"GE/"+prefix+"/"
        output_pkl_dir=output_dir+"GE_pkl/"+prefix+"/"
        makedirs_ok(output_vol_dir, output_pkl_dir)
        parallel_filter_files(ge_pir_files, mask_data, output_vol_dir, output_pkl_dir, target_indices)


        ###4. pairwise distance matrix between all source/target indices - see calculate_distance_mat_source_target_200um.py