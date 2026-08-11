"""
Vikram Nathan, 03/30/2026

Find downsampled connectivity + distance to each epicentre AFTER filtering symmmetric indices

Prerequisites: 
- find_target_voxels_symmetric_to_source_voxels.py
- calculate_distance_mat_source_target_200um.py
- filter_all_SIR_inputs_symmetric_indices.py
"""
import os
import numpy as np
from pyminc.volumes.factory import *
import pyminc.volumes.factory as pyminc
import pandas as pd
import pickle
import time

working_dir="../"
downsampled_conn_filt_dir="../derivatives/SIR_inputs_200um/symmetric_source_target_masked/connectivity/"

def conn_to_volume(conn, template_path, output_path, integer_label_vol_200um, target_200um_indices):
    """
    Inputs: 
    - conn: connectome (n_vox x n_vox; 
    each row/column index is an integer whose 1D "flat" coordinate
     needs to be converted to a 3D coordinate)

    - template_path: template file whose format will be copied in pyminc.volumeLikeFile
    - output_path: path for the output volume
    - integer_label_vol_200um: integer labels for all 200um voxels 
    - target_200um_indices: integer labels representing points with connectivity in the brain


    """
    vol = pyminc.volumeLikeFile(template_path, output_path)
    
    x_dim, y_dim, z_dim = vol.data.shape
    total_voxels = x_dim * y_dim * z_dim

    # Build a flat array of zeros, then fill in only the indices present in conn
    flat = np.zeros(total_voxels).flatten()
    # Create a boolean mask for each target in conn's column order
    # assumption: target_200um_indices is SORTED
    flattened_int_vol=integer_label_vol_200um.flatten()
    mask = np.array([np.where(flattened_int_vol == target_200um_indices[i]) for i in range(len(target_200um_indices))]).flatten()
    
    print("Columns aligned?", np.all(target_200um_indices == conn.columns))
    flat[mask] = np.array(conn.values).flatten()  # or conn.value

    vol.data = flat.reshape(x_dim, y_dim, z_dim)
    vol.writeFile()
    vol.closeVolume()

if __name__ == "__main__":
    os.chdir(working_dir)

    start_time = time.time()
    source_target_indices=pd.read_pickle("../derivatives/downsampled_connectome/source_target_indices_filt_overlap.pkl")
    source_200um_indices = np.array(source_target_indices['source_indices']).flatten().astype(int)
    target_200um_indices = np.array(source_target_indices['target_indices']).flatten().astype(int)
    #print(target_200um_indices)
    integer_label_vol_200um = np.array(pyminc.volumeFromFile("../derivatives/downsampled_connectome/target_indices_filt_overlap_source_200um.mnc", dtype='int').data)
    connectivity_matrix = pd.read_pickle(downsampled_conn_filt_dir + "voxel_model_200um_full_source_target_filt.pkl")

    ##convert to int, if not already int
    connectivity_matrix.columns = np.array(connectivity_matrix.columns).astype(int)
    connectivity_matrix.index = np.array(connectivity_matrix.index).astype(int)

    ###assert columns and targets are equal
    #print(np.all(connectivity_matrix.columns == target_200um_indices))
    ##load source and target masks for indexing purposes
    end_time = time.time()

    print("Elapsed time to load connectivity matrix:", end_time - start_time)

    ##load epicentre connectivity; find 200um "label" overlaying each 100um epicentre
    convert_100um_200um_label="../derivatives/convert_100um_to_200um_mnc_labels/average_template_200_int_upsampled_100_nn.mnc"
    cp_epi_conn_file="../derivatives/epi_connections/pir_03202026/cp_label_allen_ccfv3_xyz.mnc"
    dg_epi_conn_file="../derivatives/epi_connections/pir_03202026/dg_label_allen_ccfv3_xyz.mnc"

    os.system("mincmath -mult -clobber " + convert_100um_200um_label + " " + cp_epi_conn_file + " /tmp/cp_epi_product.mnc")
    os.system("mincmath -mult -clobber " + convert_100um_200um_label + " " + dg_epi_conn_file + " /tmp/dg_epi_product.mnc")


    ##these integer indices are derived by overlaying the integer 200um file over the epicentre file
    ##repeat for DG epicentre
    
    cp_prod_vol=volumeFromFile("/tmp/cp_epi_product.mnc")
    cp_epi_index=np.round(np.mean(cp_prod_vol.data[cp_prod_vol.data != 0]))
    print("CP Epi Index (200um):", cp_epi_index)

    dg_prod_vol=volumeFromFile("/tmp/dg_epi_product.mnc")
    dg_epi_index=np.round(np.mean(dg_prod_vol.data[dg_prod_vol.data != 0]))
    print("DG Epi Index (200um):", dg_epi_index)

    cp_epi_conn = connectivity_matrix.loc[source_200um_indices == cp_epi_index, :]
    dg_epi_conn = connectivity_matrix.loc[source_200um_indices == dg_epi_index, :]

    template = "../preprocessed/templates/average_template_200_pir.mnc"

    conn_to_volume(cp_epi_conn, template, "../derivatives/SIR_inputs_200um/cp_epi_knox_downsample_200um.mnc", integer_label_vol_200um, target_200um_indices)
    conn_to_volume(dg_epi_conn, template, "../derivatives/SIR_inputs_200um/dg_epi_knox_downsample_200um.mnc", integer_label_vol_200um, target_200um_indices)


    ####repeat for distance
    connectivity_matrix_distance = pd.read_pickle(downsampled_conn_filt_dir + "voxel_model_200um_full_pairwise_distance_source_target_filt.pkl")
    cp_epi_conn_dist = connectivity_matrix_distance.loc[source_200um_indices == cp_epi_index, :]
    dg_epi_conn_dist = connectivity_matrix_distance.loc[source_200um_indices == dg_epi_index, :]

    conn_to_volume(cp_epi_conn_dist, template, "../derivatives/SIR_inputs_200um/cp_epi_dist_knox_downsample_200um.mnc", integer_label_vol_200um, target_200um_indices)
    conn_to_volume(dg_epi_conn_dist, template, "../derivatives/SIR_inputs_200um/dg_epi_dist_knox_downsample_200um.mnc", integer_label_vol_200um, target_200um_indices)

