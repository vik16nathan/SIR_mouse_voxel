"""
By: Vikram Nathan, 03/20/2026

Goal: Find connectivity to epicentres, with PIR voxel coordinates determined using 
- ../derivatives/epi_connnections/dg_label_allen_ccfv3_xyz.mnc
- ../derivatives/epi_connections/cp_label_allen_ccfv3_xyz.mnc

The voxel coordinates are outputs in Allen PIR coordinates following multi-stage registration in ./paxinos_to_DSURQE_to_allen/resample.sh.

Usage: find_knox_conn_to_epi.py z y x outfile_strength (NOTE: default coordinate order using mincinfo/Display is z, y x)
Usage for CP:

python find_knox_conn_to_epi.py 51 34 79 \
    cp_epi_conn_strength_100um.pkl
     
Usage for DG:

python find_knox_conn_to_epi.py 78 21 79 \
     dg_epi_conn_strength_100um.pkl

"""
import numpy as np
import sys
import os
import matplotlib.pyplot as plt
import SimpleITK as sitk
import pickle
from allensdk.core.mouse_connectivity_cache import MouseConnectivityCache
from mcmodels.core import VoxelModelCache
from mcmodels.models.voxel import VoxelConnectivityArray
from allensdk.core import json_utilities
import argparse


##prerequisites: module load minc-toolkit-v2
##path variables - CHANGE THESE based on environment ###########
#should contain nodes.csv.gz and weights.csv.gz obtained by rerunning build_model.py from mouse_connectivity_models
nodes_weights_dir="/scratch/vnathan/copy_of_allen_qc_reproducible_02212026/allen_connectome_qc/mouse_connectivity_models/paper/connectivity/voxel-standard-model/"

##Working directory for all input/output files. Should contain: 
# -  ./connectivity/voxel_model_manifest.json
# -  ./inputs/templates/average_template_100_ccfv3.mnc
voxel_model_manifest_dir="/scratch/vnathan/copy_of_allen_qc_reproducible_02212026/allen_connectome_qc/mouse_connectivity_models/paper/connectivity/"
output_dir='../derivatives/epi_connections/pir_03202026/'

#figure out where minc_toolkit is installed 
os.environ["MINC_TOOLKIT"] = "/opt/quarantine/software/minc-toolkit-v2/1.9.18.2/install"

#################################################################

def find_indices_target(coord_x,coord_y,coord_z, mask):
    indices_x=np.array(np.where(np.isin(mask[:,0],[coord_x])))
    indices_y=np.array(np.where(np.isin(mask[:,1],[coord_y])))
    indices_z=np.array(np.where(np.isin(mask[:,2],[coord_z])))
    indices_one=np.intersect1d(indices_x,indices_y)
    indices_both=np.intersect1d(indices_one,indices_z)
    return(list(indices_both))

if __name__ == "__main__": 
    ##parse arguments
    parser = argparse.ArgumentParser(description='Find modified Knox connectivity to epicentre coordinate (in CCFv3)')
    parser.add_argument('z_epi', type=int, help='Z voxel coordinate of epicentre (CCFv3, PIR, 100 microns)')
    parser.add_argument('y_epi', type=int, help='Y voxel coordinate of epicentre (CCFv3, PIR, 100 microns)')
    parser.add_argument('x_epi', type=int, help='X voxel coordinate of epicentre (CCFv3, PIR, 100 microns)')
    parser.add_argument('outfile_strength', type=str, help='Output file (conn. strength)')

    args = parser.parse_args()

    z_epi = args.z_epi
    y_epi = args.y_epi
    x_epi = args.x_epi
    outfile = args.outfile_strength

    ##load connectivity (old)
    mcc = MouseConnectivityCache(resolution=100) ###cannot change to anything but 100 
    annot, annot_info = mcc.get_annotation_volume()

    cache = VoxelModelCache(manifest_file=voxel_model_manifest_dir+'voxel_model_manifest.json')

    _, source_mask, target_mask = cache.get_voxel_connectivity_array()
    source_mask_local=source_mask.coordinates
    target_mask_local=target_mask.coordinates

    ###substitute in rebuilt nodes/weights with "corrected" number of experiments/centroids
    nodes_standard = np.loadtxt(nodes_weights_dir+'nodes_rebuilt.csv.gz', delimiter=',')
    weights_standard = np.loadtxt(nodes_weights_dir+'weights_rebuilt.csv.gz', delimiter=',')

    voxel_array = VoxelConnectivityArray(weights_standard, nodes_standard)
    knox_conn_index=find_indices_target(z_epi, y_epi, x_epi, source_mask_local) ###note: this is the SOURCE index (be careful with indexing)
    ###note that the coordinates are swapped in the Knox connectome compared to our CCFv3 epicentre (double-check by visualizing outputs)

    conn_strength = voxel_array[knox_conn_index,:]

    ###write out the connection strength as an intermediate file
    with open(output_dir+outfile, 'wb') as f:
        pickle.dump(conn_strength, f)
    
