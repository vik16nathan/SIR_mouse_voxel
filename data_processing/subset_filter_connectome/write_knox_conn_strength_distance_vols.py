"""
By: Vikram Nathan, 03/20/2026

Goal: Write connectivity to epicentres, given .pkl files with connection strength to each epicentre
Prerequisite: find_knox_conn_to_epi.py (and epicentre x,y,z coords)

Usage: write_knox_conn_strength_distance_vols.py z y x infile_strength outvol_strength outvol_distance

Example Usage for DG:
 
python write_knox_conn_strength_distance_vols.py 78 21 79 \
     dg_epi_conn_strength_100um.pkl \
     dg_epi_conn_strength_100um_ccfv3_annot.mnc \
     dg_epi_conn_distance_100um_ccfv3_annot.mnc

For CP:
python write_knox_conn_strength_distance_vols.py 51 34 79 \
     cp_epi_conn_strength_100um.pkl \
     cp_epi_conn_strength_100um_ccfv3_annot.mnc \
     cp_epi_conn_distance_100um_ccfv3_annot.mnc


"""

import numpy as np
import sys
import pandas as pd
import os
import matplotlib.pyplot as plt
import SimpleITK as sitk
import pickle
import pyminc.volumes.factory as pyminc
from pyminc.volumes.factory import volumeFromFile
import argparse

epi_conn_dir = "../derivatives/epi_connections/pir_03202026/"
output_dir='../derivatives/epi_connections/pir_03202026/'

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Write volumes for voxel-level connectivity to epicentre coordinate (in CCFv3)')
    parser.add_argument('z_epi', type=int, help='Z voxel coordinate of epicentre (CCFv3, PIR, 100 microns)')
    parser.add_argument('y_epi', type=int, help='Y voxel coordinate of epicentre (CCFv3, PIR, 100 microns)')
    parser.add_argument('x_epi', type=int, help='X voxel coordinate of epicentre (CCFv3, PIR, 100 microns)')
    parser.add_argument('infile_strength', type=str, help='Input file (conn. strength)')
    parser.add_argument('outvol_strength', type=str, help='Output volume (conn. strength)')
    parser.add_argument('outvol_distance', type=str, help='Output volume (conn. distance)')

    args = parser.parse_args()

    z_epi = args.z_epi
    y_epi = args.y_epi
    x_epi = args.x_epi
    infile_strength = args.infile_strength
    outvol_strength = args.outvol_strength
    outvol_dist = args.outvol_distance


    with open(epi_conn_dir + infile_strength, 'rb') as file:
        conn_strength=pickle.load(file)
    
    ##load source and target masks (lists of voxel coordinates, in 100 um PIR space, for each source/target index)
    source_target_masks_local = pd.read_pickle("../derivatives/source_target_masks_local/source_target_masks_local_100um.pkl")
    source_mask_local = source_target_masks_local['source_mask_local']
    target_mask_local = source_target_masks_local['target_mask_local']

    
    epi_conn_strength_vol = pyminc.volumeLikeFile("../preprocessed/templates/average_template_100_ccfv3.mnc", output_dir+outvol_strength)
    conn_strength = conn_strength.flatten()
    for i in range(target_mask_local.shape[0]):
        x=target_mask_local[i,0]
        y=target_mask_local[i,1]
        z=target_mask_local[i,2]
        epi_conn_strength_vol.data[z, y, x] = conn_strength[i]

    epi_conn_strength_vol.writeFile()

    ###Repeat for connection distances; start with Euclidean distance###
    epi_conn_distance_vol = pyminc.volumeLikeFile("../preprocessed/templates/average_template_100_ccfv3.mnc", output_dir+outvol_dist)
    epi_world_coords = epi_conn_distance_vol.convertVoxelToWorld([z_epi, y_epi, x_epi]) ###epi coords are flipped relative to x,y,z target coords - see prerequisite
    
    for i in range(target_mask_local.shape[0]):
        x=target_mask_local[i,0]
        y=target_mask_local[i,1]
        z=target_mask_local[i,2]

        ###get world coordinates
        world_coords = epi_conn_distance_vol.convertVoxelToWorld([x, y, z])
        epi_conn_distance_vol.data[z, y, x] = np.linalg.norm(world_coords - epi_world_coords)

    epi_conn_distance_vol.writeFile()
