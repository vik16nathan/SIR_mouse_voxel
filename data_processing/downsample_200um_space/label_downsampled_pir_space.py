"""
GOAL: create integer labels for each voxel of the 200 um downsampled PIR/CCFv3 space, 
map 100 um --> 200 um voxels 

Prerequisites: pickle_source_target_masks_local.py (need to run on compute environment w/ WiFi to load VoxelModelCache)
Create average_template_200_int_upsampled_100_nn.mnc using commented-out lines in this script

"""
import os
import numpy as np
from pyminc.volumes.factory import *
import pyminc.volumes.factory as pyminc
import pickle

BASE_DIR="/scratch/vnathan/sir_voxel/"
if __name__ == "__main__":

    os.chdir(BASE_DIR)

    ###create volume like 200 um template RESAMPLED TO INT
    infile="preprocessed/templates/average_template_200_pir.mnc"
    output="preprocessed/templates/average_template_200_pir_int.mnc"

    os.system("mincresample -int " + infile + " -clobber " + output)

    template_ccfv3_int_vol = pyminc.volumeLikeFile("preprocessed/templates/average_template_200_pir_int.mnc", "derivatives/convert_100um_to_200um_mnc_labels/average_template_200_ccfv3_int_labels.mnc")
    template_ccfv3_int_vol.data = template_ccfv3_int_vol.data.astype(int)

    ###write ascending integer elements
    x_dim = template_ccfv3_int_vol.data.shape[0]
    y_dim = template_ccfv3_int_vol.data.shape[1]
    z_dim = template_ccfv3_int_vol.data.shape[2]

    i = 0
    for x in range(x_dim):
        for y in range(y_dim):
            for z in range(z_dim):
                template_ccfv3_int_vol.data[x,y,z] = i
                i+=1

    ###write output file
    template_ccfv3_int_vol.data = template_ccfv3_int_vol.data.astype(int)
    template_ccfv3_int_vol.writeFile()
        
    ##############################################################################
    ###run mincresample -nearest_neighbour to upsample and create a "label file" that converts 100um voxels --> 200 um voxels.
    ##load in this label file and write out source/target label dictionaries based on masks
    ##############################################################################

    infile="derivatives/convert_100um_to_200um_mnc_labels/average_template_200_ccfv3_int_labels.mnc"
    like="preprocessed/templates/average_template_100_ccfv3.mnc"
    output="derivatives/convert_100um_to_200um_mnc_labels/average_template_200_int_upsampled_100_nn.mnc"
    os.system("mincresample -int " + infile + " -clobber -like " + like + " -nearest_neighbour " + output)
    ###load upsampled space 

    upsample_vol = volumeFromFile(output, dtype="int")

    with open('derivatives/source_target_masks_local/source_target_masks_local_100um.pkl', 'rb') as f:
        source_target_masks=pickle.load(f)
    
    source_mask_local=source_target_masks["source_mask_local"]
    target_mask_local=source_target_masks["target_mask_local"]

    ##come up with "key" arrays to convert source/target mask indices --> corresponding integer representing a 200um voxel
    target_mask_key_200um = np.zeros(target_mask_local.shape[0]) 
    for i in range(target_mask_local.shape[0]):
        x=target_mask_local[i,0]
        y=target_mask_local[i,1]
        z=target_mask_local[i,2]
        target_mask_key_200um[i] = upsample_vol.data[z,y,x]

    ##generate key array for source mask (right hemisphere)
    source_mask_key_200um = np.zeros(source_mask_local.shape[0])
    for i in range(source_mask_local.shape[0]):
        x=source_mask_local[i,0]
        y=source_mask_local[i,1]
        z=source_mask_local[i,2]
        source_mask_key_200um[i] = upsample_vol.data[z,y,x]

    with open('derivatives/downsampled_connectome/source_target_mask_keys_100um_to_200um.pkl','wb') as f:
            output_dict={"source_mask_key_200um":source_mask_key_200um, "target_mask_key_200um":target_mask_key_200um}
            pickle.dump(output_dict, f)
