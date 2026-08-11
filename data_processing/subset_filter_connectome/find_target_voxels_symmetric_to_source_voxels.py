"""
Vikram Nathan, 03/25/2026

Goal: Write out a subset of the "integer labels" (denoting 200um voxels for downsampling;
#see label_downsampled_pir_space.py and regionalize_downsampled_pir_space.py )
##with voxels that are either:
# 1. Within "source voxels" in the Allen Institute's mouse_connectivity_models (right hemisphere only)
##2. A subset of the "target voxels" (RH+LH voxels) that are direct reflections of the "source voxels" in (1)

###The end goal is a ~30k x ~60k matrix with EXACTLY twice as many target voxels as source voxels, which
###allows for easy extraction/"stacking" of ipsilateral/contralateral predictions for SIR modelling purposes. 

"""
import numpy as np
import pandas as pd
from pyminc.volumes.factory import *
import pyminc.volumes.factory as pyminc
import subprocess
import sys
import pickle

label_file_200um_dir="../derivatives/convert_100um_to_200um_mnc_labels/"
subset_label_dir="../derivatives/SIR_inputs_200um/subsetted_masks/"
downsampled_conn_dir="../derivatives/downsampled_connectome/"

if __name__ == "__main__":

    ##1. load 100 --> 200um label file and source voxels
    ###repeat for STR and HIP major division voxels

    prefixes = ['','str_', 'hip_']
    int_label_file_paths = [label_file_200um_dir + "average_template_200_ccfv3_int_labels.mnc",
                            subset_label_dir + "str_major_div_int_200um_labels.mnc",
                            subset_label_dir + "hip_major_div_int_200um_labels.mnc"]

    int_label_file_path_full = label_file_200um_dir + "average_template_200_ccfv3_int_labels.mnc"
    int_label_file_full = pyminc.volumeFromFile(int_label_file_path_full, dtype='int')        
    for i in range(len(prefixes)):
        prefix = prefixes[i]
        int_label_file_path = int_label_file_paths[i]     
        int_voxel_label_file_200um = pyminc.volumeFromFile(int_label_file_path, dtype='int')
        #print(sum(np.array(int_voxel_label_file_200um.data).flatten() > 0))
        #print(len(np.unique(np.array(int_voxel_label_file_200um.data).flatten())))

        ###ensure we have no repeated indices in the label file regions
        #int_voxel_label_file_200um.data[int_voxel_label_file_200um.data > 0] = int_label_file_full.data[int_voxel_label_file_200um.data > 0]
        source_voxel_labels = np.array(pd.read_pickle(downsampled_conn_dir + "voxel_model_200um_row_labels.pkl")).astype(int).flatten()

        ##2. write out source and target voxel volumes
        label_file_voxels_in_source = np.isin(int_voxel_label_file_200um.data, source_voxel_labels)
        source_vox_vol_path=label_file_200um_dir+ prefix + "source_voxels_200um_downsampled_conn_PIR_CCFv3.mnc"
        source_output_vol = pyminc.volumeLikeFile(int_label_file_path, source_vox_vol_path)
        source_output_vol.data = np.zeros(source_output_vol.data.shape)
        source_output_vol.data[label_file_voxels_in_source] = int_voxel_label_file_200um.data[label_file_voxels_in_source]

        source_output_vol.writeFile()
        source_output_vol.closeVolume()

        target_voxel_labels = np.array(pd.read_pickle(downsampled_conn_dir + "voxel_model_200um_column_labels.pkl")).astype(int).flatten()

        ##repeat for target
        label_file_voxels_in_target = np.isin(int_voxel_label_file_200um.data, target_voxel_labels)
        target_vox_vol_path=label_file_200um_dir+ prefix + "target_voxels_200um_downsampled_conn_PIR_CCFv3.mnc"
        target_output_vol = pyminc.volumeLikeFile(int_label_file_path, target_vox_vol_path)
        target_output_vol.data = np.zeros(target_output_vol.data.shape)
        target_output_vol.data[label_file_voxels_in_target] = int_voxel_label_file_200um.data[label_file_voxels_in_target]

        target_output_vol.writeFile()
        target_output_vol.closeVolume()

        ###3. convert source voxel volume to RAS/MICe coordinates (for reflection across midline)
        source_vox_vol_path_transformed=label_file_200um_dir+ prefix + "source_voxels_200um_downsampled_conn_RAS_MICe.mnc"
        subprocess.run("python transform_space_gabe.py --clobber -v \"RAS\" -w \"MICe\" " + source_vox_vol_path + " " + source_vox_vol_path_transformed, shell=True)

        ##4. reflect source voxel volume across midline using mincreshape -xdimension; reshape to PIR/CCFv3 coordinates
        source_vox_vol_path_transformed_flipped=label_file_200um_dir+ prefix + "source_voxels_200um_downsampled_conn_RAS_MICe_flipped.mnc"
        
        subprocess.run("param2xfm -scales -1 1 1 flip_x.xfm", shell=True)
        subprocess.run("mincresample -use_input_sampling -int -clobber -transformation flip_x.xfm "+ source_vox_vol_path_transformed + " " + source_vox_vol_path_transformed_flipped, shell=True)

        source_vox_vol_path_untransformed_flipped=label_file_200um_dir+ prefix + "source_voxels_200um_downsampled_conn_PIR_CCFv3_flipped.mnc"
        subprocess.run("python transform_space_gabe.py --clobber -v \"PIR\" -w \"CCFv3\" " + source_vox_vol_path_transformed_flipped + " " + source_vox_vol_path_untransformed_flipped, shell=True)

        ##5. figure out which target voxels align with reflected source voxels; make sure this equals the number of source voxels!! 
        flipped_source_file_vol_200um = pyminc.volumeFromFile(source_vox_vol_path_untransformed_flipped, dtype='int')
        target_vol = pyminc.volumeFromFile(target_vox_vol_path, dtype='int')

        flattened_source_voxels_flipped = np.array(flipped_source_file_vol_200um.data).flatten()
        flattened_target_voxels = np.array(target_vol.data).flatten()
        nonzero_flipped_source_voxel_indices=np.where(flattened_source_voxels_flipped != 0)[0]
        nonzero_target_voxel_indices = np.where(flattened_target_voxels !=0)[0]

        flipped_source_voxel_indices_in_target=np.intersect1d(nonzero_flipped_source_voxel_indices, nonzero_target_voxel_indices) 
        #print(len(source_voxel_labels))
        #print(len(flipped_source_voxel_indices_in_target))

        ##6: ENSURE THE ORDERING OF THE TARGET VOXELS IS THE SAME AS THE SOURCE VOXELS; write the ORDERED voxel lists as .pkl files!!! 
        source_voxels_filt = flattened_source_voxels_flipped[flipped_source_voxel_indices_in_target]
        #print(len(source_voxels_filt))
        #print(len(np.unique(source_voxels_filt)))
        target_voxels_filt = flattened_target_voxels[flipped_source_voxel_indices_in_target]

        ###---> concatenate source voxels + LH voxels 
        final_200um_voxel_labels_source_target = np.array([source_voxels_filt, target_voxels_filt]).flatten()

        with open(downsampled_conn_dir + prefix + "source_target_indices_filt_overlap.pkl", 'wb') as f:
            output_dict={"source_indices":source_voxels_filt, "target_indices":final_200um_voxel_labels_source_target}
            #print(np.all(output_dict['target_indices'][:len(source_voxels_filt)] == output_dict['source_indices']))
            pickle.dump(output_dict, f)

        ###double-check that lengths are symmetric
        #print(len(source_voxels_filt))
        #print(len(final_200um_voxel_labels_source_target))
        #print(len(np.unique(source_voxels_filt)))
        #print(len(np.unique(final_200um_voxel_labels_source_target)))
        
        ###write .mnc file to visualize distribution of source/target voxels
        source_target_outfile= prefix + "target_indices_filt_overlap_source_200um.mnc"
        output_path = downsampled_conn_dir + source_target_outfile
        output_vol = pyminc.volumeLikeFile(int_label_file_path, output_path)
        output_vol.data = np.zeros(output_vol.data.shape)
        #print(final_200um_voxel_labels_source_target)
        source_target_voxels_indices_int=np.isin(int_voxel_label_file_200um.data, final_200um_voxel_labels_source_target)
        #print(source_target_voxels_indices_int)
        output_vol.data[source_target_voxels_indices_int] = int_voxel_label_file_200um.data[source_target_voxels_indices_int]

        ##Write output volume
        output_vol.writeFile()
        output_vol.closeVolume()
