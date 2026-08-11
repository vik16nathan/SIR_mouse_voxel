import numpy as np
import pandas as pd
import pickle
import pyminc.volumes.factory as pyminc
from pyminc.volumes.factory import volumeFromFile

downsampled_conn_dir="../derivatives/downsampled_connectome/"
int_200um_voxel_label_file_full_path = downsampled_conn_dir+"target_indices_filt_overlap_source_200um.mnc"

if __name__== "__main__":

    source_target_volume = pyminc.volumeFromFile(int_200um_voxel_label_file_full_path , dtype='int')

    source_target_indices = pd.read_pickle(downsampled_conn_dir  + "source_target_indices_filt_overlap.pkl")
    ###get world coordinates for all source/target indices
    source_vol_indices = source_target_indices['source_indices']
    
    ##get correct ordering of voxel indices

    source_volume_coords = np.array([np.array(np.where(source_target_volume.data == source_vol_indices[i])).flatten() for i in range(len(source_vol_indices))])
    source_world_coords = np.array([source_target_volume.convertVoxelToWorld(source_volume_coords[i,:]) for i in range(source_volume_coords.shape[0])])
    print(source_world_coords)

    ###repeat for target indices
    target_vol_indices=source_target_indices['target_indices']
    target_volume_coords = np.array([np.array(np.where(source_target_volume.data == target_vol_indices[i])).flatten() for i in range(len(target_vol_indices))])
    target_world_coords = np.array([source_target_volume.convertVoxelToWorld(target_volume_coords[i,:]) for i in range(target_volume_coords.shape[0])])
    print(target_world_coords)

    ####calculate pairwise distance b/t each source coordinate and target coordinate
    pairwise_distance = np.array([[np.linalg.norm(source_world_coords[i] - target_world_coords[j]) for j in range(len(target_vol_indices))] for i in range(len(source_vol_indices)) ])
    pairwise_distance_df = pd.DataFrame(pairwise_distance, index=source_vol_indices, columns=target_vol_indices)
    
    ####write out an ordered source x target matrix
    with open(downsampled_conn_dir + "voxel_model_200um_full_pairwise_distance.pkl", 'wb') as f:
        pickle.dump(pairwise_distance_df, f)

    ###AFTER: visualize and double-check that values equal the distance values from the epicentres!! 