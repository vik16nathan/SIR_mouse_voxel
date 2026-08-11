from mcmodels.core import VoxelModelCache
from mcmodels.models.voxel import RegionalizedModel
import numpy as np
import pickle

########################CHANGE 04/13#########################
nodes_weights_dir="/scratch/vnathan/copy_of_allen_qc_reproducible_02212026/allen_connectome_qc/mouse_connectivity_models/paper/connectivity/voxel-standard-model/"
BASE_DIR="/scratch/vnathan/sir_voxel/"

if __name__ == "__main__":
        with open(BASE_DIR+"derivatives/downsampled_connectome/source_target_mask_keys_100um_to_200um.pkl", 'rb') as file:
                source_target_keys=pickle.load(file)

        source_mask_key_200um = source_target_keys["source_mask_key_200um"]
        target_mask_key_200um = source_target_keys["target_mask_key_200um"]
        ##load nodes and weights matrices (source voxel x region, region x target_voxel)
        nodes_standard = np.loadtxt(nodes_weights_dir+'nodes_rebuilt.csv.gz', delimiter=',')
        weights_standard = np.loadtxt(nodes_weights_dir+'weights_rebuilt.csv.gz', delimiter=',')   

        ###try without chunking on Trillium (more RAM)
        voxel_model_200um_full = RegionalizedModel(weights_standard, nodes_standard, source_mask_key_200um, target_mask_key_200um, dataframe=True)

        ###load both normalized connection density and strength; see which metric is necessary
        normalized_connection_density = getattr(voxel_model_200um_full, "normalized_connection_density")

        ##load normalized connection strength and(region size normalization is no longer necessary; will divide by rowSums anyways)
        normalized_connection_strength = getattr(voxel_model_200um_full, "normalized_connection_strength")

        with open(BASE_DIR+'derivatives/downsampled_connectome/voxel_model_200um_full_density.pkl','wb') as f:
                pickle.dump(normalized_connection_density, f)

        with open(BASE_DIR+'derivatives/downsampled_connectome/voxel_model_200um_full_strength.pkl','wb') as f:
                pickle.dump(normalized_connection_strength, f)

        ###row-normalize the normalized connection strength so the connectivity in each row sums to 1 (for consistency with previous SIR modelling approach)
        row_norm_conn_strength = np.array(normalized_connection_strength) / np.sum(np.array(normalized_connection_strength), axis=1)[:, np.newaxis]
        row_norm_conn_dens = np.array(normalized_connection_density) / np.sum(np.array(normalized_connection_density), axis=1)[:, np.newaxis]

        with open(BASE_DIR+'derivatives/downsampled_connectome/voxel_model_200um_full_density_rownorm.pkl','wb') as f:
                pickle.dump(row_norm_conn_dens, f)

        with open(BASE_DIR+'derivatives/downsampled_connectome/voxel_model_200um_full_strength_rownorm.pkl','wb') as f:
                pickle.dump(row_norm_conn_strength, f)

        print(row_norm_conn_strength == row_norm_conn_dens)

        ####save row and column indices, which are all defined by the integer 200um label file 
        with open(BASE_DIR+'derivatives/downsampled_connectome/voxel_model_200um_row_labels.pkl','wb') as f:
                pickle.dump(normalized_connection_density.index, f)

        with open(BASE_DIR+'derivatives/downsampled_connectome/voxel_model_200um_column_labels.pkl','wb') as f:
                pickle.dump(normalized_connection_density.columns, f)

