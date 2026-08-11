import numpy as np
import pandas as pd
import pickle
from pyminc.volumes.factory import *
import pyminc.volumes.factory as pyminc

connectome_dir="../derivatives/downsampled_connectome/"
subset_index_dir="../derivatives/SIR_inputs_200um/subsetted_masks/"
output_dir="../derivatives/SIR_inputs_200um/subsetted_connectivity/"

if __name__ == "__main__":

    ##0. load full connectome
    full_connectome = pd.read_pickle(connectome_dir+"voxel_model_200um_full.pkl")
    full_connectome.index = np.array(full_connectome.index).astype(int)n
    full_connectome.columns = np.array(full_connectome.columns).astype(int)


    for major_div in ["str", "hip"]:

        ##1. load minc file (AS INTEGER)
        major_div_subset_mnc_file = pyminc.volumeFromFile(subset_index_dir+major_div+"_major_div_int_200um_labels.mnc", dtype='int')

        ##2. get connectome indices by extracting all unique, nonzero values 
        major_div_subset_indices = np.unique(np.array(major_div_subset_mnc_file.data).flatten())

        ##3. subset connectome by indices
        subset_connectome=full_connectome.loc[np.isin(full_connectome.index, major_div_subset_indices), np.isin(full_connectome.columns, major_div_subset_indices)]


        print(major_div, ": ", subset_connectome.shape)

        outfile=major_div+"_subset_connectome.pkl"
        with open(output_dir+outfile, 'wb') as f:
            pickle.dump(subset_connectome, f)
