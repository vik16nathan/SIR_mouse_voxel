import numpy as np
import pandas as pd
import pickle 
from pyminc.volumes.factory import *
import pyminc.volumes.factory as pyminc
from pathlib import Path

atrophy_dir="../derivatives/SIR_inputs_200um/symmetric_source_target_masked/atrophy/"
Snca_dir="../derivatives/SIR_inputs_200um/symmetric_source_target_masked/GE/"
subset_index_dir="../derivatives/SIR_inputs_200um/subsetted_masks/"
output_dir="../derivatives/SIR_inputs_200um/subsetted_GE_atrophy/"

if __name__ == "__main__":

    
    ##0. load full PIR files for all inputs (AFTER filtering by symmetric source/target)
    atrophy_pir_files_to_subset=[atrophy_dir+"M83PBS_vs_M83HuPff_t_stats_linear_allen200_PIR_CCFv3_target_filt.mnc",
                                atrophy_dir+"M83PBS_vs_M83MsPff_t_stats_linear_allen200_PIR_CCFv3_target_filt.mnc",
                                atrophy_dir + "M83PBS_vs_M83HuPff_hipp_t_stats_linear_allen200_PIR_CCFv3_target_filt.mnc"]
    ge_pir_files_to_subset=[Snca_dir+"Snca_coronal_target_filt.mnc"]

    for file_set in [atrophy_pir_files_to_subset, ge_pir_files_to_subset]:
        for file in file_set:
            input_mnc_vol=pyminc.volumeFromFile(file)
            basename=Path(file).stem
            for major_div in ["str", "hip"]:
                
                ##1. load subset minc file (AS INTEGER)
                major_div_subset_mnc_file = pyminc.volumeFromFile(subset_index_dir+major_div+"_major_div_int_200um_labels.mnc", dtype='int')

                ##2. Create output volume
                outfile=major_div+"_"+basename+".mnc"
                output_path = output_dir + outfile
                output_vol = pyminc.volumeLikeFile(file, output_path)
                output_vol.data = np.zeros(output_vol.data.shape)

                ##3. subset using binary mask of subset_file > 0
                output_vol.data[major_div_subset_mnc_file.data > 0] = input_mnc_vol.data[major_div_subset_mnc_file.data > 0]

                ##Write output volume
                output_vol.writeFile()
                output_vol.closeVolume()