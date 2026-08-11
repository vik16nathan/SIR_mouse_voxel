"""
Vikram Nathan, 03/21/2026

Converts kNN interpolated GE values from Beauchamp et al., 2022 back into .mnc volumes

note: all data is in RAS orientation with MICe world coords; need to convert back to PIR/CCFv3
using transform_space_gabe.py to be able to use with connectivity data defined for 200um "labels" in SIR modelling


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


working_dir="../"
output_dir="derivatives/gene_expression_kNN/RAS_MICe/"
ge_dir="../MouseHumanTranscriptomicSimilarity/AMBA/data/"



def flat_ge_row_to_vol(ge, template_path, mask_path, output_path):
    """
    Inputs: 
    - ge: a single row of gene expression, with columns representing voxel indices (0-61,314 for coronal)
    - template_path: template file whose format will be copied in pyminc.volumeLikeFile
    - mask_path: binary mask file; "1" values define the voxels that correspond to the indices in ge
    - output_path: path for the output volume


    """
    vol = pyminc.volumeLikeFile(template_path, output_path)
    mask_vol = pyminc.volumeFromFile(mask_path)
    mask_vol_flattened = np.array(mask_vol.data).flatten().astype(int)
    
    x_dim, y_dim, z_dim = vol.data.shape
    total_voxels = x_dim * y_dim * z_dim

    # Build a flat array of zeros, then fill in GE only the indices present in the coverage mask
    flat = np.zeros(total_voxels)
    flat[mask_vol_flattened == 1] = np.array(ge.values).flatten() 

    vol.data[...] = flat.reshape(x_dim, y_dim, z_dim)
    vol.writeFile()
    vol.closeVolume()

def process_gene(args):
    gene, ge_values, template_path, coverage_path, output_dir = args
    output_path = output_dir + gene + "_coronal.mnc"
    flat_ge_row_to_vol(ge_values, template_path, coverage_path, output_path)
    return gene

def convert_gene_to_pir(gene):
    
    subprocess.run("python analysis/transform_space_gabe.py -v \"PIR\" -w \"CCFv3\" -x 1 derivatives/gene_expression_kNN/RAS_MICe/"+gene+"_coronal.mnc " + "derivatives/gene_expression_kNN/CCFv3_PIR/"+gene+"_coronal.mnc", shell=True)

if __name__ == "__main__":
    os.chdir(working_dir)

    ge_coronal = pd.read_csv(ge_dir + "MouseExpressionMatrix_voxel_coronal_maskcoronal_grouped_imputed.csv")
    ge_coronal.set_index('Gene', inplace=True)

    ##write out gene list
    gene_list_coronal = pd.DataFrame(ge_coronal.index).to_csv(ge_dir+"gene_names_coronal_mask_0.8_allgene_filt.csv", index=False)
    print("Wrote gene list")

    template_path = "preprocessed/templates/ccfv3_converted_RAS_MICe/average_template_200.mnc"
    coverage_path = ge_dir + "imaging/coronal_200um_coverage_bin0.8.mnc"

    args_list = [
        (gene, ge_coronal.loc[gene, :], template_path, coverage_path, output_dir)
        for gene in ge_coronal.index
    ]

    with ProcessPoolExecutor() as executor:
        futures = {executor.submit(process_gene, args): args[0] for args in args_list}
        for future in as_completed(futures):
            gene = futures[future]
            try:
                future.result()
                print(f"Done: {gene}")
            except Exception as e:
                print(f"Error processing {gene}: {e}")
    
    ###convert GE to PIR space
    with ProcessPoolExecutor() as executor:
        futures = {executor.submit(convert_gene_to_pir, gene): gene for gene in ge_coronal.index}

        for future in as_completed(futures):
            gene = futures[future]
            try:
                future.result()
                print(f"Done: {gene}")
            except Exception as e:
                print(f"Error processing {gene}: {e}")