library("pacman")
library("RMINC")
pacman::p_load(dplyr, readr, readxl, stringr, ggrepel, tidyr, ggplot2, ggalluvial, tidyr)
setwd(".")
source("./ggslicer/plotting_functions/plotting_functions.R")
source("./ggslicer/plotting_functions/plotting_functions_labels.R")

allen_input_dir <- "../preprocessed/templates/"

###organize region numbers --> names --> major subdivisions
#load the dictionary from label numbers <-- --> region acronyms
aba_region_filepath <- paste0(allen_input_dir, "allen_ccfv3_tree_wang_2020_s2.xlsx")
aba_region_labels <- as.data.frame(read_excel(aba_region_filepath))  
colnames(aba_region_labels) <- aba_region_labels[1,]
aba_region_labels <- aba_region_labels[2:nrow(aba_region_labels),]


##200 microns - NEEDS TO BE IN PIR SPACE
aba_label_filepath <- paste0(allen_input_dir, "AMBA_relabeled_25um_resampled_200um_PIR_CCFv3.mnc")
aba_label_file <- round(mincGetVolume(aba_label_filepath))

###major division dictionary (depth 4/5 in Allen region tree)
major_division_dict <- data.frame(Isocortex=315, OLF=698, HPF=1089, CTXsp=703, STR=477, PAL=803,
                                  Thal=549, Hypothal=1097, Midbrain=313, Pons=771,
                                  Medulla=354, CB=512)

###helper function to find major division from aba_region_labels (table converting from region --> region hierarchy)
find_major_division <- function(label, major_division_dict, aba_region_labels) {
  sample_string <- aba_region_labels[which(aba_region_labels[,"structure ID"] == label),"structure_id_path"]
  region_hierarchy_list <- as.numeric(unlist(strsplit(sample_string, "/")))
  colnames(major_division_dict)[which(major_division_dict %in% region_hierarchy_list)]
}


#######pseudocode:########
###1. load actual label file with anatomical labels for 200um voxels --> aba_label_file (see above)
int_200um_asc_label_file_path <- "../derivatives/convert_100um_to_200um_mnc_labels/average_template_200_ccfv3_int_labels.mnc"

###2. come up with binary masks for which voxels are in str, which voxels are in hip.
major_divisions_for_all_labels <- lapply(aba_label_file, find_major_division, major_division_dict=major_division_dict, aba_region_labels=aba_region_labels)
str_major_div_mask <- ifelse(major_divisions_for_all_labels == "STR", 1, 0)
hip_major_div_mask <- ifelse(major_divisions_for_all_labels == "HPF", 1, 0)

print(sum(str_major_div_mask))
print(sum(hip_major_div_mask))

###3. write out these masks
str_mask_path <- "../derivatives/SIR_inputs_200um/subsetted_masks/str_mask_bin.mnc"
hip_mask_path <- "../derivatives/SIR_inputs_200um/subsetted_masks/hip_mask_bin.mnc"
mincWriteVolume(str_major_div_mask, clobber=TRUE, dtype='int', like=int_200um_asc_label_file_path, str_mask_path)
mincWriteVolume(hip_major_div_mask, clobber=TRUE, dtype='int', like=int_200um_asc_label_file_path, hip_mask_path)

###4. Use mincmath -mult to perform element-wise multiplication of these masks with the file from (1) (ensures integer volumes)
system(paste0("mincmath -clobber -mult -int ", int_200um_asc_label_file_path, " ", str_mask_path, " ../derivatives/SIR_inputs_200um/subsetted_masks/str_major_div_int_200um_labels.mnc"))
system(paste0("mincmath -clobber -mult -int ", int_200um_asc_label_file_path, " ", hip_mask_path, " ../derivatives/SIR_inputs_200um/subsetted_masks/hip_major_div_int_200um_labels.mnc"))