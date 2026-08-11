###taken from read-me.txt in Steph's DBM2Allen directory under SIR derivatives
module load minc-toolkit-v2 ANTs anaconda pyminc
export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREAD=8 # Set the number of threads for ANTs, to speed things up - change this to suit your ma>


##########STEP 1: ALIGN THE DBM TEMPLATES##################################
####for CP injection (Stephanie Tullo) template######
##already ran
#./align_to_Allen.sh /data/chamal/projects/stephanie/long-asyn-PFF-project/derivatives/2level_model_niagara_202202/secondlevel_template0.mnc ../../../long-asyn-PFF-project/derivatives/2level_model_niagara_202202/mask/mask_on_template_corrected.mnc Allen_aligned_ output 100
cp  -r /data/chamal/projects/stephanie/SIR-modelling-project/derivatives/DBM-template-2Allenspace/output ../derivatives/dbm_templates_realigned_to_allen/cp_steph/

####for DG injection (Janice Park) template##########
./align_to_Allen.sh /data/scratch2/janice/alt_inj_site/analysis/final-dbm/template_sharpen_shapeupdate.mnc /data/scratch2/janice/alt_inj_site/analysis/final-dbm/new_mask.mnc Allen_aligned_janice  ../derivatives/dbm_templates_realigned_to_allen/dg_janice/ 100


########STEP 2: USE THE TRANSFORMS TO ALIGN DBM RESULTS --> ALLEN TEMPLATE####

# Define a list of strings
PFFs=("HuPff" "MsPff")

# Iterate over the list
for PFF in "${PFFs[@]}"; do

    mincresample -transformation ../derivatives/dbm_templates_realigned_to_allen/cp_steph/Allen_aligned_.xfm -like ../derivatives/dbm_templates_realigned_to_allen/cp_steph/average_template.mnc \
                -clobber "../derivatives/voxel_atrophy_maps/CP/M83${PFF}_vs_M83PBS_t_stats_linear.mnc" "../derivatives/voxel_atrophy_maps/CP/M83${PFF}_vs_M83PBS_t_stats_linear_allen100.mnc"

    ##and then 100 --> 200 microns

    mincresample -2 -like ../preprocessed/templates/ccfv3_converted_RAS_MICe/average_template_200.mnc \
               -clobber "../derivatives/voxel_atrophy_maps/CP/M83${PFF}_vs_M83PBS_t_stats_linear_allen100.mnc" "../derivatives/voxel_atrophy_maps/CP/M83${PFF}_vs_M83PBS_t_stats_linear_allen200.mnc"
done


###repeat for DG
PFF="HuPff"
mincresample -transformation ../derivatives/dbm_templates_realigned_to_allen/dg_janice/Allen_aligned_janice_.xfm -like ../derivatives/dbm_templates_realigned_to_allen/dg_janice/average_template.mnc \
            "../derivatives/voxel_atrophy_maps/DG/M83${PFF}_vs_M83PBS_t_stats_linear.mnc" "../derivatives/voxel_atrophy_maps/DG/M83${PFF}_vs_M83PBS_t_stats_linear_allen100.mnc"

##and then 100 --> 200 microns

mincresample -2 -like ../preprocessed/templates/ccfv3_converted_RAS_MICe/average_template_200.mnc \
            "../derivatives/voxel_atrophy_maps/DG/M83${PFF}_vs_M83PBS_t_stats_linear_allen100.mnc" "../derivatives/voxel_atrophy_maps/DG/M83${PFF}_vs_M83PBS_t_stats_linear_allen200.mnc"

