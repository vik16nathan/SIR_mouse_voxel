#!/bin/bash

source ../.venv/bin/activate


#############BASELINES, BILATERAL#########################
###ABM, CP epicentre, retrograde (HIGHER ACCURACY)####
python ../algorithm/abm_clearance_genes.py \
  -t 1000 \
  -p "../../derivatives/SIR_inputs_200um/symmetric_source_target_masked/params/source_target_indices_filt_overlap_str_top30_retro.pkl" \
  -o "../simulations/" \
  -r "True" \
  -d 0.1 \
  -S 111045 \
  -v 0.0017003872488696233 \
  -s 0.9971421696522822 \
  -i 57.90868321082842 \
  -k1 0.9345289210280755 \
  -k2 0.5583323404362099

###ABM, CP epicentre, anterograde (LOWER ACCURACY)####
python ../algorithm/abm_clearance_genes.py \
  -t 1000 \
  -p "../../derivatives/SIR_inputs_200um/symmetric_source_target_masked/params/source_target_indices_filt_overlap_str_top30.pkl" \
  -o "../simulations/" \
  -r "False" \
  -d 0.1 \
  -S 111045 \
  -v 0.7130185134308772 \
  -s 0.8800327881355907 \
  -i 56.16471050842456 \
  -k1 0.8091228425272163 \
  -k2 0.9871793893663503


####ABM, DG epicentre, anterograde (HIGHER ACCURACY)#####
####choose an optimal parameter setting from sir_result_csvs/ ###
python3 ../algorithm/abm_clearance_genes.py \
  -t 1000 \
  -p "../../derivatives/SIR_inputs_200um/symmetric_source_target_masked/params/source_target_indices_filt_overlap_hip_top30.pkl" \
  -o "../simulations/" \
  -r "False" \
  -d 0.1 \
  -S 110656 \
  -v 0.97071891251077 \
  -s 0.6670313303245704 \
  -i 84.65707027363796 \
  -k1 0.922920498771323 \
  -k2 0.42191190675374046

####ABM, DG epicentre, retrograde (LOWER ACCURACY)#####
####choose an optimal parameter setting from sir_result_csvs/ ###
python3 ../algorithm/abm_clearance_genes.py \
  -t 1000 \
  -p "../../derivatives/SIR_inputs_200um/symmetric_source_target_masked/params/source_target_indices_filt_overlap_hip_top30_retro.pkl" \
  -o "../simulations/" \
  -r "True" \
  -d 0.1 \
  -S 110656 \
  -v 0.011944797633727307 \
  -s 0.9832548416531878 \
  -i 1.1848216069603854 \
  -k1 0.08188919228370305 \
  -k2 0.8618764079000241


#############BASELINES, IPSILATERAL#########################

###CP, retro#######
python ../algorithm/abm_clearance_genes.py \
  -l "True" \
  -t 1000 \
  -p "../../derivatives/SIR_inputs_200um/symmetric_source_target_masked/params/source_target_indices_filt_overlap_str_top30_retro.pkl" \
  -o "../simulations/" \
  -r "True" \
  -d 0.1 \
  -S 111045 \
  -v 0.37707557027361915 \
  -s 0.7819432344877069 \
  -i 55.92661754954737 \
  -k1 0.09022400830018472 \
  -k2 0.6313247056289887

###DG, antero######
python3 ../algorithm/abm_clearance_genes.py \
  -l "True" \
  -t 1000 \
  -p "../../derivatives/SIR_inputs_200um/symmetric_source_target_masked/params/source_target_indices_filt_overlap_hip_top30.pkl" \
  -o "../simulations/" \
  -r "False" \
  -d 0.1 \
  -S 110656 \
  -v 0.008746234591760752 \
  -s 0.9666841451733642 \
  -i 99.88133696594801 \
  -k1 0.6870137681953692 \
  -k2 0.012801212760166067


#############CLEARANCE, IPSILATERAL#########################
###CP, retro, top gene#######
python ../algorithm/abm_clearance_genes.py \
  -c "Sh3bgr" \
  -g "../../derivatives/SIR_inputs_200um/symmetric_source_target_masked/GE_pkl/str/" \
  -l "True" \
  -t 1000 \
  -p "../../derivatives/SIR_inputs_200um/symmetric_source_target_masked/params/source_target_indices_filt_overlap_str_top30_retro.pkl" \
  -o "../simulations/" \
  -r "True" \
  -d 0.1 \
  -S 111045 \
  -v 0.9203497196853531 \
  -s 0.817754260835934 \
  -i 93.28826921803088 \
  -k1 0.10895190647204175 \
  -k2 0.799803880680589

###DG, antero, top gene######
python3 ../algorithm/abm_clearance_genes.py \
  -c "Efna3" \
  -g "../../derivatives/SIR_inputs_200um/symmetric_source_target_masked/GE_pkl/hip/" \
  -l "True" \
  -t 1000 \
  -p "../../derivatives/SIR_inputs_200um/symmetric_source_target_masked/params/source_target_indices_filt_overlap_hip_top30.pkl" \
  -o "../simulations/" \
  -r "False" \
  -d 0.1 \
  -S 110656 \
  -v 0.8968132793840868 \
  -s 0.33268056212600106 \
  -i 71.10174286797589 \
  -k1 0.9887085449902345 \
  -k2 0.0048143225873943105