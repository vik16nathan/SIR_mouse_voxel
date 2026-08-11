###ABM, CP epicentre####
python ../algorithm/abm_clearance_genes.py \
  -t 1000 \
  -p "../../derivatives/SIR_inputs_200um/symmetric_source_target_masked/params/source_target_indices_filt_overlap_.pkl" \
  -o "../simulations_full/" \
  -r "False" \
  -d 0.1 \
  -S 111045 \
  -v 0.1 \
  -s 0.1 \
  -i 100 \
  -k1 0.5 \
  -k2 0.5



####ABM, DG epicentre#####
python3 ../algorithm/abm_clearance_genes.py \
  -t 1000 \
  -p "../../derivatives/SIR_inputs_200um/symmetric_source_target_masked/params/source_target_indices_filt_overlap_.pkl" \
  -o "../simulations_full/" \
  -r "False" \
  -d 0.1 \
  -S 110656 \
  -v 0.1 \
  -s 0.1 \
  -i 100 \
  -k1 0.5 \
  -k2 0.5