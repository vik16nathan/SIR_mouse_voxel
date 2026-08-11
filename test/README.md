# test

A minimal, self-contained usage example for `algorithm/abm_clearance_genes.py`, using the sample connectome/atrophy data checked into this repository under `SIR_inputs_200um/` (see [Data](#data) below) instead of the full external `derivatives/` tree. No clearance gene is used — this demonstrates the model's default uniform-clearance behaviour (see `SIR_mouse_voxel` background below).

## Background

This repository (`SIR_mouse_voxel`) is a voxelwise implementation of [`SIR_mouse_regional_tuned`](https://github.com/vik16nathan/SIR_mouse_regional_tuned) for predicting MRI atrophy data in mice, following the injection of aSyn into CP and HIP epicentres. We use a QC'd, voxelwise structural connectome (see Nathan et al., 2026, *Imaging Neuroscience*) and restrict our analyses to capture local aSyn propagation within the CP and HIP voxels, allowing us to tune the SIR model to capture how differences in the simulated spreading patterns predict our empirical data.

Minor algorithmic changes from [`SIR_mouse`](https://github.com/srahayel/SIR_mouse) in `algorithm/ABM_voxel_model/model/AgentBasedModel.py`:

* Since the voxels are isotropic, the inverse dependency on region size when calculating `self.trans_rate` (the parameter controlling the rate of conversion of S → I agents) was removed and fixed to an arbitrarily small constant (0.005; see `AgentBasedModel.py` line 60).
* The model's internal update of `s_edge_history`/`i_edge_history` was dropped to save space/time.

## Data

To keep this repository lean, only a small subset of `derivatives/SIR_inputs_200um/symmetric_source_target_masked/` is checked in here, under `SIR_inputs_200um/symmetric_source_target_masked/`:

* `params/source_target_indices_filt_overlap_hip.pkl` — one connectome variant (HIP epicentre, anterograde, no top-N gene filtering) out of the ~28 hip/str variants available in the full `derivatives/` tree. **This is not the parameter-tuned/optimized connectome** used in `batch_run/run_abm_voxel_subset_atrophy.sh` (those use the `_top30` variants) — it's included purely to make this example runnable without the external data tree.
* `atrophy/{hip,str}/` and `atrophy_pkl/{hip,str}/` — the empirical MRI atrophy maps (`.mnc` and `.pkl`), tiny (~2MB total).

Per-gene clearance expression data (`GE/`, `GE_pkl/`, ~8GB combined) is **not** included, since this example doesn't use a clearance gene.

## Running the example

From `algorithm/`:

```bash
python3 abm_clearance_genes.py \
  -t 50 \
  -p "../SIR_inputs_200um/symmetric_source_target_masked/params/source_target_indices_filt_overlap_hip.pkl" \
  -o "../test/" \
  -r "False" \
  -d 0.1 \
  -S 110656 \
  -v 0.1 \
  -s 0.1 \
  -i 50 \
  -x "no_clearance_sample"
```

Notes on the arguments:

* No `-c`/`-g` — clearance defaults to `None`, which `abm_clearance_genes.py` treats as a uniform clearance rate of 0.5 everywhere (the `-c` flag also explicitly accepts the string `"None"` for the same effect).
* `-S 110656` is the HIP epicentre voxel index (`epi_num_dict['DG']` in `batch_run/make_optuna_clearance_commands.py`; the CP epicentre is `111045`).
* `-r "False"` matches the non-retrograde connectome variant checked in here (`source_target_indices_filt_overlap_hip.pkl`, no `_retro` suffix).
* `-t 50` is intentionally small for a quick, lean demo (~5MB output, ~3.5 min runtime). Real runs use `-t 1000` (see `batch_run/run_abm_voxel_subset_atrophy.sh` for calibrated examples) — the growth-process convergence phase (~200s) dominates either way and doesn't scale with `-t`, only the spreading-phase loop does.
* `-v`, `-s`, `-i` are illustrative values, not fitted parameters — see `batch_run/run_abm_voxel_subset_atrophy.sh` for actual Optuna-tuned parameter sets per epicentre/direction.

## Output

```
abm_spread_v.0.1.spread_rate.0.1.dt.0.1.seed.110656.injection_amount.50.0.clearance_gene.None..no_clearance_sample.csv
```

5040 rows (target voxels) × 50 columns (timesteps) — the fraction of infected (`I`) agents per voxel at each timestep. The `clearance_gene.None.` in the filename confirms the run used the default uniform clearance rate.

Actual timing from this run (Python 3.10.13, this repo's `.venv`):

```
Begin protein growth process....
entered
Protein growth time: 201.1159 seconds
Begin protein spreading process...
Inject infectious proteins into 110656...
100%|██████████| 50/50 [00:10<00:00, 4.62it/s]
Protein spreading time: 10.8309 seconds
```
