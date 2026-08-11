# SIR_mouse_voxel
The following repository is a voxelwise implementation of https://github.com/vik16nathan/SIR_mouse_regional_tuned for predicting MRI atrophy data in mice, following the injection of aSyn into CP and HIP epicentres. Notably, we use a QC'd, voxelwise structural connectome (see Nathan et al., 2026, _Imaging Neuroscience_) and restrict our analyses to capture local aSyn propagation within the CP and HIP voxels, allowing us to tune the SIR model to capture how differences in the simulated spreading patterns predict our empirical data. 

Minor algorithmic changes from https://github.com/srahayel/SIR_mouse to the current implementation of AgentBasedModel.py:

* Since the voxels are isotropic, I removed the inverse dependency on region size when calculating self.trans_rate, the parameter controlling the rate of conversion of S --> I agents. I instead fixed this to be an arbitrarily small constant (0.005; see line 60 of AgentBasedModel.py within SIR_mouse_voxel). 
* To save space/time, I got rid of the model's internal update of s_edge_history/i_edge_history.

