# SIR_mouse_voxel
The following repository is a voxelwise implementation of ... for predicting atrophy data using the rebuilt connectome from Nathan et al., 2026, _Imaging Neuroscience_.

Minor algorithmic changes from https://github.com/srahayel/SIR_mouse to the current implementation of AgentBasedModel.py:

* Since the voxels are isotropic, I removed the inverse dependency on region size when calculating self.trans_rate, the parameter controlling the rate of conversion of S --> I agents. I instead fixed this to be an arbitrarily small constant (0.005; see line 60 of AgentBasedModel.py within SIR_mouse_voxel). 
* To save space/time, I got rid of the model's internal update of s_edge_history/i_edge_history.

