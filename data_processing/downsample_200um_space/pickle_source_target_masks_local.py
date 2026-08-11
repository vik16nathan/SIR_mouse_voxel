import pickle
import numpy as np
from mcmodels.core import VoxelModelCache


BASE_DIR="/scratch/vnathan/sir_voxel/"
if __name__ == "__main__":
    ##load source and target masks/coordinates
    cache = VoxelModelCache(manifest_file='analysis/voxel_model_manifest.json')

    _, source_mask, target_mask = cache.get_voxel_connectivity_array()
    source_mask_local=source_mask.coordinates
    target_mask_local=target_mask.coordinates

    with open('derivatives/source_target_masks_local/source_target_masks_local_100um.pkl', 'wb') as f:
        output_dict={"source_mask_local":source_mask_local, "target_mask_local":target_mask_local}
        pickle.dump(output_dict, f)
