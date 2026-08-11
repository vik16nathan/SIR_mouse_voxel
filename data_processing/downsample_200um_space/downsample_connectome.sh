#!/bin/bash
source .venv/bin/activate
python label_downsampled_pir_space.py
deactivate

source mouse_connectivity_models/.venv/bin/activate
python regionalize_downsampled_pir_space.py
