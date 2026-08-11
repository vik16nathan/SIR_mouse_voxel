#!/bin/bash

##### 1: Use mouse_connectivity_models to load connectivity to each epicentre ####
source mouse_connectivity_models/.venv/bin/activate
z_epi_cp=51
y_epi_cp=34
x_epi_cp=79

z_epi_dg=78
y_epi_dg=21
x_epi_dg=79

python find_knox_conn_to_epi.py ${z_epi_cp} ${y_epi_cp} ${x_epi_cp} \
    cp_epi_conn_strength_100um.pkl
     
python find_knox_conn_to_epi.py ${z_epi_dg} ${y_epi_dg} ${x_epi_dg} \
     dg_epi_conn_strength_100um.pkl

#### 2: write output volumes for strength/distance from epicentre using pyminc ###
deactivate
source .venv/bin/activate
python write_knox_conn_strength_distance_vols.py ${z_epi_cp} ${y_epi_cp} ${x_epi_cp} \
     cp_epi_conn_strength_100um.pkl \
     cp_epi_conn_strength_100um_ccfv3_annot.mnc \
     cp_epi_conn_distance_100um_ccfv3_annot.mnc

python write_knox_conn_strength_distance_vols.py ${z_epi_dg} ${y_epi_dg} ${x_epi_dg} \
     dg_epi_conn_strength_100um.pkl \
     dg_epi_conn_strength_100um_ccfv3_annot.mnc \
     dg_epi_conn_distance_100um_ccfv3_annot.mnc


