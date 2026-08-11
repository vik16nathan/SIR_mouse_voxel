import numpy as np
import pandas as pd
import os
from pathlib import Path

os.chdir(".")

###dictionary to convert from epicentre --> number within .pkl file sources
epi_num_dict = {'CP': 111045, 'DG': 110656}
num_reps_baseline = 40

clearance_gene_dir_orig = "../../derivatives/SIR_inputs_200um/symmetric_source_target_masked/GE_pkl/"
clearance_list_path = "../../../MouseHumanTranscriptomicSimilarity/AMBA/data/gene_names_coronal_mask_0.8_allgene_filt.csv"
sir_commands_dir = "../sir_command_files/baselines/"
results_dir = "../sir_result_csvs/baselines/"
atrophy_maps_dir = "../../derivatives/SIR_inputs_200um/symmetric_source_target_masked/atrophy_pkl/"
connectome_params_dir = "../../derivatives/SIR_inputs_200um/symmetric_source_target_masked/params/"

cp_hipp_atrophy_maps_comparison = {
    'CP': [
        "str/M83PBS_vs_M83MsPff_t_stats_linear_allen200_PIR_CCFv3_target_filt.pkl"
    ],
    'DG': [
        "hip/M83PBS_vs_M83HuPff_hipp_t_stats_linear_allen200_PIR_CCFv3_target_filt.pkl"
    ]
}

cp_hipp_atrophy_maps_comparison_full = {
    'CP': [
        "M83PBS_vs_M83MsPff_t_stats_linear_allen200_PIR_CCFv3_target_filt.pkl"
    ],
    'DG': [
        "M83PBS_vs_M83HuPff_hipp_t_stats_linear_allen200_PIR_CCFv3_target_filt.pkl"
    ]
}

###epsilon for convergence
eps = 1e-5

###list of connectome parameters
connectome_params_dict_rvm = {
    'CP': {'antero': "source_target_indices_filt_overlap_str_top30.pkl",
           'retro':  "source_target_indices_filt_overlap_str_top30_retro.pkl"},
    'DG': {'antero': "source_target_indices_filt_overlap_hip_top30.pkl",
           'retro':  "source_target_indices_filt_overlap_hip_top30_retro.pkl"}
}

###whole-brain
connectome_params_dict_rvm_full = {
    'CP': {'antero': "source_target_indices_filt_overlap_top30.pkl",
           'retro':  "source_target_indices_filt_overlap_top30_retro.pkl"},
    'DG': {'antero': "source_target_indices_filt_overlap_top30.pkl",
           'retro':  "source_target_indices_filt_overlap_top30_retro.pkl"}
}

if __name__ == "__main__":

    conn_suff = "voxel_model_top30"  # FIX 1: moved outside inner loops (was mis-indented as a block opener)

    for epicentre in ['CP', 'DG']:
        epi = str(epi_num_dict[epicentre])

        for direction in ['antero', 'retro']:
            for atrophy_map in cp_hipp_atrophy_maps_comparison[epicentre]:

                atrophy_map_path = atrophy_maps_dir + atrophy_map
                atrophy_map_stem = Path(atrophy_map).stem

                atrophy_results_dir = results_dir + "atrophy/" + atrophy_map_stem + "/" + direction + "/"
                os.makedirs(sir_commands_dir, exist_ok=True)
                os.makedirs(atrophy_results_dir, exist_ok=True)

                params_file = connectome_params_dict_rvm[epicentre][direction]  # FIX 2: epicentre not epi
                params_file_path = connectome_params_dir + params_file

                command_filename = (
                    sir_commands_dir
                    + f"{atrophy_map_stem}_{direction}_{conn_suff}_baseline.txt"
                )

                with open(command_filename, "w") as f:

                    for rep in range(num_reps_baseline):

                        cmd = (
                            f"python3 ../algorithm/abm_optuna_general.py "
                            f"-a True "
                            f"-g {clearance_gene_dir_orig} "
                            f"-m {atrophy_map_path} "
                            f"-p {params_file_path} "
                            f"-r {'True' if direction == 'retro' else 'False'} "
                            f"-t 1000 "
                            f"-d 0.1 "
                            f"-x {atrophy_map_stem}_{rep}_{conn_suff} "
                            f"-e {eps} "
                            f"-S {epi} "
                            f">> {atrophy_results_dir}{atrophy_map_stem}_{conn_suff}_{direction}.csv"
                            f"\n"
                        )

                        f.write(cmd)
    

    #########REPEAT BASELINES FOR IPSILATERAL ONLY##########
    conn_suff = "voxel_model_top30_ipsi"  # FIX 1: moved outside inner loops (was mis-indented as a block opener)

    for epicentre in ['CP', 'DG']:
        epi = str(epi_num_dict[epicentre])

        for direction in ['antero', 'retro']:
            for atrophy_map in cp_hipp_atrophy_maps_comparison[epicentre]:

                atrophy_map_path = atrophy_maps_dir + atrophy_map
                atrophy_map_stem = Path(atrophy_map).stem

                atrophy_results_dir = results_dir + "atrophy/" + atrophy_map_stem + "/" + direction + "/"
                os.makedirs(sir_commands_dir, exist_ok=True)
                os.makedirs(atrophy_results_dir, exist_ok=True)

                params_file = connectome_params_dict_rvm[epicentre][direction]  # FIX 2: epicentre not epi
                params_file_path = connectome_params_dir + params_file

                command_filename = (
                    sir_commands_dir
                    + f"{atrophy_map_stem}_{direction}_{conn_suff}_baseline.txt"
                )

                with open(command_filename, "w") as f:

                    for rep in range(num_reps_baseline):

                        cmd = (
                            f"python3 ../algorithm/abm_optuna_general.py "
                            f"-i True "
                            f"-a True "
                            f"-g {clearance_gene_dir_orig} "
                            f"-m {atrophy_map_path} "
                            f"-p {params_file_path} "
                            f"-r {'True' if direction == 'retro' else 'False'} "
                            f"-t 1000 "
                            f"-d 0.1 "
                            f"-x {atrophy_map_stem}_{rep}_{conn_suff} "
                            f"-e {eps} "
                            f"-S {epi} "
                            f">> {atrophy_results_dir}{atrophy_map_stem}_{conn_suff}_{direction}.csv"
                            f"\n"
                        )

                        f.write(cmd)
    
    #########REPEAT BASELINES FOR WHOLE-BRAIN PREDICTION####
    #conn_suff = "voxel_model_top30_wholebrain" 
    #for epicentre in ['CP', 'DG']:
    #    epi = str(epi_num_dict[epicentre])

    #    for direction in ['antero', 'retro']:
    #        for atrophy_map in cp_hipp_atrophy_maps_comparison_full[epicentre]:

    #            atrophy_map_path = atrophy_maps_dir + atrophy_map
    #            atrophy_map_stem = Path(atrophy_map).stem

    #            atrophy_results_dir = results_dir + "atrophy/" + atrophy_map_stem + "/" + direction + "/"
    #            os.makedirs(sir_commands_dir, exist_ok=True)
    #            os.makedirs(atrophy_results_dir, exist_ok=True)

    #            params_file = connectome_params_dict_rvm_full[epicentre][direction]  # FIX 2: epicentre not epi
    #            params_file_path = connectome_params_dir + params_file

    #            command_filename = (
    #                sir_commands_dir
    #                + f"{atrophy_map_stem}_{direction}_{conn_suff}_baseline_wholebrain.txt"
    #            )

    #            with open(command_filename, "w") as f:

    #                for rep in range(num_reps_baseline):

    #                    cmd = (
    #                        f"python3 ../algorithm/abm_optuna_general.py "
    #                        f"-a True "
    #                        f"-g {clearance_gene_dir_orig} "
    #                        f"-m {atrophy_map_path} "
    #                        f"-p {params_file_path} "
    #                        f"-r {'True' if direction == 'retro' else 'False'} "
    #                        f"-t 200 " ###reduce number of timesteps for infection to spread without reducing correlation!!! 
    #                        f"-y 20 "
    #                        f"-d 0.1 "
    #                        f"-x {atrophy_map_stem}_{rep}_{conn_suff} "
    #                        f"-e {eps} "
    #                        f"-S {epi} "
    #                        f">> {atrophy_results_dir}{atrophy_map_stem}_{conn_suff}_{direction}.csv"
    #                        f"\n"
    #                    )

    #                    f.write(cmd)
    
    #########REPEAT BASELINES FOR WHOLE-BRAIN; IPSILATERAL ########
    #conn_suff = "voxel_model_top30_wholebrain_ipsi" 
    #for epicentre in ['CP', 'DG']:
    #    epi = str(epi_num_dict[epicentre])

    #    for direction in ['antero', 'retro']:
    #        for atrophy_map in cp_hipp_atrophy_maps_comparison_full[epicentre]:

    #            atrophy_map_path = atrophy_maps_dir + atrophy_map
    #            atrophy_map_stem = Path(atrophy_map).stem

    #            atrophy_results_dir = results_dir + "atrophy/" + atrophy_map_stem + "/" + direction + "/"
    #            os.makedirs(sir_commands_dir, exist_ok=True)
    #            os.makedirs(atrophy_results_dir, exist_ok=True)

    #            params_file = connectome_params_dict_rvm_full[epicentre][direction]  # FIX 2: epicentre not epi
    #            params_file_path = connectome_params_dir + params_file

    #            command_filename = (
    #                sir_commands_dir
    #                + f"{atrophy_map_stem}_{direction}_{conn_suff}_baseline.txt"
    #            )

    #            with open(command_filename, "w") as f:

    #                for rep in range(num_reps_baseline):

    #                    cmd = (
    #                        f"python3 ../algorithm/abm_optuna_general.py "
    #                        f"-a True "
    #                        f"-g {clearance_gene_dir_orig} "
    #                        f"-m {atrophy_map_path} "
    #                        f"-p {params_file_path} "
    #                        f"-r {'True' if direction == 'retro' else 'False'} "
    #                        f"-t 100 " ###reduce number of timesteps for infection to spread without reducing correlation!!! 
    #                        f"-y 10 "
    #                        f"-d 0.1 "
    #                        f"-x {atrophy_map_stem}_{rep}_{conn_suff} "
    #                        f"-e {eps} "
    #                        f"-S {epi} "
    #                        f">> {atrophy_results_dir}{atrophy_map_stem}_{conn_suff}_{direction}.csv"
    #                        f"\n"
    #                    )

    #                    f.write(cmd)

    #################CLEARANCE GENES, IPSILATERAL#########################
    sir_commands_dir = "../sir_command_files/clearance/"
    results_dir = "../sir_result_csvs/clearance/"

    clearance_list = pd.read_csv(clearance_list_path)['Gene']
    conn_suff = "voxel_model_top30"  # FIX 1: moved outside inner loops (was mis-indented as a block opener)

    for epicentre in ['CP', 'DG']:
        epi = str(epi_num_dict[epicentre])

        ###set spreading directions based on highest baseline prediction accuracy
        if epicentre == 'CP':
            direction = 'retro'
            append = "str/"
        if epicentre == 'DG':
            direction = 'antero'
            append = "hip/"
        
        clearance_gene_dir = clearance_gene_dir_orig + append ###note: epicentre-specific spreading within SINGLE MAJOR DIVISION

        for atrophy_map in cp_hipp_atrophy_maps_comparison[epicentre]:

            atrophy_map_path = atrophy_maps_dir + atrophy_map
            atrophy_map_stem = Path(atrophy_map).stem

            atrophy_results_dir = results_dir + "atrophy/" + atrophy_map_stem + "/" + direction + "/"
            os.makedirs(sir_commands_dir, exist_ok=True)
            os.makedirs(atrophy_results_dir, exist_ok=True)

            params_file = connectome_params_dict_rvm[epicentre][direction]  
            params_file_path = connectome_params_dir + params_file

            command_filename = (
                sir_commands_dir
                + f"{atrophy_map_stem}_{direction}_{conn_suff}_clearance_ipsi.txt"
            )

            with open(command_filename, "w") as f:
                for gene in clearance_list:
                    cmd = (
                        f"python3 ../algorithm/abm_optuna_general.py "
                        f"-i True "
                        f"-a True "
                        f"-g {clearance_gene_dir} "
                        f"-c {gene} "
                        f"-m {atrophy_map_path} "
                        f"-p {params_file_path} "
                        f"-r {'True' if direction == 'retro' else 'False'} "
                        f"-t 1000 "
                        f"-d 0.1 "
                        f"-e {eps} "
                        f"-S {epi} "
                        f">> {atrophy_results_dir}{gene}.csv"
                        f"\n"
                    )

                    f.write(cmd)
