import sys
import re
import numpy as np
import pandas as pd
import pickle
sys.path.insert(1, './ABM_voxel_model/model/')
from AgentBasedModel import AgentBasedModel
from scipy.stats import zscore, norm
from tqdm import tqdm
import argparse
import optuna
import scipy.stats
from pathlib import Path

V_MAX = 1
SIR_input_dir="../../derivatives/SIR_inputs_200um/symmetric_source_target_masked/params/"

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--ipsilateral", type=str2bool, default=False, dest="ipsilateral",
                        nargs='?', help="Boolean for whether we want ipsilateral only or whole-brain predictions")
    parser.add_argument("-a", "--atrophy", type=str2bool, default=False, dest="atrophy",
                        nargs='?', help="Boolean for whether we want to simulate atrophy or stop at I fraction")
    parser.add_argument("-m", "--map_to_predict", dest="map",
                        nargs='?', help="Pathological map to compare to simulated atrophy")
    parser.add_argument("-x", "--suffix", dest="suffix", default=None,  nargs='?',
                        help="Suffix (if you want a .csv file outputted)")
    parser.add_argument(
        "-g", "--GEdir", default=None, dest="clearance_gene_dir",
         nargs='?', help="Clearance Gene Directory"
    )
    parser.add_argument(
        "-c", "--clearance", default=None, dest="clearance_gene", type=none_or_str,
        help="Specify the gene modulating clearance (omit, or pass 'None', for uniform clearance)"
    )
    
    parser.add_argument(
        "-r", "--retro", default=False, dest="retro",
        type=str2bool, nargs='?', help="Retrograde spreading (default False) "
    )

    parser.add_argument(
        "-t", "--time", default=30000, type=int, nargs='?',
        dest="total_time", help="Total spreading time"
    )
    parser.add_argument(
        "-d", "--delta-t", default=0.1, type=float, nargs='?',
        dest="dt", help="Size of time increment"
    )
    parser.add_argument(
        "-e", "--epsilon", default=1e-5, type=float, nargs='?',
        dest="eps", help="Convergence threshold"
    )
    parser.add_argument(
        "-p", "--params", default="../../derivatives/SIR_inputs/params_nature_yohan.pkl",
        dest="params", nargs='?', help="Connectome (pkl file w/ weights, distance, region size, sources, targets)"
    )
    parser.add_argument(
        "-S", "--seed", type=int, nargs='?',
        dest="seed", help="Simulated seeding site of misfolded alpha-syn"
    )
    parser.add_argument(
        "-y", "--num_trials", type=int, nargs='?', default=50,
        dest="num_trials", help="Number of optuna trials"
    )
    parser.add_argument(
        "-n", dest="study_name",
         nargs='?', help="Study name")
    args = parser.parse_args()
    return args

def str2bool(v):
    if isinstance(v, bool):
       return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


def none_or_str(v):
    if v is None or v.lower() == 'none':
        return None
    return v

def load_params(connectome):
    with open(connectome, "rb") as f:
        params = pickle.load(f)
        weights = params['weights']
        distance = params['distance']
        sources = params['sources']
        targets = params['targets']
        Snca = params['Snca']

    syngene = np.array(Snca).flatten()
    syngene = norm.cdf(zscore(syngene))

    return (
            weights, distance,
            sources, targets, syngene
        )
         
def load_clearance(clearance_gene=None, ipsi=False): 
    if clearance_gene is None or clearance_gene == "None":
        return 0.5
    else:
            
        ge = pd.read_pickle(clearance_gene_dir + clearance_gene + '_coronal_target_filt.pkl')
        if ipsi:
            ge = ge.iloc[:, :np.int64(len(ge.columns)/2)] ###ipsilateral only; set up so first half of indices are ipsi
        
        epr = np.array(ge).flatten()

        ###no need to append anymore; kNN interpolation smooths GE across all target voxels!! 
        #epr = np.append(epr, epr)
        return norm.cdf(zscore(epr))

def simulate_pathology(ipsilateral, atrophy, connectome, v, spread_rate, injection_amount, clearance_gene, dt, eps, k1_atrophy=None, k2_atrophy=None):
    weights, distance, sources, targets, syngene = load_params(connectome)
    clearance_rate = load_clearance(clearance_gene, ipsilateral)

    if ipsilateral:
        targets = targets[:len(sources)]
        weights = weights.iloc[:, :len(sources)]
        distance = distance.iloc[:, :len(sources)]
        syngene = syngene[:len(sources)]

    abm = AgentBasedModel(
        weights=weights, distance=distance,
        sources=sources, targets=targets, dt=1
    )

    abm.set_growth_process(growth_rate=syngene)
    abm.set_clearance_process(clearance_rate=clearance_rate)
    abm.set_spread_process(v=v)
    abm.update_spread_process(spread_scale=spread_rate)

    # Growth process — stage 1
    for t in range(30000):
        prev = np.copy(abm.s_voxel)
        abm.growth_step()
        abm.clearance_step()
        abm.s_spread_step()
        if np.any(abm.s_voxel < 0):
            break
        if np.where(np.abs(prev - abm.s_voxel) / abm.s_voxel > eps, 1, 0).sum() == 0:
            break

    # Growth process — stage 2
    abm.dt = 0.1
    for t in range(500000000):
        prev = np.copy(abm.s_voxel)
        abm.growth_step()
        abm.clearance_step()
        abm.s_spread_step()
        if np.any(abm.s_voxel < 0):
            break
        if np.where(np.abs(prev - abm.s_voxel) / abm.s_voxel > eps, 1, 0).sum() == 0:
            break

    # Growth process — stage 3
    abm.dt = dt
    for t in range(1000000000):
        prev = np.copy(abm.s_voxel)
        abm.growth_step()
        abm.clearance_step()
        abm.s_spread_step()
        if np.any(abm.s_voxel < 0):
            break
        if np.where(np.abs(prev - abm.s_voxel) / abm.s_voxel > eps, 1, 0).sum() == 0:
            break

    # Spread process
    injection_location = np.where(abm.targets == seed)[0]
    abm.injection(seed=injection_location, amount=injection_amount)

    for t in tqdm(range(total_time)):
        try:
            abm.growth_step()
            abm.clearance_step()
            abm.trans_step()
            abm.s_spread_step()
            abm.i_spread_step()
            abm.record_history_voxel()
        except:
            print("Numeric exception (see above); moving onto next Optuna trial")
            break

    infected_voxel_time = abm.i_voxel_history

    if atrophy:
        normal_voxel_time = abm.s_voxel_history
        dt = abm.dt

        total_proteins_voxel_time = normal_voxel_time + infected_voxel_time
        infected_proteins_ratio_voxel_time = infected_voxel_time / total_proteins_voxel_time

        sum_strength = np.sum(weights, axis=1)
        ratio_weights = np.array(weights) / np.repeat(sum_strength[:, np.newaxis], len(targets), axis=1)

        k1 = k1_atrophy
        k2 = k2_atrophy

        if ipsilateral:
            new_ratio_weights = ratio_weights
        else:
            new_ratio_weights = np.zeros((len(targets), len(targets)))
            new_ratio_weights[:len(sources), :] = ratio_weights
            new_ratio_weights[len(sources):, :len(sources)] = ratio_weights[:, len(sources):]
            new_ratio_weights[len(sources):, len(sources):] = ratio_weights[:, :len(sources)]

        ratio_cum = np.dot(new_ratio_weights, (1 - np.exp(-infected_proteins_ratio_voxel_time.T * dt)))
        ratio_cum = np.hstack([np.zeros((len(targets), 1)), ratio_cum[:, :-1]])
        ratio_cum = (k2 * ratio_cum) + (k1 * (1 - np.exp(-infected_proteins_ratio_voxel_time.T * dt)))
        simulated = np.cumsum(ratio_cum, axis=1)
    else:
        simulated = infected_voxel_time.T

    return simulated
    
def find_peak_fit(simulated, map, ipsi):
    input_dir = SIR_input_dir
    pathology = pd.read_pickle(map)

    n_sources = np.int64(len(np.array(pathology).flatten())/2)
    if ipsi:
        pathology = np.array(pathology).flatten()[:n_sources]

    
    sim_pathology_data = pd.DataFrame(simulated)
    results = pathology
    ###ASSUME ROW NAMES ARE ALIGNED: this is ensured when processing inputs!!!###
    empirical_common = np.array(results).flatten()
    simulated_common = sim_pathology_data

    nsteps = simulated_common.shape[1]
    correlations = np.zeros(nsteps)

    for i in range(nsteps):
        correlations[i] = scipy.stats.spearmanr(empirical_common, np.array(simulated_common.loc[:, i]).flatten()).correlation
    
    peak_corr = np.max(correlations)
    peak_timestep = np.argmax(correlations)
    return (peak_corr, peak_timestep)

class Objective:
    def __init__(self, ipsilateral, atrophy, dt, eps, seed, total_time, clearance_gene, map, connectome, clearance_gene_dir=None):
        self.ipsilateral = ipsilateral
        self.atrophy = atrophy
        self.dt = dt
        self.eps = eps
        self.seed = seed
        self.total_time = total_time
        self.clearance_gene = clearance_gene
        self.clearance_gene_dir = clearance_gene_dir
        self.map = map
        self.connectome = connectome

    def __call__(self, trial):
            v = trial.suggest_float("v", 1e-3, V_MAX, log=True)
            spread_rate = trial.suggest_float("spread_rate", 0.0001, 1)
            injection_amount = trial.suggest_float("injection_amount", 1, 100)

            if self.atrophy:
                k1_atrophy = trial.suggest_float("k1_atrophy", 0.001, 1)
                k2_atrophy = trial.suggest_float("k2_atrophy", 0.001, 1)
                simulated = simulate_pathology(self.ipsilateral, self.atrophy, self.connectome, v, spread_rate, injection_amount, clearance_gene, dt, eps, k1_atrophy, k2_atrophy)
            else:
                simulated = simulate_pathology(self.ipsilateral, self.atrophy, self.connectome, v, spread_rate, injection_amount, clearance_gene, dt, eps)

            peak_corr, peak_timestep = find_peak_fit(simulated, map, self.ipsilateral)
            self.peak_timestep = peak_timestep
            trial.set_user_attr("peak_timestep", int(self.peak_timestep))
            return peak_corr

# Callback to store the best Objective instance
best_objective = {"instance": None}

def store_best_objective(study, trial):
    if study.best_trial == trial:
        best_objective["instance"] = trial.user_attrs["peak_timestep"]

###RUN OPTUNA###
if __name__ == "__main__":

    args = parse_arguments()
    ipsilateral = args.ipsilateral
    atrophy = args.atrophy
    retro = args.retro
    dt = args.dt
    eps = args.eps
    seed = args.seed
    total_time = args.total_time
    clearance_gene = args.clearance_gene
    clearance_gene_dir = args.clearance_gene_dir
    map = args.map
    connectome = args.params

    study_name = args.study_name
    suffix = args.suffix
    num_trials = args.num_trials
    
    if study_name is None:
        study = optuna.create_study(direction="maximize")
    else:
        study = optuna.create_study(study_name=study_name, storage="sqlite:///../optuna_dbs/" + re.sub(" ", "_", study_name) + ".sqlite3", direction="maximize")

    
    study.optimize(Objective(ipsilateral, atrophy, dt, eps, seed, total_time, clearance_gene, map, connectome, clearance_gene_dir), n_trials=num_trials, callbacks=[store_best_objective])

    if best_objective["instance"] is not None:
        peak_timestep = best_objective["instance"]

    if suffix is None and study_name is None:
        pass
    else:
        df = study.trials_dataframe()
        Path('../optuna_csvs/' + str(seed) + '/').mkdir(parents=True, exist_ok=True)
        if study_name is None:
            df.to_csv('../optuna_csvs/' + str(seed) + '/ret' + str(retro) + '_' + suffix + '.csv')
        else:
            df.to_csv('../optuna_csvs/' + str(seed) + '/' + study_name + '_ret' + str(retro) + '.csv')

        df.to_csv('../optuna_csvs/' + str(seed) + '/ret' + str(retro) + "_" + str(suffix) + '.csv')

    if atrophy:
        print(clearance_gene, study.best_trial.value, study.best_trial.params['v'], study.best_trial.params['spread_rate'],
                        study.best_trial.params['injection_amount'], study.best_trial.params['k1_atrophy'],
                        study.best_trial.params['k2_atrophy'], peak_timestep, sep=",")
    else: 
        print(clearance_gene, study.best_trial.value, study.best_trial.params['v'], study.best_trial.params['spread_rate'],
                        study.best_trial.params['injection_amount'], peak_timestep, sep=",")