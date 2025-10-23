from collections import defaultdict
import pandas as pd
import os
import json
import ast

from source.src2_utils.ut0_random_manager import np, random
import config_templates.conf0_model_parameters as conf

def calculate_sm(all_cells_csv:str, synapses_json:str, save_dir:str = None) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Calculates number of synapses for each connection.

    Args:
        all_cells_csv (str): Path to CSV with 'cell_id' and 'desc' (m-type) columns.
        synapses_json (str): Path to JSON with synapse data (pre_id, post_id).
        save_dir (str): Path, where 'sm' directory will be created.

    Returns:
        connections_sm (pd.DataFrame): 2D matrix, where (i, j) cell is the number of synapses between the i-th cell acting as presynaptic and
          j-th cell acting as postsynaptic.
        type_sm (pd.DataFrame): 2D matrix, where (i, j): 'i' is the m-type of the model cells acting as presynaptic, 'j' is the m-type of postsynaptic,
          and values in the matrix represent a tuple with:
            - the mean number of synapses made in the connections between this m-type combination
            - the number of such connections
    """

    # Load cell info
    cell_df = pd.read_csv(all_cells_csv)
    cell_id_to_mtype = dict(zip(cell_df['cell_id'], cell_df['desc']))

    # sort them
    sorted_cell_ids = sorted(cell_id_to_mtype.keys())

    # load synapses
    with open(synapses_json, 'r') as f:
        syn_data = json.load(f)

    #Fetch all synapses as (pre_id, post_id)
    syn_list = []
    for post_id_str, synapses in syn_data.items():
        post_id = int(post_id_str)
        for syn in synapses:
            pre_id = syn['pre_id']
            syn_list.append((pre_id, post_id))

    # count sm for each connection
    connections_sm = pd.DataFrame(0, index=sorted_cell_ids, columns=sorted_cell_ids)
    syn_count_per_pair = defaultdict(int)

    for pre_id, post_id in syn_list:
        syn_count_per_pair[(pre_id, post_id)] += 1
        if pre_id in connections_sm.index and post_id in connections_sm.columns:
            connections_sm.loc[pre_id, post_id] += 1

    # sm for different mtypes
    mtypes = sorted(set(cell_id_to_mtype.values()))
    type_sm = pd.DataFrame(index=mtypes, columns=mtypes, dtype=object)

    # group synapse counts by m-type combinations
    grouped = defaultdict(list)
    for (pre_id, post_id), count in syn_count_per_pair.items():
        pre_m = cell_id_to_mtype.get(pre_id)
        post_m = cell_id_to_mtype.get(post_id)
        if pre_m is not None and post_m is not None:
            grouped[(pre_m, post_m)].append(count)

    # Fill type_sm with (mean synapses per connection, number of connections)
    for pre_m in mtypes:
        for post_m in mtypes:
            counts = grouped.get((pre_m, post_m), [])
            if counts:
                mean_syn = round(sum(counts) / len(counts), 4)
                type_sm.loc[pre_m, post_m] = (mean_syn, len(counts))
            else:
                type_sm.loc[pre_m, post_m] = (0.0, 0)

    if not save_dir is None:

        conn_path = f"{save_dir}/conn_sm.csv"
        connections_sm.to_csv(conn_path)

        types_path = f"{save_dir}/types_sm.csv"
        type_sm_str = type_sm.map(lambda x: f"{x[0]},{x[1]}")
        type_sm_str.to_csv(types_path)

    return connections_sm, type_sm

def calc_b_int(
        full_population_csv:str,
        synapse_json:str,
        savepath:str) -> pd.DataFrame:
    """
    Calculates interbouton density used during 03.Plasticity-reserve pruning. Extends full_population_csv with 'b_int' column.
    This takes two files:
    full_population_csv which has cell_ids and axon_len,
    post02_synapse_json which has synapses after 02 multisyn pruning
    This function calculates interbouton interval for each cell, and adds column b_int within full_population_csv

    Args:
        full_population_csv (str):
        synapse_json (str):
        savepath (str):
    Returns:
        pd.DataFrame: Extended cortical cells DataFrame, with interbouton density.
    """
    # Load population position data
    df = pd.read_csv(full_population_csv)
    if 'cell_id' not in df.columns or 'axon_len' not in df.columns:
        raise KeyError("CSV must contain 'cell_id' and 'axon_len' columns.")

    # Load synapse data
    with open(synapse_json, "r") as f:
        synapses = json.load(f)

    # Count number of synapses made by each pre-synaptic cell
    synapse_count_by_cell = defaultdict(int)
    for post_id_str, syn_list in synapses.items():
        for syn in syn_list:
            pre_id = syn["pre_id"]
            synapse_count_by_cell[pre_id] += 1

    # Compute interbouton interval
    interb_int_list = []
    for _, row in df.iterrows():
        cell_id = row["cell_id"]
        axon_len = row["axon_len"]
        syn_count = synapse_count_by_cell.get(cell_id, 0)

        if syn_count == 0:
            interb_int = float("inf")  # Or use np.nan
        else:
            interb_int = axon_len / syn_count

        interb_int_list.append(interb_int)

    df["b_int"] = interb_int_list

    # Save result
    df.to_csv(savepath, index=False)

    return df

def setup_01_pruning(model_type_sm_csv:str, goal_type_sm_json:str,
                     save_dir:str = None, verbose:bool = False) -> pd.DataFrame:
    """
    Calculates f - survival rate for synapses during the 01 general pruning step. Based on S_m values created by
    calculate_sm.
    Args:
        model_type_sm_csv (str): Path to *csv with sm_values for each me_type:me_type combination prior to any pruning.
        goal_type_sm_json (str): Path to *csv with sm_values for each me_type:me_type combination GIVEN BY BBP.
        save_dir (str): Path where pruning parameters for 01.general pruning are saved.
        verbose (bool): Prints pruning parameters.
    Returns:
        pd.DataFrame: Matrix where each column and row is me_type. Each cell has three parameters (p_goal, p_mod, f)
    """

    model_sm = pd.read_csv(model_type_sm_csv, index_col=0)
    with open(goal_type_sm_json, 'r') as f:
        goal_sm = json.load(f)


    mtypes = model_sm.index
    pruning_df = pd.DataFrame(index=mtypes, columns=mtypes, dtype=object)

    for pre_m in mtypes:
        for post_m in mtypes:
            try:
                # converting string like "(1.23, 42)" into tuple (1.23, 42)
                sm_tuple_str = model_sm.loc[pre_m, post_m]
                sm_model, _ = ast.literal_eval(sm_tuple_str)

                pre_m_short = "_".join(pre_m.split('_')[:2])
                post_m_short = "_".join(post_m.split('_')[:2])
                key = f"{pre_m_short}:{post_m_short}"
                goal_entry = goal_sm.get(key)

                if goal_entry is None:
                    pruning_df.loc[pre_m, post_m] = (-1, -1, -1)
                    continue

                sm_goal = goal_entry["mean_number_of_synapse_per_connection"]

                p_goal = 1 / (sm_goal + conf.NORM_BIAS_01)
                p_goal = min(p_goal, 1.0)

                p_mod = 1 / (sm_model + conf.NORM_BIAS_01)
                p_mod = min(p_mod, 1.0)

                f = (p_goal / (1 - p_goal)) * ((1 - p_mod) / p_mod)
                f = min(f, 1.0)  # clamp to [0,1]

                full_param = (round(p_goal, 4), round(p_mod, 4), round(f, 4))
                if verbose:
                    print(f"{pre_m} → {post_m}: {full_param}")
                pruning_df.loc[pre_m, post_m] = full_param

            except Exception as e:
                print(f"Error for {pre_m}:{post_m} — {e}")
                pruning_df.loc[pre_m, post_m] = (-1, -1, -1)

    if not save_dir is None:
        os.makedirs(save_dir, exist_ok=True)
        prune_01_path = f"{save_dir}/prune_01_params.csv"
        pruning_df_str = pruning_df.map(lambda x: f"{x[0]},{x[1]},{x[2]}")
        pruning_df_str.to_csv(prune_01_path)

    return pruning_df

def apply_general_pruning_01(
    prune_01_params_csv:str,
    init_synapse_json:str,
    output_json:str,
    verbose:bool = True) -> dict:
    """
    Applies general pruning to synapse data using survival factors (f1) per m-type:m-type pair.

    Args:
        prune_01_params_csv (str): CSV with columns ['pre_me_type', 'post_me_type', 'f1']
        init_synapse_json (str): Input JSON with unpruned synapse data
        output_json (str): Output JSON path for pruned synapses
        seed (int): Random seed for reproducibility (optional)
        verbose (bool): Whether to print pruning summary

    Returns:
        dict: Pruned synapse dictionary (post_id -> list of synapses)
    """

    # Load pruning factors DataFrame (cell values are stringified tuples)
    pruning_df = pd.read_csv(prune_01_params_csv, index_col=0)
    pruning_df = pruning_df.map(lambda s: ast.literal_eval(s) if isinstance(s, str) else s)

    # Load initial synapse data
    with open(init_synapse_json, "r") as f:
        synapses = json.load(f)

    # Group synapses by m-type pair
    synapse_groups = defaultdict(list)  # (pre_m, post_m) -> list of (post_id_str, syn)
    for post_id_str, syn_list in synapses.items():
        for syn in syn_list:
            key = (syn["pre_me_type"], syn["post_me_type"])
            synapse_groups[key].append((post_id_str, syn))

    # Apply pruning based on survival factor f from pruning_df
    pruned_synapses_by_post = defaultdict(list)
    for (pre_m, post_m), group in synapse_groups.items():
        try:
            _, _, f = pruning_df.loc[pre_m, post_m]
        except KeyError:
            f = 1.0  # keep all if (pre_m, post_m) not in pruning_df

        num_keep = int(round(f * len(group)))
        if num_keep < len(group):
            group_sample = random.sample(group, num_keep)
        else:
            group_sample = group

        for post_id_str, syn in group_sample:
            pruned_synapses_by_post[post_id_str].append(syn)

    # Sort by post_id for consistent output
    pruned_synapses_by_post = dict(sorted(pruned_synapses_by_post.items(), key=lambda x: int(x[0])))

    # Save result
    with open(output_json, "w") as f:
        json.dump(pruned_synapses_by_post, f, indent=2)

    if verbose:
        total_syn = sum(len(v) for v in pruned_synapses_by_post.values())
        print(f"\t General pruning complete. Kept {total_syn} synapses.")

    return pruned_synapses_by_post

def setup_02_pruning(model_type_sm_csv:str, goal_type_sm_json:str,
                     save_dir:str = None) -> pd.DataFrame:
    """
    Calculates mu (sigmoid cutoff threshold) for each m_type:m_type based on biological data.
    Saves it into 2d Dataframe and *.csv file.

    I alternated it a bit from (Reimann) as it is really badly explained.
    mu = sm_emp - lambda*sm_emp_std

    Args:
        model_type_sm_csv (str): Pruning metrics coming from the model so far.
        goal_type_sm_json (str): Goal pruning metrics derived from BBP.
        save_dir (str): Directory in which params for 02.Multisynaptic pruning should be saved.
    Returns:
        pd.DataFrame: Parameters used for 02.Multisynaptic pruning
    """

    model_sm = pd.read_csv(model_type_sm_csv, index_col=0)
    with open(goal_type_sm_json, 'r') as f:
        goal_sm = json.load(f)

    mtypes = model_sm.index
    mu_df = pd.DataFrame(index=mtypes, columns=mtypes, dtype=float)

    for pre_m in mtypes:
        for post_m in mtypes:
            pre_m_short = "_".join(pre_m.split('_')[:2])
            post_m_short = "_".join(post_m.split('_')[:2])
            key = f"{pre_m_short}:{post_m_short}"
            bbp_stats = goal_sm.get(key)

            if bbp_stats:
                mean_emp = bbp_stats.get("mean_number_of_synapse_per_connection", 0)
                std_emp = bbp_stats.get("number_of_synapse_per_connection_std", 0)
                mu = mean_emp - conf.LAMBDA_02 * std_emp
                mu_df.at[pre_m, post_m] = round(max(mu, 0.0), 4)
            else:
                mu_df.at[pre_m, post_m] = None  # Or np.nan if preferred

    if save_dir is not None:
        mu_df.to_csv(os.path.join(save_dir, "prune_02_params.csv"))

    return mu_df

def apply_multisyn_pruning_02(
    prune_02_params_csv:str,
    all_cells_csv:str,
    synapses_json:str,
    save_path:str = None) -> defaultdict:
    """
    Applies multi-synaptic pruning using sigmoid cutoff for each individual connection.

    Pruning probability:
        p_prune = 1 / (1 + exp(alpha * (S_conn - mu)))

    Where:
        - S_conn: synapse count of a connection (pre_id, post_id)
        - mu: threshold for (pre_m, post_m) from mu_csv

    Args:
        prune_02_params_csv (str): Path to μ-thresholds (pre_m × post_m).
        all_cells_csv (str): Path to cell_id → m_type mapping CSV.
        synapses_json (str): Path to original synapse JSON (pre-pruning).
        seed (int or None): Random seed.
        save_path (str or None): If given, pruned synapse JSON is saved here.

    Returns:
        defaultdict: Pruned synapse data (same structure as input JSON).
    """

    # Load cell → m_type mapping
    cell_df = pd.read_csv(all_cells_csv)
    cell_id_to_mtype = dict(zip(cell_df['cell_id'], cell_df['desc']))

    # Load μ threshold matrix
    mu_df = pd.read_csv(prune_02_params_csv, index_col=0)

    # Load original synapse JSON
    with open(synapses_json, 'r') as f:
        syn_data = json.load(f)

    # Step 1: Count synapses per (pre_id, post_id)
    conn_to_syns = defaultdict(list)
    for post_id_str, syn_list in syn_data.items():
        post_id = int(post_id_str)
        for syn in syn_list:
            pre_id = syn['pre_id']
            conn_to_syns[(pre_id, post_id)].append(syn)

    # Step 2: For each connection, compute p_prune and decide to prune or not
    pruned_syn_data = defaultdict(list)
    syn_cnt = 0
    for (pre_id, post_id), syn_list in conn_to_syns.items():
        pre_m = cell_id_to_mtype.get(pre_id)
        post_m = cell_id_to_mtype.get(post_id)
        if pre_m is None or post_m is None:
            continue

        mu = mu_df.at[pre_m, post_m]
        if pd.isna(mu):
            mu = 0.0  # fallback if missing data

        S_conn = len(syn_list)
        p_prune = 1.0 / (1 + np.exp(conf.ALPHA_02 * (S_conn - mu)))
        if random.random() > p_prune:
            # Keep the connection (copy all synapses)
            for syn in syn_list:
                syn_cnt += 1
                pruned_syn_data[str(post_id)].append(syn)
        # else: connection is pruned entirely (i.e., no synapses retained)

    # Save pruned JSON if requested
    if save_path is not None:
        with open(save_path, 'w') as f:
            json.dump(pruned_syn_data, f, indent=2)
    print(f"\t Multisynaptic pruning complete. Kept {syn_cnt} synapses.")
    return pruned_syn_data

def setup_03_pruning(examined_pop_csv:str, save_dir:str) -> dict:
    """
    Again, (Riemann) describes this step sooo vaguely, that I am using my own interpretation.
    Calculates a = B_d_bio / B_d_model, and uses it as a uniform sampling prob. to prune WHOLE connection.
    I still wonder whether B_d_model should be CELL-SPECIFIC or M_TYPE-SPECIFIC.
    For now, I will implement it as M_TYPE-SPECIFIC.

    Calculates B_d_model for each m_type based on examined_pop.csv, where b_int (bouton interval) is averaged and
    inversed to get B_d_model.
    B_d_bio comes from BD_EMPIRICAL_DICT.

    Then, saves the results as json at specified path. Key is m_type, within is dict where keys are: bd_bio, bd_mod, retain_ratio

    Args:
        examined_pop_csv (str): Cortical population *csv with calculated b_int.
        save_path (str): Resulting pruned synapses *json.
    """
    df = pd.read_csv(examined_pop_csv)

    if "lay_m_type" not in df.columns or "b_int" not in df.columns:
        raise ValueError("CSV must contain 'lay_m_type' and 'b_int' columns.")

    result = {}

    for lay_m_type, group in df.groupby("lay_m_type"):

        if lay_m_type not in conf.BD_EMPIRICAL_DICT:
            print(f"\t Warning: No empirical B_d for m-type '{lay_m_type}', assigning B_d={conf.REP_BD_EMP}.")
            continue

        b_int_values = group["b_int"].replace([np.inf, -np.inf], np.nan).dropna()
        if len(b_int_values) == 0:
            print(f"\t Warning: No bouton interval values for '{lay_m_type}', skipping.")
            continue

        b_int_mean = b_int_values.mean()
        bd_model = 1.0 / b_int_mean if b_int_mean != 0 else 0.0
        bd_bio = conf.BD_EMPIRICAL_DICT[lay_m_type]
        retain_ratio = bd_bio / bd_model if bd_model > 0 else 0.0

        result[lay_m_type] = {
            "bd_bio": round(bd_bio, 5),
            "bd_model": round(bd_model, 5),
            "retain_ratio": round(min(retain_ratio, 1.0), 5)  # cap at 1.0
        }

    # Save to JSON
    save_path = os.path.join(save_dir, "prune_03_params.json")
    with open(save_path, 'w') as f:
        json.dump(result, f, indent=4)

    return result

def apply_plastres_pruning_03(
    post02_synapse_json:str,
    prune_03_params_json:str,
    save_path:str) -> defaultdict:
    """
    Applies 03.plasticity-reserve pruning. Saves resulting synapses in json.

    Args:
        post02_synapse_json (str): Synapses in json after second pruning.
        prune_03_params_json (str): Calculated parameters in json for third pruning.
        save_path (str): Save path.
    """
    # Load pruning parameters
    with open(prune_03_params_json, "r") as f:
        retain_info = json.load(f)

    # Map: pre_m_type → retain_ratio
    retain_ratios = {
        mtype: info.get("retain_ratio", 1.0)
        for mtype, info in retain_info.items()
    }

    # Load synapses
    with open(post02_synapse_json, "r") as f:
        synapse_dict = json.load(f)

    pruned_synapses_by_post = defaultdict(list)
    syn_cnt = 0

    for post_id, syn_list in synapse_dict.items():
        for syn in syn_list:
            syn_cnt += 1
            pre_type = syn.get("pre_me_type")
            retain_prob = retain_ratios.get(pre_type, 1.0)  # default: retain
            if random.random() < retain_prob:
                pruned_synapses_by_post[post_id].append(syn)

    # Save result
    with open(save_path, "w") as f:
        json.dump(pruned_synapses_by_post, f, indent=2)

    print(f"\t Plasticity-reserve pruning complete. Kept {syn_cnt} synapses.")
    return pruned_synapses_by_post