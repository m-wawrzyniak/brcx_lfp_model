import csv, json, os
from neuron import h

from source.src0_core.cr1_simulation_run.s02_cells_init.CxCell import CxCell
from source.src0_core.cr1_simulation_run.s03_conn_init.CxSynapse import CxSynapse

from source.src2_utils.ut0_random_manager import np

import config_templates.conf01_simulation_parameters as conf1

def get_synapse_parameters(cx_syns_summary_csv:str ,
                           bbp_physio_file:str) -> dict[str, dict[str, float]]:
    """
    Based on 'synapse_summary.json', fetches all the parameters included in 'pathways_physiology_factsheets_simplified.json'
    and returns it as dict.

    Args:
        cx_syns_summary_csv (str): Path to *.csv with summary of existing intracortical synapses.
        bbp_physio_file (str): Path to BBP file with physiology of the connections.

    Returns:
        dict : <conn_type: dict_with_params>, where conn_type e.g. L4_SS_cAD:L23_LBC_cNAC.

    """
    # Load files
    with open(cx_syns_summary_csv, 'r') as f:
        summary_data = json.load(f)

    with open(bbp_physio_file, 'r') as f:
        physio_data = json.load(f)

    # Relevant parameter keys (mean and std)
    relevant_keys = [
        'u_mean', 'u_std',
        'd_mean', 'd_std',
        'f_mean', 'f_std',
        'risetime_mean', 'risetime_std',
        'decay_mean', 'decay_std',
        'gsyn_mean', 'gsyn_std',
        'latency_mean', 'latency_std'
    ]

    # Final output dictionary
    res_dict = {}

    for conn_type in summary_data:
        if conn_type in physio_data:
            physio = physio_data[conn_type]
            filtered_params = {k: physio[k] for k in relevant_keys if k in physio}
            res_dict[conn_type] = filtered_params

            if physio_data[conn_type]['synapse_type'].split(',')[0] == 'Excitatory':
                syn_type = 'e'
            else:
                syn_type = 'i'
            res_dict[conn_type]['type'] = syn_type
        else:
            print(f"Warning: No physiology data for connection {conn_type}")

    return res_dict

def _create_cxcx_synapse(pre_cell:CxCell, post_cell:CxCell,
                         post_loc_rel:tuple[float,float,float],
                         syn_type:str, syn_params:dict, syn_id:str=None,
                         deterministic:bool=False) -> CxSynapse:
    """
    Uses post_cell.segs_meta for O(1) synapse placement.
    """
    pre_hcell = pre_cell.h_cell
    post_abs = np.array(post_cell.center) + np.array(post_loc_rel)

    # Find closest segment
    min_dist = float("inf")
    target_sec = None
    target_loc = 0.5
    target_xyz = None

    for seg_meta in post_cell.segs_meta:
        x, y, z = seg_meta["xyz_um"]
        sec = seg_meta["sec"]
        dist = np.linalg.norm(post_abs - np.array([x, y, z]))
        if dist < min_dist:
            min_dist = dist
            target_sec = sec
            target_loc = seg_meta["x"]
            target_xyz = seg_meta["xyz_um"]

    if target_sec is None:
        raise RuntimeError("Could not find a section for synapse placement.")

    syn_template = 'ProbAMPANMDA_EMS' if syn_type == 'e' else 'ProbGABAAB_EMS'
    syn = getattr(h, syn_template)(target_loc, sec=target_sec)

    def pick_and_clip(val_key, min_val=1e-3):
        if deterministic:
            return max(syn_params[val_key], min_val)

        mean_val = syn_params[val_key]
        std_val = syn_params.get(val_key.replace("mean", "std"), 0)

        if mean_val <= 0 or std_val <= 0:
            # Degenerate or invalid values — fallback to min_val
            return min_val

        # Convert to log-normal parameters
        variance = std_val ** 2
        mu = np.log((mean_val ** 2) / np.sqrt(variance + mean_val ** 2))
        sigma = np.sqrt(np.log(1 + (variance / (mean_val ** 2))))

        val = np.random.lognormal(mean=mu, sigma=sigma)
        return max(val, min_val)

    # For u, d, f (Use, Dep, Fac) values, clamp minimum to 0 (can be zero or small positive)
    u_val = pick_and_clip('u_mean', min_val=0)
    d_val = pick_and_clip('d_mean', min_val=0)
    f_val = pick_and_clip('f_mean', min_val=0)

    # For time constants and conductances, clamp minimum to 1 (or small positive)
    risetime_val = pick_and_clip('risetime_mean', min_val=1)
    decay_val = pick_and_clip('decay_mean', min_val=1)
    gsyn_val = pick_and_clip('gsyn_mean', min_val=1e-3)

    late_comp_ratio = 0.0 if syn_type == 'e' else 0.5 #syn_params.get('late_comp_ratio', 0.0 if syn_type == 'e' else 0.0)

    if syn_type == 'e':
        syn.Use, syn.Dep, syn.Fac = u_val, d_val, f_val
        syn.tau_r_AMPA, syn.tau_d_AMPA = risetime_val, decay_val
        syn.gmax = gsyn_val/6  # TODO what if gsyn is cumulative for a WHOLE CONNECTION, not singular synapse
        syn.NMDA_ratio = late_comp_ratio
    else:
        syn.Use, syn.Dep, syn.Fac = u_val, d_val, f_val
        syn.tau_r_GABAA, syn.tau_d_GABAA = risetime_val, decay_val
        syn.gmax = gsyn_val/6
        syn.GABAB_ratio = late_comp_ratio

    nc = h.NetCon(pre_hcell.soma[0](0.5)._ref_v, syn, sec=pre_hcell.soma[0])
    nc.threshold = -20
    nc.delay = syn_params.get('latency_mean', 1.0)

    # Choose weight
    # TODO: NEW WEIGHTS SYSTEM
    pre_index = conf1.CELL_ORDER_MAP[pre_cell.cell_name]
    post_index = conf1.CELL_ORDER_MAP[post_cell.cell_name]
    nc.weight[0] = conf1.WEIGHT_MATRIX[pre_index][post_index]


    return CxSynapse(
        syn_obj=syn,
        netcon=nc,
        pre_id=pre_cell.cell_id,
        post_id=post_cell.cell_id,
        pre_me_type=pre_cell.cell_name,
        post_me_type=post_cell.cell_name,
        post_loc_rel=post_loc_rel,
        syn_type=syn_type,
        u_val=u_val,
        d_val=d_val,
        f_val=f_val,
        risetime_val=risetime_val,
        decay_val=decay_val,
        gsyn_val=gsyn_val,
        late_comp_ratio=late_comp_ratio,
        syn_id = syn_id,
        xyz=target_xyz
    )

def create_cxcx_synapses(cx_cells:dict[str, CxCell], syn_params:dict,
                         cxcx_synapses_data:str) -> dict[str, CxSynapse]:
    """
    Based on specified intracortical synapses (cxcx_synapses_data), creates them using _create_cxcx_synapse().
    Returns the synapses as a dictionary.

    Args:
        cx_cells (dict): <cell_id: CxCell()>
        syn_params (dict): <param: value>
        cxcx_synapses_data (str): Path to *json with all specified intracortical synapses to be made.

    Returns:
        dict : Synapse dictionary <syn_id: CxSynapse()>

    """
    with open(cxcx_synapses_data, 'r') as f:
        syn_defs = json.load(f)

    synapse_map = {}
    n_post = len(syn_defs)
    post_int = 0

    for post_id, pre_list in syn_defs.items():
        if (post_int+1)%3 == 0 or (post_int+1)==n_post:
            print(f'\t\t IntraCX post-cell connected: {post_int / n_post:.1%}')
        post_int += 1

        for pre_dict in pre_list:
            pre_id = pre_dict['pre_id']
            pre_type = pre_dict['pre_me_type']
            post_type = pre_dict['post_me_type']
            post_loc_rel = pre_dict['post_loc']

            pre_short = '_'.join(pre_type.split('_')[:-1])
            post_short = '_'.join(post_type.split('_')[:-1])

            conn_key = f"{pre_short}:{post_short}"
            if conn_key not in syn_params:
                print(f"\t\t [Warning] No synaptic parameters found for {conn_key}. Skipping synapse {pre_id}:{post_id}.")
                continue

            params = syn_params[conn_key]

            try:
                pre_cell = cx_cells[str(pre_id)]
                post_cell = cx_cells[str(post_id)]
                syn_type = params['type']  # 'e' or 'i'

                syn_obj = _create_cxcx_synapse(pre_cell, post_cell, post_loc_rel, syn_type, params)
                synapse_map[syn_obj.syn_id] = syn_obj

            except Exception as e:
                print(f"[Error] Failed to create synapse {pre_id}:{post_id}: {e}")

    return synapse_map

def save_synapses_to_csv(synapse_map:dict[str, CxSynapse], save_path:str):
    """
    Saves parameters of all intracortical synapses to CSV.

    Args:
        synapse_map (dict): <syn_id: CxSynapse()>
        save_path (str): Where the resulting *csv should be saved.

    """
    fieldnames = [
        "syn_id",
        "pre_id",
        "post_id",
        "pre_me_type",
        "post_me_type",
        "post_loc_x",
        "post_loc_y",
        "post_loc_z",
        "abs_x",
        "abs_y",
        "abs_z",
        "syn_type",
        "u_val",
        "d_val",
        "f_val",
        "risetime_val",
        "decay_val",
        "gsyn_val",
        "late_comp_ratio",
        "weight"
    ]

    with open(save_path, mode='w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for syn_id, syn_obj in synapse_map.items():
            # Extract parameters
            params = syn_obj.get_params()

            # Recover IDs and types (you might want to store these in syn_obj on creation)
            pre_id = None
            post_id = None
            pre_me_type = None
            post_me_type = None
            post_loc = None
            syn_type = None

            # Example: assuming you stored these in syn_obj or can infer from your data structure
            # You can extend CxSynapse to store these at creation, or store in syn_obj.extra dict

            if hasattr(syn_obj, 'pre_id'):
                pre_id = syn_obj.pre_id
            if hasattr(syn_obj, 'post_id'):
                post_id = syn_obj.post_id
            if hasattr(syn_obj, 'pre_me_type'):
                pre_me_type = syn_obj.pre_me_type
            if hasattr(syn_obj, 'post_me_type'):
                post_me_type = syn_obj.post_me_type
            if hasattr(syn_obj, 'post_loc_rel'):
                post_loc = syn_obj.post_loc_rel
            if hasattr(syn_obj, 'xyz'):
                xyz = syn_obj.xyz
            if hasattr(syn_obj, 'syn_type'):
                syn_type = syn_obj.syn_type
            if hasattr(syn_obj, 'nc'):
                syn_weight = syn_obj.nc.weight[0]

            row = {
                "syn_id": syn_id,
                "pre_id": pre_id,
                "post_id": post_id,
                "pre_me_type": pre_me_type,
                "post_me_type": post_me_type,
                "post_loc_x": post_loc[0] if post_loc else None,
                "post_loc_y": post_loc[1] if post_loc else None,
                "post_loc_z": post_loc[2] if post_loc else None,
                "abs_x": xyz[0],
                "abs_y": xyz[1],
                "abs_z": xyz[2],
                "syn_type": syn_type,
                "u_val": params["u_val"],
                "d_val": params["d_val"],
                "f_val": params["f_val"],
                "risetime_val": params["risetime_val"],
                "decay_val": params["decay_val"],
                "gsyn_val": params["gsyn_val"],
                "late_comp_ratio": params["late_comp_ratio"],
                "weight": syn_weight
            }
            writer.writerow(row)


