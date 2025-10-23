import json, csv, os
from neuron import h

from source.src2_utils.ut0_random_manager import np

from source.src0_core.cr1_simulation_run.s02_cells_init.CxCell import CxCell
from source.src0_core.cr1_simulation_run.s02_cells_init.VPMCell import VPMCell

from source.src0_core.cr1_simulation_run.s03_conn_init.TcCxSynapse import TcCxSynapse

import config_templates.conf01_simulation_parameters as conf01

def _create_tccx_synapse(pre_vpm_cell: VPMCell, post_cell: CxCell,
                         post_loc_rel: tuple[float,float,float] = (0,0,0),
                         syn_type: str = 'e', syn_params: dict = conf01.TCCX_SYNAPSE_PARAMS,
                         deterministic: bool = False) -> TcCxSynapse:
    """
    Place a TC->CX synapse using post_cell.segs_meta (O(1) lookup) without relying on section names.
    """

    post_abs = np.array(post_cell.center) + np.array(post_loc_rel)

    # Find closest segment
    min_dist = float("inf")
    target_sec = None
    target_loc = 0.5

    for seg_meta in post_cell.segs_meta:
        x, y, z = seg_meta["xyz_um"]
        sec = seg_meta["sec"]
        dist = np.linalg.norm(post_abs - np.array([x, y, z]))
        if dist < min_dist:
            min_dist = dist
            target_sec = sec
            target_loc = seg_meta.get("x", 0.5)

    if target_sec is None:
        raise RuntimeError(f"Could not find a segment for TC->CX synapse on post cell {post_cell.cell_id}")

    # Create synapse mechanism
    syn_template = 'ProbAMPANMDA_EMS'
    syn = getattr(h, syn_template)(target_loc, sec=target_sec)

    # Assign parameters
    u_val = np.random.normal(syn_params['u_mean'], syn_params['u_std'])
    d_val = np.random.normal(syn_params['d_mean'], syn_params['d_std'])
    f_val = np.random.normal(syn_params['f_mean'], syn_params['f_std'])
    risetime_val = np.random.normal(syn_params['risetime_mean'], syn_params['risetime_std'])
    decay_val = np.random.normal(syn_params['decay_mean'], syn_params['decay_std'])
    gsyn_val = np.random.normal(syn_params['gsyn_mean'], syn_params['gsyn_std'])
    late_comp_ratio = conf01.TCCX_NMDA_COMPONENT

    syn.Use = u_val
    syn.Dep = d_val
    syn.Fac = f_val
    syn.tau_r_AMPA = risetime_val
    syn.tau_d_AMPA = decay_val
    syn.gmax = gsyn_val
    syn.NMDA_ratio = late_comp_ratio

    # Create NetCon
    nc = h.NetCon(pre_vpm_cell.soma(0)._ref_v, syn, sec=pre_vpm_cell.soma)
    nc.threshold = conf01.TCCX_NC_THRESH
    nc.delay = np.random.normal(syn_params['latency_mean'], syn_params['latency_std'])

    pre_index = conf01.CELL_ORDER_MAP['VPM']
    post_index = conf01.CELL_ORDER_MAP[post_cell.cell_name]
    nc.weight[0] = conf01.WEIGHT_MATRIX[pre_index][post_index]


    return TcCxSynapse(
        syn_obj=syn,
        netcon=nc,
        pre_id=pre_vpm_cell.cell_id,
        post_id=post_cell.cell_id,
        pre_me_type=pre_vpm_cell.cell_name,
        post_me_type=post_cell.cell_name,
        syn_type=syn_type,
        u_val=u_val,
        d_val=d_val,
        f_val=f_val,
        risetime_val=risetime_val,
        decay_val=decay_val,
        gsyn_val=gsyn_val,
        late_comp_ratio=late_comp_ratio,
        post_loc_rel=post_loc_rel
    )

def create_tccx_synapses(tc_cells: dict[str, VPMCell], cx_cells: dict[str, CxCell],
                         tccx_synapses_data: str ) -> dict[str, TcCxSynapse]:
    """
    Optimized creation of all TC-CX synapses using segs_meta lookup.

    Args:
        tc_cells: <cell_id: VPMCell>
        cx_cells: <cell_id: CxCell>
        tccx_synapses_data: path to JSON with synapse definitions
    Returns:
        dict: <syn_id: TcCxSynapse>
    """
    print("[inits02] Creating tccx synapses.")

    with open(tccx_synapses_data, 'r') as f:
        syn_defs = json.load(f)

    synapse_map = {}
    n_post = len(syn_defs)
    post_int = 0

    for post_id, pre_list in syn_defs.items():
        if post_int % 3 == 0 or post_int == n_post:
            print(f'\t tccx cells connected: {post_int / n_post:.1%}')
        post_int += 1

        post_cell = cx_cells[str(post_id)]

        for pre_dict in pre_list:
            pre_id = pre_dict['pre_id']
            post_loc_raw = pre_dict.get('post_loc', (0.0, 0.0, 0.0))

            # Ensure post_loc is a numeric tuple (x, y, z)
            if isinstance(post_loc_raw, str):
                # Parse string like "[1.0, 2.0, 3.0]"
                post_loc = tuple(map(float, post_loc_raw.strip("[]").split(",")))
            elif isinstance(post_loc_raw, list) or isinstance(post_loc_raw, tuple):
                post_loc = tuple(float(x) for x in post_loc_raw)
            else:
                # Single numeric value fallback
                post_loc = (float(post_loc_raw), 0.0, 0.0)

            pre_cell = tc_cells[str(pre_id)]

            try:
                syn_obj = _create_tccx_synapse(
                    pre_vpm_cell=pre_cell,
                    post_cell=post_cell,
                    post_loc_rel=post_loc
                )
                synapse_map[syn_obj.syn_id] = syn_obj
            except Exception as e:
                print(f"\t ERR: Failed to create synapse {pre_id}:{post_id}: {e}")

    print("[inits02] SUCCESS: Created all tccx synapses.")
    return synapse_map

def save_synapses_to_csv(synapse_map:dict[str, TcCxSynapse], save_path:str):
    """
    Save synapse parameters to CSV.

    Args:
        synapse_map (dict): <syn_id, TcCxSynapse>. Result from create_tccx_synapses().
        save_path (str): Path where resulting *csv should be saved.
    """
    print("[inits02] Saving tccx synapses with electrophysiological parameters.")

    fieldnames = [
        "syn_id",
        "pre_id",
        "post_id",
        "pre_me_type",
        "post_me_type",
        "post_loc_x",
        "post_loc_y",
        "post_loc_z",
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
        cnt = 0
        for syn_id, syn_obj in synapse_map.items():
            # Extract parameters
            cnt+=1
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

    print(f"[inits02] SUCCESS: Saved tccx synapses with electrophysiological parameters. Count = {cnt}")