import json
from neuron import h

import config_templates.conf01_simulation_parameters as conf01

from source.src0_core.cr1_simulation_run.s03_conn_init.PrvTcSynapse import PrvTcSynapse

def create_prvtc_synapses(tc_cells:dict, prv_json_path:str):
    """
    Create PrV->VPM synapses from JSON specification.

    Parameters
    ----------
    tc_cells : dict
        Mapping {cell_id: VPMCell_object} with `soma` section.
    prv_json_path : str
        Path to JSON file containing PrV cell spike_times and target_vpm.
    """
    print("[inits03] Creating prvtc synapses.")

    # Load PrV JSON
    with open(prv_json_path, 'r') as f:
        prv_cells_data = json.load(f)

    prv_vpm_synapses = {}
    syn_cnt=0
    for prv_id, prv_info in prv_cells_data.items():
        syn_cnt += 1
        target_vpm_id = prv_info['target_vpm']
        spike_times = prv_info['spike_times']

        # Get the target VPM cell
        if target_vpm_id not in tc_cells:
            raise KeyError(f"Warning: Target VPM cell {target_vpm_id} not found in vpm_cells")

        vpm_cell = tc_cells[target_vpm_id]

        # Create synapse at soma(0.5)
        syn = h.Exp2Syn(vpm_cell.soma(0.5))

        # Create VecStim for PrV spikes
        vec = h.Vector(spike_times)
        stim = h.VecStim()
        stim.play(vec)

        # NetCon from VecStim to Exp2Syn
        nc = h.NetCon(stim, syn)
        nc.delay = 0
        pre_index = conf01.CELL_ORDER_MAP['PRV']
        post_index = conf01.CELL_ORDER_MAP['VPM']
        nc.weight[0] = conf01.WEIGHT_MATRIX[pre_index][post_index]


        # Wrap into PrVTcSynapse
        prvtc_syn = PrvTcSynapse(syn_obj=syn, nc=nc, vestim=stim, spike_vec=spike_times)

        prv_vpm_synapses[prv_id] = prvtc_syn

    print(f"[inits03] SUCCESS: Created prvtc synapses. Count = {syn_cnt}")

    return prv_vpm_synapses
