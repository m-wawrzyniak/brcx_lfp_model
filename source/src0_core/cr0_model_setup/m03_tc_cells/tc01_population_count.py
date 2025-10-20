import os, json, math, csv

import config_templates.conf0_model_parameters as conf0

def setup_vpm_pop(tccx_synapses_path:str, tc_cells_path:str, syn_per_tc_cell: int = conf0.SYN_PER_TC_CELL):
    """
    Calculates and creates VPM cells csv file, which will later be used to initialize them during the simulation.
    Resulting csv file has only one column: cell_id (str)

    Args:
        tccx_synapses_path (str): Path to all designated TC-CX synapses. Json format.
        vpm_cells_path (str): Path where chosen VPM cells will be saved . Csv format.
        syn_per_tc_cell (int): Number of synapses a single TC cell fiber creates within cortical column.
    """
    if not os.path.exists(tccx_synapses_path):
        raise FileNotFoundError(f"TC-CX synapse file not found at: {tccx_synapses_path}")

    with open(tccx_synapses_path, 'r') as f:
        tc_syn_data = json.load(f)

    total_syn = sum(len(syn_list) for syn_list in tc_syn_data.values())
    n_vpm_cells = math.ceil(total_syn / syn_per_tc_cell)

    with open(tc_cells_path, mode='w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["cell_id"])  # Header
        for i in range(1, n_vpm_cells + 1):
            writer.writerow([f"v_{i}"])

    print(f"\t Based on number of TC-CX synapses,  {n_vpm_cells} VPM cells will be created.")