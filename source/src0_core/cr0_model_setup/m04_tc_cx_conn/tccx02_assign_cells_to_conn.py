import json, csv

from source.src2_utils.ut0_random_manager import random

def assign_tc_synapses(preassignment_synapses_path:str, assigned_synapses_path:str, tc_pop_path:str):
    """
    When both the number of VPM cells and TCCX synapses are decided, this function assigns VPM cell to each of the TCCX synapse,
    which already has its postcell target.

    Args:
        preassignment_synapses_path (str): Input json file containing synapses with no assigned precell.
        assigned_synapses_path (str): Output json file with assigned precell.

    """
    # Load VPM cells
    with open(tc_pop_path, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        vpm_ids = [row['cell_id'] for row in reader]

    if not vpm_ids:
        raise ValueError("No VPM cells found in CSV.")

    # Load initial synapses
    with open(preassignment_synapses_path, 'r') as f:
        syn_data = json.load(f)

    # Assign synapses
    updated_syn_data = {}

    for post_cell_id, syn_list in syn_data.items():
        updated_syn_list = []
        for syn in syn_list:
            vpm_cell_id = random.choice(vpm_ids)
            syn['pre_id'] = vpm_cell_id
            syn['pre_me_type'] = "VPM_standard_standard"
            updated_syn_list.append(syn)
        updated_syn_data[post_cell_id] = updated_syn_list

    # Save updated synapses
    with open(assigned_synapses_path, 'w') as f:
        json.dump(updated_syn_data, f, indent=2)

    print(f"\t Saved updated VPM-assigned synapses")