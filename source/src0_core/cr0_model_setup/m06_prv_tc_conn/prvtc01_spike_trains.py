import pandas as pd
import os, json

from source.src2_utils.ut0_random_manager import rng, np
import config_templates.conf0_model_parameters as conf0

from source.src0_core.cr0_model_setup.m05_prv_cells.PrvCell import PrvCell

def generate_prv_spikes(stim_paradigm_type, stim_paradigm_subtype, tc_cells_path, prvtc_save_path):
    print('p01: Setup of PrV cells.')
    ### SETUP
    # Instantiating PrV cells
    print('\t Instantiating PrV cells and creating deflection-response spike trains')
    vpm_cells = pd.read_csv(tc_cells_path)
    vpm_ids = vpm_cells['cell_id'].tolist()

    n_prv = conf0.PRV_PER_VPM_CELL * len(vpm_ids)
    prv_cells = {}
    cell_gen = rng
    for i in range(n_prv):
        prv_cell = PrvCell(cell_gen)
        prv_cells[prv_cell.cell_id] = {'cell_obj':prv_cell}

    # Creating deflection-response spike trains
    for prv_id, prv_dict in prv_cells.items():
        c = prv_dict['cell_obj']
        spike_times = c.set_spike_sequence(
            paradigm_type=stim_paradigm_type,
            paradigm_subtype=stim_paradigm_subtype)
        prv_dict['spike_times'] = list(spike_times)

    # Assigning VPM cells to PrV cell
    print('\t Assigning VPM cells to each PrV cell')
    assignments = np.repeat(vpm_ids, conf0.PRV_PER_VPM_CELL)
    cell_gen.shuffle(assignments)

    for prv_id, target_vpm in zip(prv_cells.keys(), assignments):
        prv_cells[prv_id]['target_vpm'] = target_vpm

    # Saving synapses parametrization
    print('\t Saving setup results.')
    prv_cells_json = {
        prv_id: {k: v for k, v in prv_dict.items() if k != 'cell_obj'}
        for prv_id, prv_dict in prv_cells.items()
    }

    with open(prvtc_save_path, 'w') as f:
        json.dump(prv_cells_json, f, indent=2)
