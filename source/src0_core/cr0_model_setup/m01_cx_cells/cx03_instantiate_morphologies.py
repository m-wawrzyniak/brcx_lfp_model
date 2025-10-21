import os
import glob
import pandas as pd

from source.src2_utils.ut0_random_manager import np
from source.src0_core.cr0_model_setup.m01_cx_cells.CellTopology import CellTopology


def load_population(csv_path:str) -> dict:
    """
    Loads *_me_pop.csv and extracts all the important information for defining the intracortical connectivity.
    Args:
        csv_path (str): Path to cortical cells *.csv file

    Returns:
        dict: Keyed with 'cell_id' , valued with (dict) containing pos, layer, e_type, m_type and rot_ang.
    """
    df = pd.read_csv(csv_path)
    all_cells = {}
    for _, row in df.iterrows():
        all_cells[int(row['cell_id'])] = {
            'pos': (row['x'], row['y'], row['z']),
            'layer': row['layer'],
            'e_type': row['e_type'],
            'm_type': row['m_type'],
            'rot_ang': (np.pi/2, 0, row['rot_ang'])
        }
    return all_cells

def instantiate_morphologies(all_cells:dict , morph_dir:str) -> dict:
    """
    Instantiates morphology of each cortical cell and extends the all_cells (dict) with axon and dendrite points later
    used for finding appositions.
    Extended dictionary is returned.
    """
    n_cells = len(all_cells)

    for c_id, c_info in all_cells.items():
        if (c_id+1)%10 == 0 or (c_id+1)==n_cells:
            print(f'\t Morphologies instantiated: {(c_id + 1) / n_cells:.1%}')

        cell_morph_path = os.path.join(
            morph_dir, f"{c_info['layer']}_{c_info['m_type']}_{c_info['e_type']}_minimal/morphology/"
        )
        asc_file = glob.glob(os.path.join(cell_morph_path, "*.asc"))
        asc_file = asc_file[0]

        topo = CellTopology(asc_file, offset=c_info['pos'], rotation=c_info['rot_ang'])


        # save both the flat points and the CellTopology object itself
        all_cells[c_id]['axon_pts'] = topo.get_axon_points()
        all_cells[c_id]['dendrite_pts'] = topo.get_dendrite_points()
        all_cells[c_id]['axon_len'] = topo.axon_len
        all_cells[c_id]['topo'] = topo   # <-- NEW LINE

    return all_cells

def save_all_cells(cell_dict:dict, output_file:str):
    """
    After apositions have been found, saves the all cells in cell_dict to new *csv file with all data that need to be
    propagated.
    Args:
        cell_dict (dict):  Dictionary of cortical cells, keyed with cell_id (str).
        output_file (str): Path where *.csv file should be saved.
    """
    rows = []
    for cell_id, cell_data in cell_dict.items():
        row = {'cell_id': cell_id}
        for key, value in cell_data.items():
            if key in ('axon_pts', 'dendrite_pts'):
                continue
            elif key == 'pos':
                row['x'], row['y'], row['z'] = value
            elif key == 'rot_ang':
                row['rot_ang'] = value[2]  # only z-axis rotation
            else:
                row[key] = value

        # Create 'desc' field
        row['lay_m_type'] = f"{row['layer']}_{row['m_type']}"
        row['desc'] = f"{row['layer']}_{row['m_type']}_{row['e_type']}"
        rows.append(row)

    df = pd.DataFrame(rows)
    df.to_csv(output_file, index=False)