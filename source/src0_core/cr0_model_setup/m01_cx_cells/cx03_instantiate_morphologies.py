import os
import glob
import pandas as pd

from source.src2_utils.ut0_random_manager import np
from source.src0_core.cr0_model_setup.m01_cx_cells.CellTopology import CellTopology


def load_population(csv_path:str) -> dict:
    print(f"[cx03] Loading layer population from {csv_path}")
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
    print(f"[cx03] SUCCESS: Layer population loaded.")
    return all_cells

def instantiate_morphologies(all_cells:dict , morph_dir:str) -> dict:
    print("[cx03] Instantiating cx morphologies")

    n_cells = len(all_cells)

    for c_id, c_info in all_cells.items():
        if (c_id+1)%10 == 0 or (c_id+1)==n_cells:
            print(f'\t Cells instantiated: {(c_id + 1) / n_cells:.1%}')

        cell_morph_path = os.path.join(
            morph_dir, f"{c_info['layer']}_{c_info['m_type']}_{c_info['e_type']}_minimal/morphology/"
        )
        asc_file = glob.glob(os.path.join(cell_morph_path, "*.asc"))
        asc_file = asc_file[0]

        topo = CellTopology(asc_file, offset=c_info['pos'], rotation=c_info['rot_ang'])

        all_cells[c_id]['axon_pts'] = topo.get_axon_points()
        all_cells[c_id]['dendrite_pts'] = topo.get_dendrite_points()
        all_cells[c_id]['axon_len'] = topo.axon_len
        all_cells[c_id]['topo'] = topo
    print("[cx03] SUCCESS: cx morphologies instantiated")

    return all_cells

def save_all_cells(cell_dict:dict, output_file:str):
    print("[cx03] Saving cx cells after morphology initialization.")

    rows = []
    for cell_id, cell_data in cell_dict.items():
        row = {'cell_id': cell_id}
        for key, value in cell_data.items():
            if key in ('axon_pts', 'dendrite_pts'):
                continue
            elif key == 'pos':
                row['x'], row['y'], row['z'] = value
            elif key == 'rot_ang':
                row['rot_ang'] = value[2]
            else:
                row[key] = value

        row['lay_m_type'] = f"{row['layer']}_{row['m_type']}"
        row['desc'] = f"{row['layer']}_{row['m_type']}_{row['e_type']}"
        rows.append(row)

    df = pd.DataFrame(rows)
    df.to_csv(output_file, index=False)

    print("[cx03] SUCCESS: cx cells saved.")