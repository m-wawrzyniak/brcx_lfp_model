import csv, json

from source.src0_core.cr1_simulation_run.s02_cells_init.CxCell import CxCell

from source.src2_utils.ut0_random_manager import np

def init_cx_cells(cx_cells_csv:str, temp_dict:dict, me_comp_summary_path) -> dict[str, CxCell]:
    """
    Initializes all cx_cells specified in cx_cells_csv. For this task, it utilizes loaded cell templates in temp_dict.

    Args:
        cx_cells_csv (str): Path to CX cells *csv.
        temp_dict (dict): Dictionary with all necessary CX cell templates. Result of b01.load_all_templates()

    Returns:
        dict: <cell_id: CxCell()>

    """
    print("[initc01] Initializing cx cells.")

    cells = {}
    with open(me_comp_summary_path, 'r') as f:
        me_summary = json.load(f)
        population_data = me_summary["population"]
        n_cells = population_data["tot_n_cells"]

    with open(cx_cells_csv, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        cnt = 0
        for row in reader:
            if (cnt + 1) % 10 == 0 or (cnt + 1) == n_cells:
                print(f'\t\t CX cells instantiated: {(cnt + 1) / n_cells:.1%}')
            cnt += 1
            cell_id = row['cell_id']
            desc = row['desc']
            x = float(row['x'])
            y = float(row['y'])
            z = float(row['z'])
            rot_ang = float(row['rot_ang'])

            abb_name = "_".join(desc.split('_')[:2])
            new_cell = CxCell(cell_id=cell_id, cell_name=desc, cell_temp=temp_dict[abb_name], offset=(x, y, z), rotation=(np.pi/2, 0, rot_ang))
            cells[cell_id] = new_cell

    print(f"[initc01] SUCCESS: Initialized {cnt} cx cells.")

    return cells