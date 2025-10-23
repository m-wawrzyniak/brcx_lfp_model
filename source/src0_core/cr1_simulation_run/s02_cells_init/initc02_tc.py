import csv

from source.src0_core.cr1_simulation_run.s02_cells_init.VPMCell import VPMCell

def init_vpm_cells(vpm_cells_path:str) -> dict[str, VPMCell]:
    """
    Initializes VPMCell() object for each cell noted in vpm_cells_path.
    Args:
        vpm_cells_path (str): Path to csv file containing data about VPM cells.
    Returns:
        dict : In format <vpm_cell_id : VPMCell()>
    """
    print("[initc02] Initializing tc cells.")

    vpm_cells = {}

    with open(vpm_cells_path, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        cnt = 0
        for row in reader:
            cnt+=1
            cell_id = row["cell_id"]
            vpm_cells[cell_id] = VPMCell(cell_id)

    print(f"[initc02] SUCCESS: Initialized {cnt} tc cells.")

    return vpm_cells