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
    vpm_cells = {}

    with open(vpm_cells_path, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            cell_id = row["cell_id"]
            vpm_cells[cell_id] = VPMCell(cell_id)  # Pass additional args if needed

    print(f"\t v03: Instantiated{len(vpm_cells)} VPMCell objects.")
    return vpm_cells