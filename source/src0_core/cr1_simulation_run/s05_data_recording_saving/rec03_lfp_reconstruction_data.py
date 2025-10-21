import h5py

from source.src3_utils.ut0_random_manager import np


def start_recording_all_imem(cells_dict):
    """
    Iterate over dict of {cell_id: CxCell} and start recording i_mem for each.
    Returns nothing; recordings are stored inside each CxCell.
    """
    for cid, cell in cells_dict.items():
        cell.start_recording_imem()


def save_all_imem_h5(cells_dict, filename):
    """
    Save transmembrane currents from all cells into an HDF5 file.
    """
    with h5py.File(filename, "w") as f:
        # take time vector from first cell
        first_cell = next(iter(cells_dict.values()))
        _, t_ms, _ = first_cell.get_imem_array()
        f.create_dataset("time", data=np.asarray(t_ms, dtype=np.float64))

        cells_grp = f.create_group("cells")

        for cid, cell in cells_dict.items():
            imem_A, t_ms, segs_meta = cell.get_imem_array()

            grp = cells_grp.create_group(str(cid))
            grp.create_dataset("imem", data=np.asarray(imem_A, dtype=np.float64))

            meta_dtype = np.dtype([
                ("x_um", "f8"), ("y_um", "f8"), ("z_um", "f8"),
                ("area_cm2", "f8"), ("length_um", "f8"),
                ("stype", h5py.string_dtype(encoding="utf-8"))
            ])

            meta_arr = np.zeros(len(segs_meta), dtype=meta_dtype)
            for i, m in enumerate(segs_meta):
                x, y, z = m["xyz_um"]
                meta_arr[i] = (x, y, z, m["area_cm2"], m["length_um"], m["stype"])
            grp.create_dataset("meta", data=meta_arr)