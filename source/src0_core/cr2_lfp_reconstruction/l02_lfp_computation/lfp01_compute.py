import h5py

from source.src2_utils.ut0_random_manager import np

def reconstruct_lfp_from_hdf(imem_hdf_path, electrode, out_path, sigma=None):
    """
    Reconstruct LFP from a single HDF5 file containing all cortical cell i_membrane currents.

    Args:
        imem_hdf_path (str): Path to HDF5 with structure: cells/{cell_id}/imem, meta
        electrode (MasmanidisEl): Electrode with site_topology
        out_dir (str): Directory to save the reconstructed LFP HDF5
        sigma (float, optional): Extracellular conductivity [S/m], defaults to electrode.sigma

    Returns:
        out_path (str): Path to saved HDF5 file containing LFP
    """
    sigma = sigma or electrode.sigma

    # Load electrode coordinates
    site_ids = list(electrode.site_topology.keys())
    electrode_xyz = np.array([electrode.site_topology[s] for s in site_ids])  # (n_sites,3)
    n_sites = electrode_xyz.shape[0]

    # Initialize storage for LFP
    lfp_dict = {}

    # Open main imem HDF5
    with h5py.File(imem_hdf_path, 'r') as f:
        tvec = f["time"][:]
        cells_grp = f["cells"]

        for cell_id in cells_grp.keys():
            cell_grp = cells_grp[cell_id]
            imem = cell_grp["imem"][:]   # (n_segments, n_time)
            meta = cell_grp["meta"][:]   # structured array with fields x_um, y_um, z_um, area_cm2, length_um, stype

            n_segments, n_time = imem.shape
            lfp_cell = np.zeros((n_sites, n_time), dtype=float)

            # Line-source approximation for each segment
            for seg_idx in range(n_segments):
                seg_xyz = np.array([meta[seg_idx]["x_um"], meta[seg_idx]["y_um"], meta[seg_idx]["z_um"]])
                r = np.linalg.norm(electrode_xyz - seg_xyz, axis=1)  # um
                r[r < 1.0] = 1.0  # prevent singularity
                r_m = r * 1e-6    # convert to meters
                # contribution: I / (4*pi*sigma*r)
                lfp_cell += imem[seg_idx, :] / (4 * np.pi * sigma * r_m[:, np.newaxis])

            lfp_dict[cell_id] = lfp_cell

    # Save to HDF5
    with h5py.File(out_path, 'w') as f_out:
        f_out.create_dataset("time", data=tvec)
        grp_lfp = f_out.create_group("lfp")
        for cell_id, lfp_cell in lfp_dict.items():
            grp_lfp.create_dataset(f"{cell_id}", data=lfp_cell)

        # Save electrode topology
        grp_el = f_out.create_group("electrode")
        grp_el.create_dataset("site_ids", data=np.array(site_ids, dtype='S'))
        grp_el.create_dataset("xyz_um", data=electrode_xyz)

    print(f"\t \t l01: LFP reconstructed and saved to {out_path}")
    return out_path

def compute_and_save_net_lfp(lfp_hdf_path, electrode_site_ids=None):
    """
    Compute net LFP by summing contributions from all cortical cells
    and save as new group `net_lfp` in the same HDF5 file.

    Args:
        lfp_hdf_path (str): Path to reconstructed LFP HDF file.
        electrode_site_ids (list of str, optional): List of electrode IDs.
            If None, they are inferred from the HDF.
    """
    with h5py.File(lfp_hdf_path, 'a') as f:
        lfp_grp = f["lfp"]
        time = f["time"][:]

        # Get electrode site IDs if not provided
        if electrode_site_ids is None:
            electrode_grp = f["electrode"]
            electrode_site_ids = [s.decode("utf-8") for s in electrode_grp["site_ids"][:]]
        n_sites = len(electrode_site_ids)
        n_time = len(time)

        # Initialize net LFP array
        net_lfp = np.zeros((n_sites, n_time), dtype=np.float64)

        # Iterate over all cells
        for cell_key in lfp_grp:
            data = lfp_grp[cell_key][:]
            mask = np.isnan(data)
            if np.any(mask):
                print(f"Warning: {cell_key} contains {mask.sum()} NaNs → replacing with 0")
                data[mask] = 0
            net_lfp += data

        # Create (or overwrite) net_lfp group
        if "net_lfp" in f:
            del f["net_lfp"]
        net_grp = f.create_group("net_lfp")
        net_grp.create_dataset("lfp", data=net_lfp)
        net_grp.create_dataset("time", data=time)
        dt = h5py.string_dtype(encoding='utf-8')
        net_grp.create_dataset("site_ids", data=np.array(electrode_site_ids, dtype=dt))

    print(f"\t \t l01: Net LFP computed and saved in group 'net_lfp' ({n_sites} sites, {n_time} time points).")