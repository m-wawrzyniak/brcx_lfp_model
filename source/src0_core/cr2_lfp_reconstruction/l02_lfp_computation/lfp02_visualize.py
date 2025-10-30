import os, h5py
import matplotlib.pyplot as plt

from source.src2_utils.ut0_random_manager import np

def _plot_cell_lfp_with_morpho(cell_id, lfp_hdf_path, imem_hdf_path, electrode,
                               save_path=None, lfp_group='lfp', cell_group='cells'):
    """
    Plot LFP traces for a specific cell alongside 2D projection of its morphology + electrode sites.

    Left subplot: x-z projection of soma, dendrite, axon segments and electrode sites.
    Right subplot: LFP traces from all electrode sites for this cell.

    Args:
        cell_id (str or int): ID of the cell
        lfp_hdf_path (str): Path to HDF5 containing reconstructed LFPs
        imem_hdf_path (str): Path to HDF5 containing imem metadata for morphology
        electrode (MasmanidisEl): Electrode object with site_topology
        save_path (str or None): If provided, save the figure to this path
        lfp_group (str): HDF5 group containing LFP data
        cell_group (str): HDF5 group containing cell morphology data
    """
    cell_key = str(cell_id)

    # --- Load morphology metadata ---
    with h5py.File(imem_hdf_path, 'r') as f:
        meta = f[f"{cell_group}/{cell_key}/meta"][:]

    soma_mask = meta['stype'] == b'soma'
    dend_mask = (meta['stype'] == b'dend') | (meta['stype'] == b'dendrite')
    axon_mask = meta['stype'] == b'axon'

    soma_coords = np.stack([meta['x_um'][soma_mask], meta['z_um'][soma_mask]], axis=1)
    dend_coords = np.stack([meta['x_um'][dend_mask], meta['z_um'][dend_mask]], axis=1)
    axon_coords = np.stack([meta['x_um'][axon_mask], meta['z_um'][axon_mask]], axis=1)

    # --- Load electrode site positions ---
    site_ids = list(electrode.site_topology.keys())
    site_positions = np.array([electrode.site_topology[s] for s in site_ids])  # (n_sites, 3)
    site_xz = site_positions[:, [0, 2]]  # project to x–z plane

    # --- Load LFP data ---
    with h5py.File(lfp_hdf_path, 'r') as f:
        tvec = f["time"][:]
        lfp_grp = f[lfp_group]
        if cell_key not in lfp_grp:
            raise ValueError(f"Cell {cell_key} not found in LFP file.")
        lfp_data = lfp_grp[cell_key][:]  # shape: (n_sites, n_time)

    # --- Plot ---
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Left: morphology + electrode sites (xz-plane)
    ax0 = axes[0]
    if dend_coords.size > 0:
        ax0.scatter(dend_coords[:,0], dend_coords[:,1], c='blue', s=5, label='Dendrite')
    if axon_coords.size > 0:
        ax0.scatter(axon_coords[:,0], axon_coords[:,1], c='red', s=5, label='Axon')
    if soma_coords.size > 0:
        ax0.scatter(soma_coords[:,0], soma_coords[:,1], c='green', s=15, label='Soma')

    ax0.scatter(site_xz[:,0], site_xz[:,1], c='black', s=40, marker='x', label='Electrode sites')
    x_off, y_off = 10, 10
    for i, pos in enumerate(site_xz):
        ax0.text(pos[0]+x_off, pos[1]+y_off, site_ids[i], fontsize=8, color='black')
    ax0.grid(
        True,  # turn grid on
        which='both',  # 'major', 'minor', or 'both'
        linestyle='--',  # dashed lines
        linewidth=0.5,  # line thickness
        color='gray',  # line color
        alpha=0.7  # transparency
    )
    ax0.set_xlabel("x [µm]")
    ax0.set_ylabel("z [µm]")
    ax0.set_title(f"Cortical cell {cell_id} morphology (xz-plane)")
    ax0.legend(fontsize=8)

    # Right: LFP traces
    ax1 = axes[1]
    for i, sid in enumerate(site_ids):
        ax1.plot(tvec[5:], lfp_data[i, 5:], label=sid)
    ax1.set_xlabel("Time [ms]")
    ax1.set_ylabel("LFP [a.u.]")
    ax1.set_title(f"LFP traces for cortical cell {cell_id}")
    ax1.legend(fontsize=7, loc='upper right')

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300)
        plt.close(fig)
    else:
        plt.show()

def export_all_cells_lfp_plots(lfp_hdf_path, imem_hdf_path, electrode, out_dir,
                               lfp_group='lfp', cell_group='cells'):
    """
    Generate and save morphology+LFP plots for all cells.

    Args:
        lfp_hdf_path (str): Path to LFP HDF5 file
        imem_hdf_path (str): Path to morphology HDF5 file
        electrode (MasmanidisEl): Electrode object with site_topology
        out_dir (str): Directory where figures will be saved
        lfp_group (str): Group in HDF containing LFP data
        cell_group (str): Group in HDF containing cell morphology data
    """

    with h5py.File(lfp_hdf_path, 'r') as f:
        lfp_grp = f[lfp_group]
        cell_ids = list(lfp_grp.keys())

    n_cells = len(cell_ids)

    for i, cid in enumerate(cell_ids, 1):
        save_path = os.path.join(out_dir, f"cell_{cid}.png")
        _plot_cell_lfp_with_morpho(cid, lfp_hdf_path, imem_hdf_path, electrode,
                                   save_path=save_path,
                                   lfp_group=lfp_group, cell_group=cell_group)

        percent = (i / n_cells) * 100
        print(f"\t \t l02: LFP component plotted: {percent:.1f}%")


def plot_net_lfp(lfp_hdf_path, offset, save_path=None, stim_paradigm=None):
    """
    Plot net LFP (summed over all cortical cells) with offset traces per site.

    Args:
        lfp_hdf_path (str): Path to reconstructed LFP HDF5 file
        offset (float): Vertical offset between traces
        save_path (str or None): If provided, save figure instead of showing
        stim_paradigm (dict or None): Optional, e.g., {'phase1': ('wk', 50), 'phase2': ('str', 100)}
    """
    STIM_COLORS = {
        'r': 'white',
        'wk': '#fff7a0',  # yellowish
        'str': '#ff9999'  # reddish
    }

    with h5py.File(lfp_hdf_path, "r") as f:
        time = f["net_lfp/time"][:]
        lfp = f["net_lfp/lfp"][:]  # shape (n_sites, n_time)
        site_ids = f["net_lfp/site_ids"][:]

    # --- Drop first 10 samples ---
    time = time[10:]
    lfp = lfp[:, 10:]

    n_sites = lfp.shape[0]

    plt.figure(figsize=(12, 7))

    # --- Plot stimulus windows ---
    if stim_paradigm:
        x_start = time[0]
        for phase, (phase_type, duration) in stim_paradigm.items():
            if duration == 0:
                continue
            x_end = x_start + duration
            plt.axvspan(x_start, x_end, color=STIM_COLORS.get(phase_type, 'gray'), alpha=0.3)
            x_start = x_end

    # --- Colormap for traces ---
    cmap = plt.get_cmap("tab20")

    # --- Plot LFP traces with offset and colors ---
    for i in range(n_sites):
        trace = lfp[i, :] + i * offset
        sid = site_ids[i].decode() if isinstance(site_ids[i], bytes) else str(site_ids[i])
        color = cmap(i % cmap.N)
        plt.plot(time, trace, color=color, label=sid)

    # --- Vertical dashed lines every 100 ms ---
    t_min, t_max = time[0], time[-1]
    for tline in range(int(t_min // 100) * 100, int(t_max) + 100, 100):
        plt.axvline(tline, linestyle='--', color='gray', linewidth=0.5)

    # --- Formatting ---
    plt.xlabel("Time [ms]")
    plt.ylabel("LFP with spatial offset")
    plt.title("Net LFP at electrode sites")
    plt.ylim(-0.005, n_sites * offset + 0.002)
    plt.yticks([])
    plt.legend(loc="upper right", fontsize=8, title="Electrodes")
    plt.tight_layout()

    # --- Save or show ---
    if save_path:
        plt.savefig(save_path, dpi=300)
        plt.close()
        print('\t \t l02: Net LFP plotted.')
    else:
        plt.show()