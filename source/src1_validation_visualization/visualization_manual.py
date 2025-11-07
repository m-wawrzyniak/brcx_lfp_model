
import pandas as pd
import matplotlib.pyplot as plt
import os

cx_v_csv = "/home/mateusz-wawrzyniak/PycharmProjects/brcx_lfp_model/data/weight_skeleton_10/recordings/cells/cx/cx_v.csv"
cx_pop_csv = "/home/mateusz-wawrzyniak/PycharmProjects/brcx_lfp_model/data/weight_skeleton_10/setup/cells/cx/cx04/cx04_pop_post3.csv"
save_dir = "/home/mateusz-wawrzyniak/PycharmProjects/brcx_lfp_model/data/weight_skeleton_10/visualizations/sim_dep/cells/cx/example_v"


def plot_cell_v(cell_ids,
                v_csv: str,
                pop_csv: str,
                save_dir: str):
    """
    Plot and save membrane potential traces for one or multiple cells.

    Parameters
    ----------
    cell_ids : int | list[int]
        Single cell ID or list of IDs to plot.
    v_csv : str
        Path to the CSV file with voltage recordings (time vs cell_id).
    pop_csv : str
        Path to the CSV file with cell metadata.
    save_dir : str
        Directory where plots will be saved. Created if it doesn't exist.
    """

    # Ensure save directory exists
    os.makedirs(save_dir, exist_ok=True)

    # Load voltage data
    df_v = pd.read_csv(v_csv)
    df_pop = pd.read_csv(pop_csv)

    # Handle single cell input
    if isinstance(cell_ids, int):
        cell_ids = [cell_ids]

    # Iterate through requested cells
    for cell_id in cell_ids:
        if str(cell_id) not in df_v.columns:
            print(f"⚠️ Cell ID {cell_id} not found in {v_csv}")
            continue

        # Try to get metadata
        if cell_id in df_pop["cell_id"].values:
            cell_info = df_pop.loc[df_pop["cell_id"] == cell_id].iloc[0]
            label = f"{cell_info['lay_m_type']} ({cell_info['e_type']}, {cell_info['m_type']})"
            desc = str(cell_info["desc"])
        else:
            label = f"Cell {cell_id}"
            desc = "unknown"

        # Plot
        plt.figure(figsize=(8, 4))
        plt.plot(df_v["time(ms)"], df_v[str(cell_id)], linewidth=1.2)
        plt.title(f"Membrane potential of cell {cell_id} [{desc}] during the simulation")
        plt.xlabel("Time [ms]")
        plt.ylabel("Membrane potential [mV]")
        plt.grid(True, alpha=0.3)
        plt.tight_layout()

        # Save
        filename = f"{cell_id}_{desc}_v.jpg".replace("/", "_")
        filepath = os.path.join(save_dir, filename)
        plt.savefig(filepath, dpi=300)
        plt.close()

        print(f"✅ Saved: {filepath}")

"""
plot_cell_v(cell_ids=[91, 92, 142, 16, 17, 178, 179],
            v_csv=cx_v_csv,
            pop_csv=cx_pop_csv,
            save_dir=save_dir)
"""

from source.src0_core.cr2_lfp_reconstruction.l02_lfp_computation import lfp02_visualize as lfp02
import config_templates.conf0_model_parameters as conf0
import config_templates.conf02_lfp_parameters as conf02

RECONSTRUCTED_LFP_PATH = "/home/mateusz-wawrzyniak/PycharmProjects/brcx_lfp_model/data/weight_skeleton_11/recordings/lfp/reconstructed/component_lfp.hdf"
LFP_NET_VIS_SAVEPATH = "/home/mateusz-wawrzyniak/PycharmProjects/brcx_lfp_model/data/weight_skeleton_11/visualizations/lfp/net_lfp/net_lfp.jpg"

stimulation_paradigm = conf0.WHISKER_STIMULATION_PARADIGMS[conf0.STIM_PARADIGM_TYPE][conf0.STIM_PARADIGM_SUBTYPE]
lfp02.plot_net_lfp(
    lfp_hdf_path=RECONSTRUCTED_LFP_PATH,
    save_path=LFP_NET_VIS_SAVEPATH,
    el_variant=conf02.ELECTRODE_VARIANT,
    el_z_offset=conf02.ELECTRODE_Z_OFFSET,
    crossspecies_corr=conf02.CROSS_SPECIES_SCALE,
    stim_paradigm=stimulation_paradigm,
    offset=0.01)
