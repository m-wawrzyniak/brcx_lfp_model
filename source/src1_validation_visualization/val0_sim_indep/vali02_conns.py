import os
from pathlib import Path

import pandas as pd
from neuron import h
import matplotlib.pyplot as plt

from source.src0_core.cr1_simulation_run.s02_cells_init.CxCell import CxCell
from source.src0_core.cr1_simulation_run.s02_cells_init.VPMCell import VPMCell

from source.src0_core.cr1_simulation_run.s03_conn_init import inits01_cxcx as inits01
from source.src0_core.cr1_simulation_run.s03_conn_init import inits02_tccx as inits02
from source.src0_core.cr1_simulation_run.s05_data_recording_saving import (rec01_setup_record as rec01,
                                                                           rec02_record_save as rec02)

from source.src2_utils.ut0_random_manager import np


def _cxcx_conn_epsp(pre_id, post_id,
                    synapses_params_csv, cell_pop_csv,
                    cell_templates, output_dir):
    """
    1. Initializes pre_cell and post_cell based on their me-types in synapse *.json file.
    2. Creates all synapses described in the *.json between these two cells.
    3. Artificially stimulates pre-cell with 30 Hz pulse train (9 pulses) and 1 recovery AP.
    4. During this, records EPSP and saves it as *.csv.
    """

    # --- Load inputs ---
    synsdf = pd.read_csv(synapses_params_csv)
    cell_df = pd.read_csv(cell_pop_csv)

    # --- Filter synapses ---
    conn_df = synsdf[(synsdf['pre_id'] == pre_id) & (synsdf['post_id'] == post_id)]
    if conn_df.empty:
        raise ValueError(f"No synapses found between {pre_id} and {post_id}")

    # --- Cell instantiation ---
    pre_me_type = conn_df.iloc[0]['pre_me_type']
    post_me_type = conn_df.iloc[0]['post_me_type']
    pre_me_abb = "_".join(pre_me_type.split('_')[:2])
    post_me_abb = "_".join(post_me_type.split('_')[:2])

    def get_pos_and_rot(cell_id):
        row = cell_df[cell_df['cell_id'] == cell_id]
        if row.empty:
            raise ValueError(f"Cell ID {cell_id} not found in {cell_pop_csv}")
        r = row.iloc[0]
        return (r['x'], r['y'], r['z']), (np.pi/2, 0.0, r['rot_ang'])

    pre_offset, pre_rotation = get_pos_and_rot(pre_id)
    post_offset, post_rotation = get_pos_and_rot(post_id)
    pre_cell = CxCell(
        cell_name=pre_me_type, cell_id=pre_id, cell_temp=cell_templates[pre_me_abb],
        offset=pre_offset, rotation=pre_rotation)
    pre_cell.populate_segs_meta_for_synapses()

    post_cell = CxCell(
        cell_name=post_me_type, cell_id=post_id,
        cell_temp=cell_templates[post_me_abb], offset=post_offset, rotation=post_rotation)
    post_cell.populate_segs_meta_for_synapses()

    # --- Create synapses ---
    synapses = {}
    for _, row in conn_df.iterrows():
        syn_id = row['syn_id']
        syn_type = row['syn_type']
        post_loc_rel = (row['post_loc_x'], row['post_loc_y'], row['post_loc_z'])
        syn_param_dict = {'late_comp_ratio': row['late_comp_ratio']}
        for param in ['u', 'd', 'f', 'risetime', 'decay', 'gsyn']:
            syn_param_dict[f'{param}_mean'] = row[f'{param}_val']
            syn_param_dict[f'{param}_std'] = 0.0
        cx_syn = inits01._create_cxcx_synapse(
            pre_cell, post_cell, post_loc_rel, syn_type, syn_param_dict,
            deterministic=True, syn_id=syn_id
        )
        synapses[cx_syn.syn_id] = cx_syn

    # --- Stimulate pre_cell directly (ad hoc VecStim) ---
    stim_times = [i * (1000 / 30) for i in range(9)] + [600.0]  # 30Hz + recovery spike
    spike_vec = h.Vector(stim_times)
    vecstim = h.VecStim()
    vecstim.play(spike_vec)

    # Create ExpSyn on pre_cell soma
    stim_syn = h.ExpSyn(pre_cell.h_cell.soma[0](0.5))

    # Set synaptic time constant (ms)
    stim_syn.tau = 1
    stim_nc = h.NetCon(vecstim, stim_syn)
    stim_nc.delay = 0
    stim_nc.weight[0] = 0.5

    # --- Record EPSP and synaptic currents ---
    cells = {pre_id: pre_cell, post_id: post_cell}

    t = h.Vector()
    t.record(h._ref_t)

    v_recordings = rec01.record_soma_v(cells)
    i_recordings = rec01.record_synapses_currents(synapses)

    print('Running simulation...')
    h.finitialize(-70)
    h.continuerun(800.0)

    voltage_file = os.path.join(output_dir, f"cellv_{pre_me_type}_{post_me_type}.csv")
    syn_currents_file = os.path.join(output_dir, f"syni_{pre_me_type}_{post_me_type}.csv")

    rec02.save_cell_v_csv(v_recordings, t, voltage_file)
    rec02.save_synapses_currents_csv(i_recordings, t, syn_currents_file)

def cxcx_conn_epsp_check(synapses_csv_path: str, cell_pop_csv,
                         cell_templates, save_dir: str,
                         ):
    """
    For each unique pre_me_type → post_me_type combination in the synapse CSV,
    simulate one representative EPSP and save results.

    Args:
        synapses_csv_path (str): Path to the synapse parameter CSV file.
        save_dir (str): Directory to store EPSP results.
    """
    # Load synapse parameter file
    df = pd.read_csv(synapses_csv_path)

    # Drop NaNs just in case
    df = df.dropna(subset=["pre_id", "post_id", "pre_me_type", "post_me_type"])

    # Extract all unique me_type pairs
    unique_pairs = df[["pre_me_type", "post_me_type"]].drop_duplicates()

    for _, row in unique_pairs.iterrows():
        pre_type = row["pre_me_type"]
        post_type = row["post_me_type"]

        # Get first matching synapse for this type pair
        match = df[
            (df["pre_me_type"] == pre_type) &
            (df["post_me_type"] == post_type)
        ].iloc[0]

        pre_id = int(match["pre_id"])
        post_id = int(match["post_id"])

        # Prepare EPSP output dir and plot path
        epsp_dirname = f"epsp_{pre_type}_{post_type}"
        epsp_dir = os.path.join(save_dir, epsp_dirname)
        os.makedirs(epsp_dir, exist_ok=True)

        print(f"[SIMULATING] {pre_type} → {post_type} (IDs: {pre_id} → {post_id})")
        _cxcx_conn_epsp(pre_id, post_id, synapses_csv_path, cell_pop_csv,
                        cell_templates=cell_templates, output_dir=epsp_dir)
        #plot_epsp(epsp_dir)

    print("[DONE] All unique EPSPs simulated.")


def conn_epsp_plot(me_pair_epsp_dir: str):
    """
    Scan all subdirectories in root_epsp_dir. Each subdir corresponds to a pre→post pair.
    Look for cellv_*.csv and syni_*.csv files and plot EPSPs + synaptic currents.

    Args:
        me_pair_epsp_dir (str): Root directory containing pair subdirectories.
    """
    me_pair_epsp_dir = os.path.abspath(me_pair_epsp_dir)
    for pair_name in os.listdir(me_pair_epsp_dir):
        pair_dir = os.path.join(me_pair_epsp_dir, pair_name)
        if not os.path.isdir(pair_dir):
            continue

        # Find CSVs
        cellv_files = [f for f in os.listdir(pair_dir) if f.startswith("cellv")]
        syni_files = [f for f in os.listdir(pair_dir) if f.startswith("syni")]
        if not cellv_files or not syni_files:
            print(f"[SKIP] {pair_name}: missing CSVs")
            continue

        cellv_file = os.path.join(pair_dir, cellv_files[0])
        syni_file = os.path.join(pair_dir, syni_files[0])

        # Load CSVs
        cellv_df = pd.read_csv(cellv_file)
        syni_df = pd.read_csv(syni_file)

        # Extract columns
        time = cellv_df.iloc[:, 0]
        pre_v = cellv_df.iloc[:, 1]
        post_v = cellv_df.iloc[:, 2]
        pre_id = cellv_df.columns[1]
        post_id = cellv_df.columns[2]

        # Plot
        fig, axes = plt.subplots(3, 1, figsize=(8, 8), sharex=True)
        fig.suptitle(pair_name.replace("_", " "), fontsize=12)

        axes[0].plot(time, pre_v, color='tab:blue')
        axes[0].set_ylabel(f"Pre-cell V (mV)\n[{pre_id}]")
        axes[0].grid(True, alpha=0.3)

        axes[1].plot(time, post_v, color='tab:orange')
        axes[1].set_ylabel(f"Post-cell V (mV)\n[{post_id}]")
        axes[1].grid(True, alpha=0.3)

        for col in syni_df.columns[1:]:
            axes[2].plot(syni_df.iloc[:, 0], syni_df[col], label=col)
        axes[2].set_xlabel("Time (ms)")
        axes[2].set_ylabel("Synaptic I (nA)")
        axes[2].legend(fontsize=8, ncol=2, loc='upper right', frameon=False)
        axes[2].grid(True, alpha=0.3)

        plt.tight_layout(rect=[0, 0, 1, 0.96])

        # Save
        out_path = os.path.join(pair_dir, "epsp_plot.png")
        plt.savefig(out_path, dpi=200)
        plt.close(fig)
        print(f"[SAVED] {out_path}")


def _tccx_conn_epsp(pre_id, post_id,
                    tccx_synapses_params_csv, cell_pop_csv,
                    cell_templates, output_dir):
    """
    1. Initializes pre_cell (VPMCell) and post_cell (CxCell) based on their types.
    2. Creates all TC->CX synapses described in tccx CSV.
    3. Artificially stimulates pre-cell with 30 Hz spike train (9 pulses) + 1 recovery AP.
    4. Records EPSPs in post-cell soma and synaptic currents.
    5. Saves both to CSV files.
    """

    # --- Load inputs ---
    synsdf = pd.read_csv(tccx_synapses_params_csv)
    cell_df = pd.read_csv(cell_pop_csv)

    # --- Filter synapses ---
    conn_df = synsdf[(synsdf['pre_id'] == pre_id) & (synsdf['post_id'] == post_id)]
    if conn_df.empty:
        raise ValueError(f"No synapses found between {pre_id} and {post_id}")

    # --- Instantiate pre-cell (VPM) ---
    pre_cell = VPMCell(cell_id=pre_id)

    # --- Instantiate post-cell (CxCell) ---
    post_row = cell_df[cell_df['cell_id'] == post_id]
    if post_row.empty:
        raise ValueError(f"Post cell {post_id} not found in {cell_pop_csv}")
    r = post_row.iloc[0]
    post_offset = (r['x'], r['y'], r['z'])
    post_rotation = (np.pi/2, 0.0, r['rot_ang'])

    post_me_type = conn_df.iloc[0]['post_me_type']
    post_me_abb = "_".join(post_me_type.split('_')[:2])

    post_cell = CxCell(
        cell_name=post_me_type, cell_id=post_id,
        cell_temp=cell_templates[post_me_abb],
        offset=post_offset, rotation=post_rotation
    )
    post_cell.populate_segs_meta_for_synapses()

    # --- Create TC->CX synapses ---
    synapses = {}
    for _, row in conn_df.iterrows():
        syn = inits02._create_tccx_synapse(
            pre_vpm_cell=pre_cell,
            post_cell=post_cell,
            post_loc_rel=(row['post_loc_x'], row['post_loc_y'], row['post_loc_z']),
            syn_type=row['syn_type']
        )
        synapses[syn.syn_id] = syn

    # --- Artificial stimulation using VecStim (30 Hz + recovery spike) ---
    stim_times = [i * (1000 / 30) for i in range(9)] + [600.0]  # ms
    spike_vec = h.Vector(stim_times)
    vecstim = h.VecStim()
    vecstim.play(spike_vec)

    # Create ExpSyn on pre_cell soma
    stim_syn = h.ExpSyn(pre_cell.soma(0.5))

    # Set synaptic time constant (ms)
    stim_syn.tau = 1
    stim_nc = h.NetCon(vecstim, stim_syn)
    stim_nc.delay = 0
    stim_nc.weight[0] = 0.5

    # --- Set up recording ---
    t_vec = h.Vector()
    t_vec.record(h._ref_t)

    v_recordings = rec01.record_soma_v({"VPM":pre_cell, post_id: post_cell})
    i_recordings = rec01.record_synapses_currents(synapses)


    # --- Run simulation ---
    print("Running simulation...")
    h.finitialize(-70)
    h.continuerun(stim_times[-1] + 200)  # simulate slightly beyond last spike

    # --- Save results ---
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    voltage_file = output_dir / f"cellv_VPM_{post_me_type}.csv"
    syn_currents_file = output_dir / f"syni_VPM_{post_me_type}.csv"

    rec02.save_cell_v_csv(v_recordings, t_vec, voltage_file)
    rec02.save_synapses_currents_csv(i_recordings, t_vec, syn_currents_file)

    print(f"EPSPs saved to {voltage_file}")
    print(f"Synaptic currents saved to {syn_currents_file}")


def tccx_conn_epsp_check(tccx_synapses_csv_path: str, cell_pop_csv,
                         cell_templates, save_dir: str):
    """
    For each unique VPM → CX me_type combination in the tccx CSV,
    simulate one representative EPSP and save results.

    Args:
        tccx_synapses_csv_path (str): Path to the tccx synapse CSV.
        cell_pop_csv (str): Path to cortical cell population CSV.
        cell_templates (dict): Dictionary of cell templates.
        save_dir (str): Directory to store EPSP results.
    """
    # Load synapse CSV
    df = pd.read_csv(tccx_synapses_csv_path)
    df = df.dropna(subset=["pre_id", "post_id", "pre_me_type", "post_me_type"])

    # Extract unique VPM → CX pairs
    unique_pairs = df[["pre_me_type", "post_me_type"]].drop_duplicates()

    for _, row in unique_pairs.iterrows():
        pre_type = row["pre_me_type"]
        post_type = row["post_me_type"]

        # Take first matching synapse for this type pair
        match = df[
            (df["pre_me_type"] == pre_type) &
            (df["post_me_type"] == post_type)
        ].iloc[0]

        pre_id = match["pre_id"]
        post_id = match["post_id"]

        # Prepare EPSP output directory
        epsp_dirname = f"epsp_{pre_type}_{post_type}"
        epsp_dir = os.path.join(save_dir, epsp_dirname)
        os.makedirs(epsp_dir, exist_ok=True)

        print(f"[SIMULATING] {pre_type} → {post_type} (IDs: {pre_id} → {post_id})")
        _tccx_conn_epsp(pre_id, post_id,
                        tccx_synapses_params_csv=tccx_synapses_csv_path,
                        cell_pop_csv=cell_pop_csv,
                        cell_templates=cell_templates,
                        output_dir=epsp_dir)

    print("[DONE] All unique VPM→CX EPSPs simulated.")
