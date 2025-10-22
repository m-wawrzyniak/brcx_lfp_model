import os
import pandas as pd

# TODO: mess

def _cxcx_conn_epsp(pre_id, post_id,
                    synapses_params_csv, cell_pop_csv,
                    output_path):
    """
    1. Initializes pre_cell and post_cell based on their me-types in synapse *.json file.
    2. Creates all synapses described in the *.json between these two cells.
    3. Artificially stimulates pre-cell with 30 Hz pulse train (9 pulses) and 1 recovery AP.
    4. During this, records EPSP and saves it as *.csv.
    :return:
    """


    # --- Load inputs ---
    synsdf = pd.read_csv(synapses_params_csv)
    cell_temps = b01.__main__()

    cell_df = pd.read_csv(cell_pop_csv)

    # --- Filter synapses ---
    conn_df = synsdf[(synsdf['pre_id'] == pre_id) & (synsdf['post_id'] == post_id)]
    if conn_df.empty:
        raise ValueError(f"No synapses found between {pre_id} and {post_id}")
    #print(conn_df)
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

    pre_cell = CxCell(pre_me_type, pre_id, cell_temp=cell_temps[pre_me_abb],
                      offset = pre_offset, rotation= pre_rotation)
    post_cell = CxCell(post_me_type, post_id, cell_temp=cell_temps[post_me_abb],
                       offset= post_offset, rotation= post_rotation)

    # --- Create synapses ---
    synapses = {}
    for _, row in conn_df.iterrows():
        syn_id = row['syn_id']
        syn_type = row['syn_type']
        post_loc_rel = (row['post_loc_x'], row['post_loc_y'], row['post_loc_z'])
        syn_param_dict = dict()
        syn_param_dict['late_comp_ratio'] = row['late_comp_ratio']

        for param in ['u', 'd', 'f', 'risetime', 'decay', 'gsyn']:
            syn_param_dict[f'{param}_mean'] = row[f'{param}_val']
            syn_param_dict[f'{param}_std'] = 0.0  # deterministic
        cx_syn = b03._create_cxcx_synapse(pre_cell, post_cell, post_loc_rel,
                                                syn_type, syn_param_dict, deterministic=True, syn_id=syn_id)
        #print(f"{cx_syn.syn_id} Rstate:{cx_syn.hsyn.Rstate}")
        #print(f"{cx_syn.syn_id} tsyn_dac:{cx_syn.hsyn.tsyn_fac}")
        #print(f"{cx_syn.syn_id} u:{cx_syn.hsyn.u}")
        synapses[cx_syn.syn_id] = cx_syn

        # --- Stimulate pre_cell using StimSynapse ---
    stim_times = [i * (1000 / 30) for i in range(9)] + [600.0]  # 30Hz + recovery spike
    stim_syn = StimSynapse()
    stim_syn.set_stimvec(stim_times)
    stim_syn.connect_synapse(pre_cell)

    # --- Record EPSP ---
    cells = {
        pre_id : pre_cell,
        post_id: post_cell
    }

    pprint(cells)

    t = h.Vector()
    t.record(h._ref_t)

    v_recordings = r01.record_cx_cells_v(cells)
    i_recordings = r01.record_synapses_currents(synapses)

    print('Running simulation')
    h.finitialize(-70)
    h.continuerun(800.0)

    pair_dir = os.path.join(output_path, f"epsp-{pre_me_type}-{post_me_type}")
    os.makedirs(pair_dir, exist_ok=True)

    voltage_file = os.path.join(pair_dir, f"cellv-{pre_me_type}-{post_me_type}.csv")
    syn_currents_file = os.path.join(pair_dir, f"syni-{pre_me_type}-{post_me_type}.csv")

    r01.save_cell_v_csv(v_recordings, t, voltage_file)
    r01.save_synapses_currents_csv(i_recordings, t, syn_currents_file)

def check_all_connection_types(synapses_csv_path: str, output_base_dir: str,
                               cell_pop_csv):
    """
    For each unique pre_me_type → post_me_type combination in the synapse CSV,
    simulate one representative EPSP and save results.

    Args:
        synapses_csv_path (str): Path to the synapse parameter CSV file.
        output_base_dir (str): Directory to store EPSP results.
    """
    # Load synapse parameter file
    df = pd.read_csv(synapses_csv_path)

    # Drop NaNs just in case
    df = df.dropna(subset=["pre_id", "post_id", "pre_me_type", "post_me_type"])

    # Extract all unique me_type pairs
    unique_pairs = df[["pre_me_type", "post_me_type"]].drop_duplicates()

    for _, row in unique_pairs.iterrows():
        #print(row)
        pre_type = row["pre_me_type"]
        post_type = row["post_me_type"]

        # Get first matching synapse for this type pair
        match = df[
            (df["pre_me_type"] == pre_type) &
            (df["post_me_type"] == post_type)
        ].iloc[0]
        #print(match)

        pre_id = int(match["pre_id"])
        post_id = int(match["post_id"])

        # Prepare EPSP output dir and plot path
        epsp_dirname = f"epsp-{pre_type}-{post_type}"
        epsp_dir = os.path.join(output_base_dir, epsp_dirname)

        os.makedirs(epsp_dir, exist_ok=True)

        print(f"[SIMULATING] {pre_type} → {post_type} (IDs: {pre_id} → {post_id})")
        #print(pre_id, post_id, synapses_csv_path, cell_pop_csv, output_base_dir, epsp_dir)
        _cxcx_conn_epsp(pre_id, post_id, synapses_csv_path, cell_pop_csv, output_base_dir)
        plot_epsp(epsp_dir)

    print("[DONE] All unique EPSPs simulated.")