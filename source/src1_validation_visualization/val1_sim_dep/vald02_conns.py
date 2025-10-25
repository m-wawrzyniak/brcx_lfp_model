from pathlib import Path

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import json
from collections import Counter
from itertools import product

from source.src2_utils.ut0_random_manager import np

def conn_matrix_emp(csv_file, save_path, layer_comp_params, unique_connections=False):
    # Read CSV
    df = pd.read_csv(csv_file)

    # If counting unique connections, drop duplicates based on pre_id + post_id
    if unique_connections:
        df = df.drop_duplicates(subset=['pre_id', 'post_id'])

    # Generate all full cell types
    cell_types = [f"{layer}_{ctype}" for layer, types in layer_comp_params.items() for ctype in types.keys()]

    # Initialize matrices
    conn_matrix_ratio = pd.DataFrame(0.0, index=cell_types, columns=cell_types)
    conn_matrix_count = pd.DataFrame(0, index=cell_types, columns=cell_types, dtype=int)

    # Compute counts and ratios
    for pre_type in cell_types:
        total_possible = len(df[df['pre_me_type'] == pre_type])
        for post_type in cell_types:
            actual = len(df[(df['pre_me_type'] == pre_type) & (df['post_me_type'] == post_type)])
            conn_matrix_count.loc[pre_type, post_type] = actual
            conn_matrix_ratio.loc[pre_type, post_type] = 100 * actual / total_possible if total_possible > 0 else 0.0

    # Create combined labels for heatmap: count \n ratio%
    labels = conn_matrix_count.astype(str) + "\n" + conn_matrix_ratio.round(2).astype(str) + "%"

    # Plot heatmap
    plt.figure(figsize=(12, 10))
    ax = sns.heatmap(
        conn_matrix_ratio,
        annot=labels,
        fmt="",
        cmap="viridis",
        cbar_kws={'label': 'Connection ratio (%)'},
        linewidths=0.5,         # adds grid lines between cells
        linecolor='black',      # darker borders
        square=True             # make cells square
    )

    # Titles and labels
    if unique_connections:
        title = "Cortico-cortical connectivity matrix: connection count and row-wise ratio (%)"
    else:
        title = "Cortico-cortical connectivity matrix: synapse count and row-wise ratio (%)"
    ax.set_title(title, fontsize=13, pad=20)
    ax.set_xlabel("Post-synaptic cell type", fontsize=12, labelpad=10)
    ax.set_ylabel("Pre-synaptic cell type", fontsize=12, labelpad=10)

    # Rotate bottom tick labels for readability
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)

    # Tight layout to keep everything visible and centered
    plt.tight_layout()

    # Save or show
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.close()
    else:
        plt.show()


def conn_matrix_theo(json_file, layer_comp_params, save_path=None):
    # Load JSON
    with open(json_file, "r") as f:
        data = json.load(f)

    # Generate all full cell types (sorted for consistent axis order)
    cell_types = [f"{layer}_{'_'.join(ctype.split('_')[:-1])}" for layer, types in layer_comp_params.items() for ctype in types.keys()]
    cell_types = sorted(cell_types)

    # Initialize matrices
    conn_matrix_ratio = pd.DataFrame(0.0, index=cell_types, columns=cell_types)
    conn_matrix_count = pd.DataFrame(0, index=cell_types, columns=cell_types, dtype=int)

    # Fill in values from JSON
    for pair, values in data.items():
        pre, post = pair.split(":")
        if pre in cell_types and post in cell_types:
            conn_matrix_count.loc[pre, post] = values.get("total_synapse_count", 0)
            conn_matrix_ratio.loc[pre, post] = 100 * values.get("connection_probability", 0.0)  # Convert to %

    # Create combined labels: count \n ratio%
    labels = conn_matrix_count.astype(str) + "\n" + conn_matrix_ratio.round(2).astype(str) + "%"

    # Plot heatmap
    plt.figure(figsize=(12, 10))
    ax = sns.heatmap(
        conn_matrix_ratio,
        annot=labels,
        fmt="",
        cmap="viridis",
        cbar_kws={'label': 'Connection probability (%)'},
        linewidths=0.5,
        linecolor='black',
        square=True
    )

    # Titles and labels
    title = "Target cortico-cortical connectivity matrix: synapses count and row-wise ratio (%)"
    ax.set_title(title, fontsize=13, pad=20)
    ax.set_xlabel("Post-synaptic cell type", fontsize=12, labelpad=10)
    ax.set_ylabel("Pre-synaptic cell type", fontsize=12, labelpad=10)

    # Rotate bottom tick labels for readability
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)

    # Center and fit layout
    plt.tight_layout()

    # Save or show
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.close()
    else:
        plt.show()
        
def conn_matrix_emp_tc(tccx_syns: str,
                       layer_comp_params,
                       unique_connections: bool = False,
                       savefig_path=None):
    # --- Load data ---
    df = pd.read_csv(tccx_syns)

    # Only VPM presynaptic entries
    df = df[df['pre_me_type'].str.startswith("VPM")]
    if df.empty:
        raise ValueError("No VPM synapses found in the given file.")

    # If counting unique connections, drop duplicates
    if unique_connections:
        df = df.drop_duplicates(subset=['pre_id', 'post_id'])

    # --- Build list of all possible postsynaptic me-types ---
    all_cell_types = [
        f"{layer}_{ctype}" for layer, ctypes in layer_comp_params.items() for ctype in ctypes.keys()
    ]
    all_cell_types = sorted(all_cell_types)

    # --- Count connections to each type ---
    counts = df['post_me_type'].value_counts()
    counts = counts.reindex(all_cell_types, fill_value=0)  # include missing types
    total = counts.sum()

    # Convert to percentage ratios
    ratios = (100 * counts / total) if total > 0 else counts.astype(float)

    # --- Build one-row matrices ---
    matrix_count = pd.DataFrame([counts], index=["VPM"], columns=all_cell_types).astype(int)
    matrix_ratio = pd.DataFrame([ratios], index=["VPM"], columns=all_cell_types).astype(float)

    # Labels: count \n ratio%
    labels = matrix_count.astype(str) + "\n" + matrix_ratio.round(2).astype(str) + "%"

    # --- Plot ---
    plt.figure(figsize=(len(all_cell_types) * 1.1, 3.2))
    ax = sns.heatmap(
        matrix_ratio,
        annot=labels,
        fmt="",
        cmap="viridis",
        cbar_kws={'label': 'Connection ratio (%)'},
        linewidths=0.5,
        linecolor='black',
        square=True
    )

    # Titles and labels
    title = "Thalamocortical connectivity matrix: synapses and row-wise ratio (%)"
    ax.set_title(title, fontsize=13, pad=20)
    ax.set_xlabel("Post-synaptic cell type", fontsize=12, labelpad=10)
    ax.set_ylabel("Pre-synaptic cell type", fontsize=12, labelpad=10)

    # Rotate and align tick labels
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)

    # Center and fit layout
    plt.tight_layout()

    # Save or show
    if savefig_path:
        plt.savefig(savefig_path, dpi=300, bbox_inches='tight')
        plt.close()
    else:
        plt.show()

    return matrix_count, matrix_ratio

def interbouton_int_histogram(cx_pop_csv:Path, bins:int = 30, save_path= None):
    """
    Creates a histogram of b_int for the cells in examined_pop_csv. Saves *.png.

    Args:
        cx_pop_csv (str): Path to csv with cortical cells with already calculated b_int.
        bins (int): Number of bins in the histogram.
    """

    df = pd.read_csv(cx_pop_csv)
    if "b_int" not in df.columns:
        raise KeyError("'interb_int' column not found in the CSV.")

    data = df["b_int"].replace([np.inf, -np.inf], np.nan).dropna()

    # map filenames to descriptive titles
    title_map = {
        "intb_dist_pre.jpg": "pre-pruning",
        "intb_dist_post1.jpg": "after General Pruning (01)",
        "intb_dist_post2.jpg": "after Multisynapse Pruning (02)",
        "intb_dist_post3.jpg": "after Plasticity-Reserve Pruning (03)"
    }

    # get the descriptive title
    file_name = cx_pop_csv.name  # or Path(cx_pop_csv).name if needed
    plot_title = title_map.get(file_name, file_name)  # fallback to filename if not in map

    plt.figure(figsize=(8, 5))
    plt.hist(data, bins=bins, edgecolor='black')
    plt.grid(True)
    plt.xlim(0, 100)  # limit x-axis to 0-100
    plt.xlabel("Interbouton interval (um)")
    plt.ylabel("No. cells")
    plt.title(f"Distribution of interbouton interval ({plot_title})")
    plt.grid(True)
    plt.tight_layout()

    plt.savefig(save_path)

def syn_per_conn_histogram(synapses_json: Path, savepath):
    """
    Creates a histogram of synapses per connection for all connections. Saves it as PNG.
    """
    # Map filenames to descriptive titles
    title_map = {
        "cxcx01_appositions.json": "pre-pruning",
        "prune01_cxcx.json": "after General Pruning (01)",
        "prune02_cxcx.json": "after Multisynapse Pruning (02)",
        "prune03_cxcx.json": "after Plasticity-Reserve Pruning (03)"
    }

    main_title = f"Corticocortical synapses per connection distribution ({title_map.get(synapses_json.name, synapses_json.name)})"

    # Load JSON
    with open(synapses_json, "r") as f:
        data = json.load(f)

    # Count synapses per pre-post pair
    synapse_pairs = []
    for post_id_str, syns in data.items():
        post_id = int(post_id_str)
        for syn in syns:
            pre_id = syn["pre_id"]
            synapse_pairs.append((pre_id, post_id))
    pair_counts = Counter(synapse_pairs)

    synapse_nums = list(pair_counts.values())

    # Plot histogram
    plt.figure(figsize=(8, 5))
    bins = np.arange(0, 21, 1)  # bin edges 0–20
    plt.hist(synapse_nums, bins=bins, edgecolor='black', density=True, align='left')
    plt.xlabel("Synapses per connection")
    plt.ylabel("Probability (PMF)")
    plt.title(main_title)

    # Axis limits and ticks
    plt.xlim(0, 20)
    plt.xticks(np.arange(0, 21, 2.5))
    plt.ylim(0, 0.15)
    plt.yticks(np.arange(0, 0.151, 0.025))
    plt.grid(True, linestyle='--', alpha=0.5)

    plt.tight_layout()
    plt.savefig(savepath, dpi=300)
    plt.close()

def _conn_prob_with_distance(cx_cells_csv, cxcx_synapses_json, bin_size=10, max_dist=None):
    """
    Compute connection probability vs intersomatic distance.
    """
    # --- Load data ---
    cells = pd.read_csv(cx_cells_csv)
    with open(cxcx_synapses_json, 'r') as f:
        syn_data = json.load(f)

    # Soma coordinates
    coords = {int(row.cell_id): np.array([row.x, row.y, row.z]) for _, row in cells.iterrows()}
    cell_ids = list(coords.keys())

    # Connected pairs (set of tuples)
    connected_pairs = set()
    for post_id, syn_list in syn_data.items():
        post_id = int(post_id)
        for s in syn_list:
            pre_id = int(s["pre_id"])
            connected_pairs.add((pre_id, post_id))

    # --- Compute all pair distances ---
    results = []
    for pre_id, post_id in product(cell_ids, cell_ids):
        if pre_id == post_id:
            continue
        p1, p2 = coords[pre_id], coords[post_id]
        dist = np.linalg.norm(p1 - p2)
        connected = (pre_id, post_id) in connected_pairs
        results.append((dist, connected))

    df = pd.DataFrame(results, columns=["dist", "connected"])

    # --- Bin distances ---
    if max_dist is None:
        max_dist = df["dist"].max()
    bins = np.arange(0, max_dist + bin_size, bin_size)
    df["bin"] = pd.cut(df["dist"], bins=bins, right=False)

    grouped = df.groupby("bin").agg(
        n_total=("connected", "count"),
        n_connected=("connected", "sum")
    ).reset_index()

    grouped["conn_prob"] = grouped["n_connected"] / grouped["n_total"]
    grouped["bin_start"] = grouped["bin"].apply(lambda b: b.left)
    grouped["bin_end"] = grouped["bin"].apply(lambda b: b.right)
    grouped = grouped[["bin_start", "bin_end", "n_total", "n_connected", "conn_prob"]]

    return grouped

def plot_conn_prob_with_distance(cx_cells_csv, cxcx_synapses_json, savepath,
                                 bin_size=10, max_dist=1000, show=False):
    """
    Compute and plot connection probability vs. intersomatic distance as a histogram.

    Args:
        cx_cells_csv (str): Path to CSV file with columns ['cell_id', 'x', 'y', 'z', ...].
        cxcx_synapses_json (str): Path to JSON file mapping post_id -> list of synapses.
        savepath: Output path for the saved plot (e.g. './plots/conn_prob_dist.png').
        bin_size (int): Bin size for distance histogram (in µm).
        max_dist (float or None): Maximum distance cutoff.
        show (bool): Whether to display the plot interactively.
    """
    # Compute connection probabilities per bin
    df = _conn_prob_with_distance(cx_cells_csv, cxcx_synapses_json,
                                 bin_size=bin_size, max_dist=max_dist)

    bins = df["bin_start"].values
    probs = df["conn_prob"].values

    # Use bar plot for histogram
    plt.figure(figsize=(8, 4))
    plt.bar(bins, probs, width=bin_size, edgecolor='black', align='edge')
    plt.xlabel("Distance between somata (µm)")
    plt.ylabel("Connection probability")
    plt.title("Connection Probability vs Distance")
    plt.xlim(0, 1000)
    plt.ylim(0, 1)
    plt.yticks(np.arange(0, 1.1, 0.1))
    plt.grid(True, linestyle='--', alpha=0.6, axis='y')
    plt.tight_layout()

    # Save figure
    plt.savefig(savepath, dpi=300)
    if show:
        plt.show()
    plt.close()
