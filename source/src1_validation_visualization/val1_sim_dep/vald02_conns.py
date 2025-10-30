from pathlib import Path

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import json
from collections import Counter
from itertools import product

from source.src2_utils.ut0_random_manager import np

def conn_matrix_mean_syn(csv_file, save_path, layer_comp_params):
    """
    Plot cortico-cortical connectivity matrix showing the MEAN number of synapses per connection type,
    aggregated by layer + m-type (ignoring e-type).

    Args:
        csv_file (str): Path to CSV file with columns ['pre_id', 'post_id', 'pre_me_type', 'post_me_type'].
        save_path: Path to save figure; if None, display interactively.
        layer_comp_params (dict): e.g. {"L23": {"exc": ..., "inh": ...}, "L5": {...}}
    """
    # Load CSV
    df = pd.read_csv(csv_file)

    # Helper: collapse full type (layer_mtype_etype) -> layer_mtype
    def collapse_type(full_type):
        parts = full_type.split("_")
        if len(parts) >= 3:
            return "_".join(parts[:-1])  # remove last (etype)
        else:
            return full_type  # fallback

    # Add collapsed columns
    df["pre_collapsed"] = df["pre_me_type"].apply(collapse_type)
    df["post_collapsed"] = df["post_me_type"].apply(collapse_type)

    # Generate all unique collapsed types
    cell_types = sorted(df["pre_collapsed"].unique().tolist() + df["post_collapsed"].unique().tolist())
    cell_types = sorted(set(cell_types))

    # Initialize matrix
    conn_matrix_mean = pd.DataFrame(0.0, index=cell_types, columns=cell_types)

    # Compute mean synapses per connection (aggregated)
    for pre_type in cell_types:
        for post_type in cell_types:
            sub = df[(df["pre_collapsed"] == pre_type) & (df["post_collapsed"] == post_type)]

            if len(sub) == 0:
                mean_syn = 0.0
            else:
                # Count unique pairs (pre→post)
                n_unique_conns = sub[["pre_id", "post_id"]].drop_duplicates().shape[0]
                mean_syn = len(sub) / n_unique_conns if n_unique_conns > 0 else 0.0

            conn_matrix_mean.loc[pre_type, post_type] = mean_syn

    # Plot heatmap
    plt.figure(figsize=(12, 10))
    ax = sns.heatmap(
        conn_matrix_mean,
        annot=True,
        fmt=".2f",
        cmap="viridis",
        cbar_kws={'label': 'Mean synapses per connection'},
        linewidths=0.5,
        linecolor='black',
        square=True
    )

    ax.set_title("Empirical cortico-cortical connectivity (mean synapses per connection, collapsed e-types)",
                 fontsize=13, pad=20)
    ax.set_xlabel("Post-synaptic cell type", fontsize=12, labelpad=10)
    ax.set_ylabel("Pre-synaptic cell type", fontsize=12, labelpad=10)

    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)
    plt.tight_layout()

    # Save or show
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.close()
    else:
        plt.show()

def conn_matrix_theo(json_file, layer_comp_params, save_path=None):
    """
    Plot theoretical cortico-cortical connectivity matrix based on
    the mean number of synapses per connection type.

    Args:
        json_file (str): Path to JSON file containing entries like:
                         { "pre:post": { "mean_number_of_synapse_per_connection": float }, ... }
        layer_comp_params (dict): e.g. {"L23": {"exc": ..., "inh": ...}, "L5": {...}}
        save_path (str | None): Path to save figure; if None, display interactively.
    """
    # Load JSON
    with open(json_file, "r") as f:
        data = json.load(f)

    # Generate all full cell types (sorted for consistent axis order)
    cell_types = [f"{layer}_{'_'.join(ctype.split('_')[:-1])}"
                  for layer, types in layer_comp_params.items()
                  for ctype in types.keys()]
    cell_types = sorted(set(cell_types))

    # Initialize matrix
    conn_matrix_mean = pd.DataFrame(0.0, index=cell_types, columns=cell_types)

    # Fill in values from JSON
    for pair, values in data.items():
        pre, post = pair.split(":")
        if pre in cell_types and post in cell_types:
            conn_matrix_mean.loc[pre, post] = values.get("mean_number_of_synapse_per_connection", 0.0)

    # Plot heatmap
    plt.figure(figsize=(12, 10))
    ax = sns.heatmap(
        conn_matrix_mean,
        annot=True,
        fmt=".2f",
        cmap="viridis",
        cbar_kws={'label': 'Mean synapses per connection'},
        linewidths=0.5,
        linecolor='black',
        square=True
    )

    # Titles and labels
    title = "Target cortico-cortical connectivity matrix: mean synapses per connection"
    ax.set_title(title, fontsize=13, pad=20)
    ax.set_xlabel("Post-synaptic cell type", fontsize=12, labelpad=10)
    ax.set_ylabel("Pre-synaptic cell type", fontsize=12, labelpad=10)

    # Rotate tick labels for readability
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)

    # Layout
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
        cbar_kws={'label': 'Connection ratio (%)' if unique_connections else "Synapses fraction (%)"},
        linewidths=0.5,
        linecolor='black',
        square=True
    )

    # Titles and labels
    if unique_connections==False:
        title = "Thalamocortical connectivity matrix: synapses and row-wise ratio (%)"
    else:
        title = "Thalamocortical connectivity matrix: connections and row-wise ratio (%)"
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

def interbouton_int_histogram(cx_pop_csv: Path, bins: int = 20, save_path=None):
    """
    Creates a normalized histogram (probability distribution) of b_int for the cells in cx_pop_csv.
    Saves *.png.

    Args:
        cx_pop_csv (Path): Path to CSV with cortical cells containing 'b_int'.
        bins (int): Number of bins in the histogram (default=20).
    """

    df = pd.read_csv(cx_pop_csv)
    if "b_int" not in df.columns:
        raise KeyError("'b_int' column not found in the CSV.")

    data = df["b_int"].replace([np.inf, -np.inf], np.nan).dropna()

    # Map filenames to descriptive titles and colors
    title_map = {
        "cx04_pop_pre.csv": "pre-pruning",
        "cx04_pop_post1.csv": "after General Pruning [01]",
        "cx04_pop_post2.csv": "after Multisynapse Pruning [02]",
        "cx04_pop_post3.csv": "after Plasticity-Reserve Pruning [03]"
    }

    color_map = {
        "cx04_pop_pre.csv": "red",
        "cx04_pop_post1.csv": "orange",
        "cx04_pop_post2.csv": "gold",    # softer yellow for visibility
        "cx04_pop_post3.csv": "green"
    }

    # Get descriptive title and color
    file_name = cx_pop_csv.name
    plot_title = title_map.get(file_name, file_name)
    color = color_map.get(file_name, "gray")

    plt.figure(figsize=(8, 5))
    bin_edges = np.linspace(0, 100, bins + 1)  # 20 bins between 0 and 100 µm
    plt.hist(
        data,
        bins=bin_edges,
        density=True,              # normalize histogram
        edgecolor="black",
        color=color,
        alpha=0.8
    )

    # Axis settings
    plt.xlim(0, 100)
    plt.xticks(np.arange(0, 101, 10))
    plt.grid(True, linestyle='--', linewidth=0.5, alpha=0.7)
    plt.xlabel("Interbouton interval (µm)")
    plt.ylabel("P")
    plt.title(f"Distribution of interbouton interval ({plot_title})")
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.close()


def syn_per_conn_histogram(synapses_json: Path, savepath):
    """
    Creates a normalized histogram (PMF) of synapses per connection for all connections.
    Saves it as PNG.
    """

    # Map filenames to descriptive titles and colors
    title_map = {
        "cxcx01_appositions.json": "pre-pruning",
        "prune01_cxcx.json": "after General Pruning [01]",
        "prune02_cxcx.json": "after Multisynapse Pruning [02]",
        "prune03_cxcx.json": "after Plasticity-Reserve Pruning [03]"
    }

    color_map = {
        "cxcx01_appositions.json": "red",
        "prune01_cxcx.json": "orange",
        "prune02_cxcx.json": "gold",
        "prune03_cxcx.json": "green"
    }

    # Determine title and color
    file_name = synapses_json.name
    main_title = f"Corticocortical synapses per connection distribution ({title_map.get(file_name, file_name)})"
    color = color_map.get(file_name, "gray")

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
    bins = np.arange(0, 21, 1)  # 1-synapse bins (0–20)
    plt.hist(
        synapse_nums,
        bins=bins,
        edgecolor='black',
        density=True,
        align='left',
        color=color,
        alpha=0.8
    )

    # Axis settings
    plt.xlabel("Synapses per connection")
    plt.ylabel("P")
    plt.title(main_title)
    plt.xlim(0, 20)
    plt.xticks(np.arange(0, 21, 2))
    plt.ylim(0, 0.15)
    plt.yticks(np.arange(0, 0.151, 0.025))
    plt.grid(True, linestyle='--', linewidth=0.5, alpha=0.7)

    plt.tight_layout()
    plt.savefig(savepath, dpi=300, bbox_inches='tight')
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

    # --- Title and color maps ---
    title_map = {
        "conn_prob_pre.jpg": "pre-pruning",
        "conn_prob_post1.jpg": "after General Pruning [01]",
        "conn_prob_post2.jpg": "after Multisynapse Pruning [02]",
        "conn_prob_post3.jpg": "after Plasticity-Reserve Pruning [03]"
    }

    color_map = {
        "conn_prob_pre.jpg": "red",
        "conn_prob_post1.jpg": "orange",
        "conn_prob_post2.jpg": "gold",
        "conn_prob_post3.jpg": "green"
    }

    savepath = Path(savepath)
    plot_title = title_map.get(savepath.name, savepath.name)
    color = color_map.get(savepath.name, "gray")

    # --- Plot ---
    plt.figure(figsize=(8, 4))
    plt.bar(
        bins, probs, width=bin_size,
        edgecolor='black', align='edge',
        color=color, alpha=0.8
    )
    plt.xlabel("Distance between somata (µm)")
    plt.ylabel("Connection probability")
    plt.title(f"Cortical cell connection probability vs. distance ({plot_title})")

    # Axis styling
    plt.xlim(0, max_dist)
    plt.xticks(np.arange(0, max_dist + 1, 100))
    plt.ylim(0, 1)
    plt.yticks(np.arange(0, 1.1, 0.1))
    plt.grid(True, linestyle='--', linewidth=0.5, alpha=0.7, axis='y')

    plt.tight_layout()
    plt.savefig(savepath, dpi=300, bbox_inches='tight')

    if show:
        plt.show()
    plt.close()