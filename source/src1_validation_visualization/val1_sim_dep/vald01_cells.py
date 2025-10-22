import json
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from collections import defaultdict

from source.src2_utils.ut0_random_manager import np

STIM_COLORS = {
    'r': 'white',
    'wk': '#fff7a0',   # yellowish
    'str': '#ff9999'   # reddish
}

def cx_plot_me_composition_emp(json_path, save_path=None):
    """
    Loads JSON with layer/population info and creates one figure with pie charts
    for L23, L4, L5, and population. Excitatory (PC, TTPC, SS) are reddish shades,
    inhibitory are bluish shades. Each slice is progressively deeper in color.
    """
    # load data
    with open(json_path, "r") as f:
        data = json.load(f)

    layers = ["L23", "L4", "L5", "population"]

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    axes = axes.flatten()

    for i, layer in enumerate(layers):
        ax = axes[i]
        layer_data = data[layer]["me_types"]

        labels = list(layer_data.keys())
        sizes = list(layer_data.values())

        # split into excitatory and inhibitory
        excitatory = [(lbl, val) for lbl, val in zip(labels, sizes)
                      if any(exc in lbl for exc in ["PC", "TTPC", "SS"])]
        inhibitory = [(lbl, val) for lbl, val in zip(labels, sizes)
                      if not any(exc in lbl for exc in ["PC", "TTPC", "SS"])]

        # create progressive red shades for excitatory
        colors = []
        ordered_labels = []
        ordered_sizes = []

        if excitatory:
            reds = cm.Reds(np.linspace(0.5, 0.9, len(excitatory)))  # lighter → deeper red
            for (lbl, val), col in zip(excitatory, reds):
                ordered_labels.append(lbl)
                ordered_sizes.append(val)
                colors.append(col)

        if inhibitory:
            blues = cm.Blues(np.linspace(0.5, 0.9, len(inhibitory)))  # lighter → deeper blue
            for (lbl, val), col in zip(inhibitory, blues):
                ordered_labels.append(lbl)
                ordered_sizes.append(val)
                colors.append(col)

        wedges, texts, autotexts = ax.pie(
            ordered_sizes,
            labels=ordered_labels,
            autopct="%.1f%%",
            startangle=90,
            colors=colors,
            textprops={"fontsize": 8}
        )
        ax.set_title(f"{layer} (n={data[layer]['tot_n_cells']})", fontsize=12)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=200)
        plt.close()
    else:
        plt.show()


def cx_plot_me_composition_theo(layer_comp_params, save_path=None):
    """
    Creates one figure with pie charts for L23, L4, L5, and population
    based on target ratios. Excitatory (PC, TTPC, SS) are reddish shades,
    inhibitory are bluish shades. Each slice is progressively deeper in color.
    """
    layers = ["L23", "L4", "L5"]
    layer_comp_params = dict(layer_comp_params)

    # add population by averaging across layers
    population = {}
    for layer in layers:
        for k, v in layer_comp_params[layer].items():
            population[k] = population.get(k, 0) + v
    # normalize to sum=1
    total = sum(population.values())
    for k in population:
        population[k] /= total
    layer_comp_params["population"] = population

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    axes = axes.flatten()

    for i, layer in enumerate(["L23", "L4", "L5", "population"]):
        ax = axes[i]
        layer_data = layer_comp_params[layer]

        labels = list(layer_data.keys())
        sizes = list(layer_data.values())

        # split into excitatory and inhibitory
        excitatory = [(lbl, val) for lbl, val in zip(labels, sizes)
                      if any(exc in lbl for exc in ["PC", "TTPC", "SS"])]
        inhibitory = [(lbl, val) for lbl, val in zip(labels, sizes)
                      if not any(exc in lbl for exc in ["PC", "TTPC", "SS"])]

        # create progressive red/blue shades
        colors, ordered_labels, ordered_sizes = [], [], []

        if excitatory:
            reds = cm.Reds(np.linspace(0.5, 0.9, len(excitatory)))
            for (lbl, val), col in zip(excitatory, reds):
                ordered_labels.append(lbl)
                ordered_sizes.append(val)
                colors.append(col)

        if inhibitory:
            blues = cm.Blues(np.linspace(0.5, 0.9, len(inhibitory)))
            for (lbl, val), col in zip(inhibitory, blues):
                ordered_labels.append(lbl)
                ordered_sizes.append(val)
                colors.append(col)

        wedges, texts, autotexts = ax.pie(
            ordered_sizes,
            labels=ordered_labels,
            autopct="%.1f%%",
            startangle=90,
            colors=colors,
            textprops={"fontsize": 8}
        )
        ax.set_title(f"{layer} target composition", fontsize=12)

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=200)
        plt.close()
    else:
        plt.show()


def cx_somata_dist_z(cx_cells_path: str, save_path, n_bins:int = 40):
    """
    Plots vertical distribution histograms (z-axis) of cell somata for each unique 'desc' type
    and saves the resulting figure.

    - One shared y-axis (z) with ticks and label.
    - Subplots close together.
    - Each histogram in a different color from one colormap.
    """
    # Load data
    df = pd.read_csv(cx_cells_path)
    if not {'z', 'desc'}.issubset(df.columns):
        raise ValueError("CSV must contain at least 'z' and 'desc' columns.")

    # Group by description
    desc_groups = df.groupby('desc')
    desc_names = sorted(desc_groups.groups.keys())
    n_desc = len(desc_names)

    # Color palette
    cmap = cm.get_cmap('viridis', n_desc)
    colors = [cmap(i) for i in range(n_desc)]

    # Global z limits for consistent scaling
    z_min, z_max = -2000, 0
    bins = np.linspace(z_min, z_max, n_bins + 1)

    # Color palette
    cmap = cm.get_cmap('viridis', n_desc)
    colors = [cmap(i) for i in range(n_desc)]

    # Figure setup
    fig, axes = plt.subplots(1, n_desc, figsize=(2.2 * n_desc, 6), sharey=True)
    if n_desc == 1:
        axes = [axes]

    # Plot
    for i, (ax, desc, color) in enumerate(zip(axes, desc_names, colors)):
        z_vals = desc_groups.get_group(desc)['z']
        ax.hist(z_vals, bins=bins, orientation='horizontal', color=color, alpha=0.8)
        ax.set_title(desc, rotation=45, fontsize=9, pad=8)
        ax.set_xlabel('Count', fontsize=9)
        ax.set_ylim(z_min, z_max)
        ax.grid(alpha=0.25)

        # Y-labels and ticks
        if i == 0:
            ax.set_ylabel('z (µm)', fontsize=10)
            ax.tick_params(axis='y', labelsize=8)
        else:
            ax.tick_params(axis='y', left=False, labelleft=False)
        ax.tick_params(axis='x', labelsize=8)

    # Layout tuning
    plt.subplots_adjust(wspace=0.05, bottom=0.1, top=0.9)
    plt.tight_layout()

    # Save figure
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved somata z-distribution figure to: {save_path}")


def plot_rasterplot(
        save_path: Path,
        paradigm:dict,
        cx_spikes_csv: str,
        vpm_spikes_csv: str,
        prv_spikes_json: str,
        cell_pop_path: str
):
    cx_spikes_raw, cx_ids = [], []
    with open(cx_spikes_csv, 'r') as f:
        next(f)
        for line in f:
            cell_id, spikes_str = line.strip().split(',', 1)
            spikes = list(map(float, spikes_str.split())) if spikes_str else []
            cx_ids.append(int(cell_id))
            cx_spikes_raw.append(spikes)

    # --- Load VPM spikes ---
    vpm_spikes, vpm_ids = [], []
    with open(vpm_spikes_csv, 'r') as f:
        next(f)
        for line in f:
            cell_id, spikes_str = line.strip().split(',', 1)
            spikes = list(map(float, spikes_str.split())) if spikes_str else []
            vpm_ids.append(cell_id)
            vpm_spikes.append(spikes)

    # --- Load PRV spikes ---
    prv_spikes, prv_ids = [], []
    with open(prv_spikes_json, 'r') as f:
        prv_data = json.load(f)
        for prv_id, data in prv_data.items():
            spikes = data.get("spike_times", [])
            prv_ids.append(prv_id)
            prv_spikes.append(spikes)

    # --- Load cortical cell population metadata ---
    cell_meta = pd.read_csv(cell_pop_path)
    cx_id_to_desc = dict(zip(cell_meta["cell_id"], cell_meta["desc"]))

    # Group CX spikes by desc
    cx_grouped = defaultdict(list)
    for cell_id, spikes in zip(cx_ids, cx_spikes_raw):
        desc = cx_id_to_desc.get(cell_id, "Unknown")
        cx_grouped[desc].append((cell_id, spikes))

    # Assign colors to desc groups
    color_map = {}
    cmap = plt.get_cmap("tab20")
    for i, desc in enumerate(sorted(cx_grouped.keys())):
        color_map[desc] = cmap(i % cmap.N)

    # --- Prepare plot ---
    total_cells = len(cx_ids) + len(vpm_ids) + len(prv_ids)
    fig_height = max(8, total_cells * 0.05)   # auto-scale height
    fig, ax = plt.subplots(figsize=(12, fig_height))

    # --- Add background for stimulation phases ---
    x_start = 0
    for phase, (phase_type, duration) in paradigm.items():
        if duration == 0:
            continue
        x_end = x_start + duration
        ax.axvspan(x_start, x_end, color=STIM_COLORS.get(phase_type, 'gray'), alpha=0.3)
        x_start = x_end

    # --- Plot CX spikes grouped by desc ---
    y_pos = 0
    yticks, ylabels = [], []
    for desc in sorted(cx_grouped.keys()):
        color = color_map[desc]
        for cell_id, spikes in sorted(cx_grouped[desc], key=lambda x: x[0]):
            ax.vlines(spikes, y_pos + 0.5, y_pos + 1.5, color=color)
            yticks.append(y_pos + 1)
            ylabels.append(str(cell_id))
            y_pos += 1

    # --- Plot VPM spikes (red) ---
    for i, spikes in enumerate(vpm_spikes):
        ax.vlines(spikes, y_pos + 0.5, y_pos + 1.5, color='red')
        yticks.append(y_pos + 1)
        ylabels.append(vpm_ids[i])
        y_pos += 1

    # --- Plot PRV spikes (green) ---
    for i, spikes in enumerate(prv_spikes):
        ax.vlines(spikes, y_pos + 0.5, y_pos + 1.5, color='green')
        yticks.append(y_pos + 1)
        ylabels.append(prv_ids[i])
        y_pos += 1

    # --- Final touches ---
    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels, fontsize=6)
    ax.set_xlabel("Time (ms)")
    ax.set_ylabel("Cell ID")
    ax.set_title("Raster plot with stimulation phases")
    ax.set_ylim(0.5, total_cells + 0.5)

    # Legend for cortical cell groups
    handles = [plt.Line2D([0], [0], color=color_map[desc], lw=2) for desc in sorted(cx_grouped.keys())]
    ax.legend(handles, sorted(cx_grouped.keys()), title="Cortical cell types", bbox_to_anchor=(1.05, 1), loc='upper left')

    plt.tight_layout()
    plt.savefig(save_path, dpi=600, format="pdf")
    plt.close()