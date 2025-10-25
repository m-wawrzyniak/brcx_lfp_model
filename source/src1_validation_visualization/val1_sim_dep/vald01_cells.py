import json
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.ticker import MaxNLocator
from collections import defaultdict
import math
from typing import Dict
from pprint import pprint

from source.src2_utils.ut0_random_manager import np

from config_templates import conf0_model_parameters as conf0

STIM_COLORS = {
    'r': 'white',
    'wk': '#fff7a0',   # yellowish
    'str': '#ff9999'   # reddish
}

def cx_plot_me_composition_emp(json_path, save_path=None):
    # load data
    with open(json_path, "r") as f:
        data = json.load(f)

    layers = ["L23", "L4", "L5", "Whole population"]

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

        # custom autopct to make percentages bold
        def bold_autopct(pct):
            return f"{pct:.1f}%"

        # adjust pctdistance for whole population to reduce overlap
        pctdist = 0.7 if layer == "Whole population" else 0.6

        wedges, texts, autotexts = ax.pie(
            ordered_sizes,
            labels=ordered_labels,
            autopct=bold_autopct,
            startangle=90,
            colors=colors,
            wedgeprops={'edgecolor':'black', 'linewidth':0.5},  # black contour
            textprops={"fontsize": 8},
            labeldistance=1.05,
            pctdistance=pctdist
        )

        # rotate labels for readability
        for t in texts:
            t.set_rotation(30)  # adjust angle as needed

        # bold percentages
        for at in autotexts:
            at.set_fontweight('bold')

        ax.set_title(f"{layer} (n_cells={data[layer]['tot_n_cells']})", fontsize=12)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300)
        plt.close()
    else:
        plt.show()


def cx_plot_me_composition_theo(layer_comp_params, save_path=None):
    """
    Creates one figure with pie charts for L23, L4, L5, and population
    based on target ratios. Excitatory (PC, TTPC, SS) are reddish shades,
    inhibitory are bluish shades. Each slice is progressively deeper in color.
    """
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib import cm

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

        # bold percentages
        def bold_autopct(pct):
            return f"{pct:.1f}%"

        # adjust pctdistance for whole population to reduce overlap
        pctdist = 0.7 if layer == "population" else 0.6

        wedges, texts, autotexts = ax.pie(
            ordered_sizes,
            labels=ordered_labels,
            autopct=bold_autopct,
            startangle=90,
            colors=colors,
            wedgeprops={'edgecolor':'black', 'linewidth':0.5},  # black contour
            textprops={"fontsize": 8},
            labeldistance=1.05,
            pctdistance=pctdist
        )

        # rotate labels for readability
        for t in texts:
            t.set_rotation(30)

        # bold percentages
        for at in autotexts:
            at.set_fontweight('bold')

        ax.set_title(f"{layer} target composition", fontsize=12)

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300)
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
    fig.suptitle("Somata location distribution for all me-types", fontsize=14, y=1.02)

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

        # Make x-axis ticks integers only
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))

    # Layout tuning
    plt.subplots_adjust(wspace=0.05, bottom=0.1, top=0.9)
    plt.tight_layout()

    # Save figure
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved somata z-distribution figure to: {save_path}")


def plot_rasterplot(
        save_path: Path,
        paradigm: dict,
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
    spike_height = 1.0  # slightly taller spikes
    for desc in sorted(cx_grouped.keys()):
        color = color_map[desc]
        for cell_id, spikes in sorted(cx_grouped[desc], key=lambda x: x[0]):
            ax.vlines(spikes, y_pos + 0.5, y_pos + 0.5 + spike_height, color=color)
            yticks.append(y_pos + 0.5 + spike_height / 2)
            ylabels.append(str(cell_id))
            y_pos += 1

    # --- Plot VPM spikes (red) ---
    for i, spikes in enumerate(vpm_spikes):
        ax.vlines(spikes, y_pos + 0.5, y_pos + 0.5 + spike_height, color='red')
        yticks.append(y_pos + 0.5 + spike_height / 2)
        ylabels.append(vpm_ids[i])
        y_pos += 1

    # --- Plot PRV spikes (green) ---
    for i, spikes in enumerate(prv_spikes):
        ax.vlines(spikes, y_pos + 0.5, y_pos + 0.5 + spike_height, color='green')
        yticks.append(y_pos + 0.5 + spike_height / 2)
        ylabels.append(prv_ids[i])
        y_pos += 1

    # --- Final touches ---
    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels, fontsize=6)
    ax.set_xlabel("Time (ms)")
    ax.set_ylabel("Cell ID")
    ax.set_title("Rasterplot based on model activity during whisker stimulation")
    ax.set_ylim(0.5, total_cells + 0.5)

    # Vertical grid at major ticks
    ax.grid(axis='x', which='major', linestyle='--', color='gray', alpha=0.5)

    # Legend for cortical cell groups
    handles = [plt.Line2D([0], [0], color=color_map[desc], lw=2) for desc in sorted(cx_grouped.keys())]
    ax.legend(handles, sorted(cx_grouped.keys()), title="Cortical cell types", bbox_to_anchor=(1.05, 1), loc='upper left')

    plt.tight_layout()
    plt.savefig(save_path, dpi=600, format="pdf")
    plt.close()


def calc_goal_densities(layer_comp_target:dict[str, list],
                        tissue_params:dict,
                        goal_n_cells:int,
                        bbp_layer_info_path:str,
                        scale_factor:float):
    """
    1. Given which layers should be included in the model, and how many cells should it have,
        calculates how many cells should be in specific cell, e.g.:
        alpha_4 = n_cells_l4_emp / sum(n_cells_emp for all layers which should be included)
        n_cells_l4_mod = goal_n_cells*alpha_4

    2. Given goal e-types existing within a layer, calculates goal density in this layer:
        mu_4 = sum(n_cells_l4_emp for e-type cells which will be in model) / n_cells_l4_emp_total
        Then goal density is:
        ro_4_mod = mu_4 * ro_4_emp

    3. Given n_cells_mod and ro_mod for specific layer, gets what the radius of the column should be e.g.:
        V_cyl_4 = n_cells_l4_mod / ro_4_mod
        R_goal = sqrt(V_cyl_4 / pi*H)

    Verification:
        alpha_4 = 5792 /
    """
    # ---- load BBP e-type counts ----
    with open(bbp_layer_info_path, "r") as f:
        bbp = json.load(f)

    # Build: layer -> {etype: count} and totals per layer
    bbp_etypes_by_layer: Dict[str, Dict[str, float]] = {}
    bbp_total_by_layer: Dict[str, float] = {}
    for layer_name, payload in bbp.items():
        et_map = payload.get("No. of neurons per electrical types", {})
        if isinstance(et_map, dict) and et_map:
            et_map_f = {str(k): float(v) for k, v in et_map.items()}
            bbp_etypes_by_layer[layer_name] = et_map_f
            bbp_total_by_layer[layer_name] = float(sum(et_map_f.values()))

    included_layers = list(layer_comp_target.keys())
    if not included_layers:
        raise ValueError("layer_comp_target is empty—nothing to compute.")

    # ---- Step 1: alpha and goal allocation from BBP totals ----
    # Use only included layers for the denominator
    emp_counts = {}
    for layer in included_layers:
        if layer not in bbp_total_by_layer:
            raise KeyError(f"Layer '{layer}' missing in BBP file (e-type totals).")
        emp_counts[layer] = bbp_total_by_layer[layer]

    total_emp_included = sum(emp_counts.values())
    if total_emp_included <= 0:
        raise ValueError("Total BBP empirical counts for included layers is non-positive.")

    provisional = []
    for layer in included_layers:
        n_emp = emp_counts[layer]
        alpha = n_emp / total_emp_included
        n_goal_float = goal_n_cells * alpha
        provisional.append((layer, alpha, n_goal_float))

    # Largest remainder to ensure integer sum equals goal_n_cells
    floor_alloc = {layer: int(math.floor(nf)) for layer, _, nf in provisional}
    assigned = sum(floor_alloc.values())
    remaining = goal_n_cells - assigned
    remainders = sorted(
        [(layer, nf - math.floor(nf)) for layer, _, nf in provisional],
        key=lambda x: x[1],
        reverse=True,
    )
    for i in range(remaining):
        floor_alloc[remainders[i % len(remainders)][0]] += 1

    # ---- Step 2 & 3: mu, ro_goal, R_goal (using tissue_params D and H) ----
    def um_to_mm(x_um: float) -> float:
        return x_um / 1000.0

    results = {}
    for layer, alpha, _ in provisional:
        chosen_etypes = layer_comp_target.get(layer, [])
        et_counts = bbp_etypes_by_layer.get(layer, {})
        total_layer_etype_emp = float(sum(et_counts.values())) if et_counts else 0.0

        if not chosen_etypes or total_layer_etype_emp <= 0:
            mu = 0.0
        else:
            chosen_sum = sum(et_counts.get(et, 0.0) for et in chosen_etypes)
            mu = chosen_sum / total_layer_etype_emp

        # Pull geometry/density for this layer
        if layer not in tissue_params:
            raise KeyError(f"Layer '{layer}' missing in tissue_params.")
        params = tissue_params[layer]
        H_um = float(params["H"])
        D_th = float(params["D"])  # thousands / mm^3
        ro_emp = D_th * 1000.0     # cells / mm^3
        ro_goal = mu * ro_emp * scale_factor

        n_goal = floor_alloc[layer]
        H_mm = um_to_mm(H_um)

        if ro_goal > 0 and H_mm > 0:
            V_needed = n_goal / ro_goal  # mm^3
            R_goal_mm = math.sqrt(V_needed / (math.pi * H_mm))
            R_goal_um = 1000.0 * R_goal_mm
        else:
            R_goal_um = 0.0

        # Save useful context (also echo the BBP empirical numbers we used)
        results[layer] = {
            "bbp_n_emp": emp_counts[layer],     # BBP empirical layer total (e-type sum)
            "alpha": alpha,
            "n_goal": int(n_goal),

            "mu": float(mu),
            "ro_emp": ro_emp,                   # cells/mm^3
            "ro_goal": ro_goal,                 # cells/mm^3
            "H_um": H_um,
            "R_goal_um": R_goal_um,
        }

    results["_totals"] = {
        "goal_n_cells": int(goal_n_cells),
        "allocated_sum": int(sum(results[l]["n_goal"] for l in included_layers)),
        "included_layers": included_layers,
        "bbp_total_included": float(total_emp_included),
    }
    return results

if __name__ == "__main__":
    res = calc_goal_densities(
        layer_comp_target=conf0.LAYER_COMP_TARGET,
        tissue_params=conf0.TISSUE_PARAMS,
        goal_n_cells=conf0.GOAL_N_CELLS,
        bbp_layer_info_path="/home/mateusz-wawrzyniak/PycharmProjects/brcx_lfp_model/source/src0_core/cr0_model_setup/m00_bbp_parameters/layer_info.json",
        scale_factor=conf0.SCALE_FACTOR
    )
    pprint(res)