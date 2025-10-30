import os
from neuron import h
import matplotlib.pyplot as plt
import glob
from matplotlib.ticker import ScalarFormatter

from source.src0_core.cr0_model_setup.m01_cx_cells.CellTopology import CellTopology

from source.src2_utils.ut0_random_manager import np

def _run_single_iclamp_cx(cell_type, cell_templates, tstop=1800, delay=600, dur=600, amp=1):
    """
    Run single-cell simulation with IClamp at soma.
    Returns time vector, voltage vector, and current vector.
    """
    cell = cell_templates[cell_type](0)
    soma = cell.soma[0]

    stim = h.IClamp(soma(0.5))
    stim.delay = delay
    stim.dur = dur
    stim.amp = amp

    # record vectors
    t_vec = h.Vector().record(h._ref_t)
    v_vec = h.Vector().record(soma(0.5)._ref_v)
    i_vec = h.Vector().record(stim._ref_i)  # record injected current

    h.finitialize(-70)
    h.continuerun(tstop)

    return list(t_vec), list(v_vec), list(i_vec)

def run_iclamp_cx(save_dir, cell_templates, tstop=1500, delay=500, dur=500, amp=1.5):
    """
    Runs IClamp simulation for all cells in cell_templates,
    saves plots with Vm and injected current.
    """
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)

    for cell_type in cell_templates.keys():
        print(f"Running {cell_type}...")
        t, v, i = _run_single_iclamp_cx(cell_type, cell_templates, tstop, delay, dur, amp)

        # Create figure with two vertically stacked plots
        fig, (ax1, ax2) = plt.subplots(
            2, 1, figsize=(14, 6), sharex=True,
            gridspec_kw={'height_ratios': [5, 1]}
        )

        # Determine tick positions every 100 ms
        xticks = np.arange(0, tstop + 100, 100)

        # --- Somatic membrane potential ---
        ax1.plot(t, v, linewidth=1, color='black')
        ax1.set_ylabel("Somata membrane potential [mV]", fontsize=11)
        ax1.set_ylim(-100, 40)
        ax1.grid(True, linestyle='--', linewidth=0.5, alpha=0.7)
        ax1.set_title(f"Cortical {cell_type} response to injected current at the somata", fontsize=13, pad=15)

        # Vertical dashed lines every 100 ms
        for x in xticks:
            ax1.axvline(x=x, color='gray', linestyle='--', linewidth=0.5, alpha=0.4)

        # Remove x-tick labels and ticks
        ax1.set_xticks([])
        ax1.tick_params(axis='x', length=0)
        ax1.yaxis.set_label_coords(-0.08, 0.5)

        # --- Injected current ---
        ax2.plot(t, i, color='red', linewidth=1)
        ax2.set_xlabel("Time (ms)", fontsize=11)
        ax2.set_ylabel("Current injected [nA]", fontsize=11)
        ax2.set_ylim(-0.5, 2.5)
        ax2.grid(True, linestyle='--', linewidth=0.5, alpha=0.7)

        # Set x-ticks every 100 ms on the bottom plot
        ax2.set_xticks(xticks)

        # Align y-label position
        ax2.yaxis.set_label_coords(-0.08, 0.5)

        # Layout and save
        plt.tight_layout()
        filename = os.path.join(save_dir, f"{cell_type}.png")
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Saved {filename}")




def _run_single_iclamp_tc(tc_cell, tstop, delay, dur, amp, v_init):
    """Run IClamp stimulation on a VPMCell (variant = follower/initiator)."""

    # Inject current at soma
    stim = h.IClamp(tc_cell.soma(0.5))
    stim.delay = delay
    stim.dur = dur
    stim.amp = amp

    # Record vectors
    t_vec = h.Vector().record(h._ref_t)
    v_vec = h.Vector().record(tc_cell.soma(0.5)._ref_v)
    i_vec = h.Vector().record(stim._ref_i)

    # Spike detector
    tc_cell.setup_spike_detector(threshold=0.0)

    # Run
    h.finitialize(v_init)
    h.continuerun(tstop)

    return list(t_vec), list(v_vec), list(i_vec), tc_cell.get_spike_times()

def run_iclamp_tc(save_dir, tc_cell, tstop=1500, delay=500, dur=500, amp=1.5, v_init=-70):
    """
    Runs IClamp simulation for a thalamic relay (TC) cell and saves
    the membrane potential and injected current plots.

    Visual formatting consistent with cortical IClamp plots.
    """
    os.makedirs(save_dir, exist_ok=True)

    t, v, i, spikes = _run_single_iclamp_tc(tc_cell=tc_cell, tstop=tstop, delay=delay, dur=dur, amp=amp, v_init=v_init)

    # Create figure
    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(12, 6), sharex=True,
        gridspec_kw={'height_ratios': [5, 1]}
    )

    # --- Somatic membrane potential ---
    ax1.plot(t, v, linewidth=1, color='black')
    ax1.set_ylabel("Somata membrane potential [mV]", fontsize=11)
    ax1.set_ylim(-130, 90)
    ax1.grid(True, linestyle='--', linewidth=0.5, alpha=0.7)
    ax1.set_xticklabels([])  # hide x-ticks for top plot
    ax1.tick_params(axis='x', length=0)
    ax1.set_title("Thalamic relay cell (VPM) response to injected current at the somata",
                  fontsize=13, pad=15)

    # Align y-label
    ax1.yaxis.set_label_coords(-0.08, 0.5)

    # --- Injected current ---

    ax2.plot(t, i, color='red', linewidth=1)
    ax2.set_xlabel("Time [ms]", fontsize=11)
    ax2.set_ylabel("Current injected [nA]", fontsize=11)
    ax2.set_ylim(-0.5, 3.5)
    ax2.grid(True, linestyle='--', linewidth=0.5, alpha=0.7)

    # Align y-label
    ax2.yaxis.set_label_coords(-0.08, 0.5)

    # Final layout adjustments
    plt.tight_layout()

    # Save figure
    filename = os.path.join(save_dir, "VPMCell.png")
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved {filename}")


def plot_dend_axon_z_dist(
    cell_templates_path: str,
    outdir: str,
    rotation=(np.pi/2, 0, 0),
    bins: int = 40,
    density: bool = True,
    mode: str = "stacked",  # or "overlay"
    verbose: bool = True
):
    """
    For each cell template in `cell_templates_path`, loads its morphology (.asc),
    extracts dendrite and axon 3D coordinates, computes z-axis distributions,
    and saves histogram plots to `outdir`.

    When `density=True`, each histogram (dendrite, axon) is normalized
    so its integral over z equals 1 — effectively representing a PMF density.
    """

    # Collect .asc files
    tasks = []
    for entry in sorted(os.listdir(cell_templates_path)):
        template_dir = os.path.join(cell_templates_path, entry)
        if not os.path.isdir(template_dir):
            continue
        morph_dir = os.path.join(template_dir, "morphology")
        if not os.path.isdir(morph_dir):
            if verbose:
                print(f"[skip] no morphology/ in {entry}")
            continue
        asc_files = sorted(glob.glob(os.path.join(morph_dir, "*.asc")))
        if not asc_files:
            if verbose:
                print(f"[skip] no .asc in {morph_dir}")
            continue
        tasks.append((entry, asc_files[0]))

    bin_range = (-1500, 1500)
    os.makedirs(outdir, exist_ok=True)
    saved = []

    for template_name, asc_path in tasks:
        try:
            # Initialize morphology
            cell = CellTopology(morph_path=asc_path, rotation=rotation)

            # Get 3D points
            dend_points = cell.get_dendrite_points()
            axon_points = cell.get_axon_points()

            # Skip empty morphologies
            if len(dend_points) == 0 and len(axon_points) == 0:
                if verbose:
                    print(f"[skip] {template_name}: no axon/dendrite points")
                continue

            # Compute z-values
            data, labels = [], []
            if len(dend_points):
                data.append(dend_points[:, 2])
                labels.append("dendrite")
            if len(axon_points):
                data.append(axon_points[:, 2])
                labels.append("axon")

            # --- Plot ---
            fig, ax = plt.subplots(figsize=(8, 4.5))
            colors = ["tab:green", "tab:blue"]

            if mode == "stacked" and len(data) > 1:
                ax.hist(
                    data,
                    bins=bins,
                    range=bin_range,
                    density=density,
                    stacked=True,
                    label=labels,
                    alpha=0.8,
                    color=colors[:len(data)],
                    edgecolor="black",
                    linewidth=0.3,
                )
            else:
                for z_values, label, color in zip(data, labels, colors):
                    ax.hist(
                        z_values,
                        bins=bins,
                        range=bin_range,
                        density=density,
                        alpha=0.5,
                        label=label,
                        histtype="stepfilled" if mode == "overlay" else "bar",
                        color=color,
                        edgecolor="black",
                        linewidth=0.4,
                    )

            # --- Axes and labels ---
            ax.set_xlabel("z [µm]", fontsize=11)
            ax.set_ylabel("PMF density", fontsize=11)
            ax.set_title(f"{template_name} cortical cell morphology distribution with respect to z-axis (depth)",
                         fontsize=12, pad=12)

            # Scientific notation for small densities
            ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
            ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))

            # Grid styling
            ax.grid(True, which='major', linestyle='--', linewidth=0.5, alpha=0.6)
            ax.legend()

            # Save
            out_path = os.path.join(outdir, f"{template_name}.jpg")
            fig.savefig(out_path, dpi=300, bbox_inches="tight")
            plt.close(fig)

            saved.append((template_name, out_path))
            if verbose:
                print(f"[ok] saved {out_path}")

        except Exception as e:
            if verbose:
                print(f"[err] {template_name}: {e}")

    return saved