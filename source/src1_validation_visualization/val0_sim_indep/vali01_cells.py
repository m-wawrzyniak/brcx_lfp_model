import os
from neuron import h
import matplotlib.pyplot as plt
import glob

from source.src0_core.cr0_model_setup.m01_cx_cells.CellTopology import CellTopology

from source.src2_utils.ut0_random_manager import np

def _run_single_iclamp_cx(cell_type, cell_templates, tstop=1800, delay=600, dur=600, amp=1.5):
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

def run_iclamp_cx(save_dir, cell_templates, tstop=1800, delay=600, dur=600, amp=1.5):
    """
    Runs IClamp simulation for all cells in cell_templates,
    saves plots with Vm and injected current.
    """
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)

    for cell_type in cell_templates.keys():
        print(f"Running {cell_type}...")
        t, v, i = _run_single_iclamp_cx(cell_type, cell_templates, tstop, delay, dur, amp)

        # create figure with two subplots (Vm big, IClamp small)
        fig, (ax1, ax2) = plt.subplots(
            2, 1, figsize=(10, 6), sharex=True,
            gridspec_kw={'height_ratios': [5, 1]}
        )

        # Vm plot
        ax1.plot(t, v, linewidth=0.5)
        ax1.set_ylabel("Vm (mV)")
        ax1.set_title(f"{cell_type} response to IClamp")

        # IClamp plot
        ax2.plot(t, i, color='red', linewidth=2)
        ax2.set_xlabel("Time (ms)")
        ax2.set_ylabel("I (nA)")

        plt.tight_layout()

        # save
        filename = os.path.join(save_dir, f"{cell_type}.png")
        plt.savefig(filename, dpi=200)
        plt.close()
        print(f"Saved {filename}")


def _run_single_iclamp_tc(tc_cell, tstop=1800, delay=600, dur=5, amp=3, v_init=-70):
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

def run_iclamp_tc(save_dir, **kwargs):
    os.makedirs(save_dir, exist_ok=True)

    t, v, i, spikes = _run_single_iclamp_tc(**kwargs)

    # Plot
    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(10, 6), sharex=True, gridspec_kw={'height_ratios': [5, 1]}
    )
    ax1.plot(t, v, linewidth=0.5)
    ax1.set_ylabel("Vm (mV)")

    ax2.plot(t, i, color="red", linewidth=2)
    ax2.set_xlabel("Time (ms)")
    ax2.set_ylabel("I (nA)")

    plt.tight_layout()
    filename = os.path.join(save_dir, f"VPMCell.png")
    plt.savefig(filename, dpi=200)
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
    and saves histogram plots to CELL_TEMPLATES_PATH/outdir_name/*.jpg.
    """

    # collect .asc files
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

            # get 3D points
            dend_points = cell.get_dendrite_points()
            axon_points = cell.get_axon_points()

            # skip empty morphologies
            if len(dend_points) == 0 and len(axon_points) == 0:
                if verbose:
                    print(f"[skip] {template_name}: no axon/dendrite points")
                continue

            # compute z-values
            data = []
            labels = []
            if len(dend_points):
                data.append(dend_points[:, 2])
                labels.append("dendrite")
            if len(axon_points):
                data.append(axon_points[:, 2])
                labels.append("axon")

            # plot histograms
            fig, ax = plt.subplots(figsize=(6, 4))
            if mode == "stacked" and len(data) > 1:
                ax.hist(
                    data,
                    bins=bins,
                    range=bin_range,
                    density=density,
                    stacked=True,
                    label=labels,
                    alpha=0.8
                )
            else:
                for z_values, label in zip(data, labels):
                    ax.hist(
                        z_values,
                        bins=bins,
                        range=bin_range,
                        density=density,
                        alpha=0.5,
                        label=label,
                        histtype="stepfilled" if mode == "overlay" else "bar"
                    )

            ax.set_xlabel("z [µm]")
            ax.set_ylabel("Density" if density else "Count")
            ax.set_title(template_name)
            ax.legend()
            ax.grid(alpha=0.3)

            out_path = os.path.join(outdir, f"{template_name}.jpg")
            fig.savefig(out_path, dpi=200, bbox_inches="tight")
            plt.close(fig)

            saved.append((template_name, out_path))
            if verbose:
                print(f"[ok] saved {out_path}")

        except Exception as e:
            if verbose:
                print(f"[err] {template_name}: {e}")

    return saved