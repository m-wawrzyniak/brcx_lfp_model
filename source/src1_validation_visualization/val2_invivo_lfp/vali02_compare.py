import os
import pandas as pd
from scipy.signal import butter, sosfiltfilt
import h5py
import matplotlib.pyplot as plt
import re

from source.src2_utils.ut0_random_manager import np

def export_net_lfp_to_csv(lfp_hdf_path, csv_path):
    """
    Export net LFP from HDF5 to CSV.

    CSV format:
        time, site_0, site_1, ...
        t0,   v00,    v10, ...
        t1,   v01,    v11, ...
        ...

    Args:
    """
    with h5py.File(lfp_hdf_path, "r") as f:
        time = f["net_lfp/time"][:]               # shape (n_time,)
        lfp = f["net_lfp/lfp"][:]                 # shape (n_sites, n_time)
        site_ids = f["net_lfp/site_ids"][:]

    # Decode site IDs if they are stored as bytes
    site_ids = [sid.decode() if isinstance(sid, bytes) else str(sid) for sid in site_ids]

    # Build dataframe: transpose lfp to (n_time, n_sites)
    df = pd.DataFrame(lfp.T, columns=site_ids)
    df.insert(0, "time", time)

    # Save to CSV
    df.to_csv(csv_path, index=False, encoding='utf-8-sig')
    print(f"\t \t l01: Net LFP exported to {csv_path}")

def plot_single_lfp(lfp_signals_csv, save_path, *args):
    """
    Reads LFP data from CSV and plots all electrode sites vs. time on a single plot.

    Args:
        lfp_signals_csv: Path to CSV file with LFP data.
        save_path : Directory in which plot should be saved.
    """
    # Load CSV
    df = pd.read_csv(lfp_signals_csv)

    # Ensure the first column is time
    time = df.iloc[:, 0]
    site_columns = df.columns[1:]

    # Plot
    plt.figure(figsize=(10, 6))
    for col in site_columns:
        plt.plot(time, df[col], label=col)

    # Labels and legend
    plt.xlabel("Time (ms)")
    plt.ylabel("LFP (mV)")
    plt.title("Local Field Potentials over Time")
    plt.legend(loc="best")
    plt.tight_layout()

    # Save figure
    plt.savefig(save_path, dpi=300)
    plt.close()

def get_model_markers(paradigms, paradigm_type, paradigm_subtype):
    # Get the sequence for the given type and subtype
    try:
        stim_sequence = paradigms[paradigm_type][paradigm_subtype]
    except KeyError:
        raise ValueError(f"Invalid type '{paradigm_type}' or subtype '{paradigm_subtype}'")

    markers = []
    current_time = 0

    for name, (stim_type, duration) in stim_sequence.items():
        if stim_type in ('str', 'wk'):
            markers.append(current_time)
        current_time += duration

    return np.array(markers, dtype=float)

def get_emp_markers(file_path, dt):
    # Read the single row as a list of floats
    with open(file_path, 'r') as f:
        line = f.readline().strip()
    markers = np.array([float(x) for x in line.split(',')])*dt
    return markers

def trim_signals(model_lfp_csv, emp_lfp_csv, emp_marker_csv,
                 paradigm_type, paradigm_subtype,
                 invivo_csv_savepath,
                 model_csv_savepath,
                 markers_csv_savepath,
                 whisker_stim_paradigms,
                 paradigm_stand,
                 dt):
    """
    Trim model and empirical LFP signals independently, then adjust markers.
    """


    model_markers = get_model_markers(whisker_stim_paradigms, paradigm_type, paradigm_subtype)

    # --------------- Load data -----------------
    df_model = pd.read_csv(model_lfp_csv)
    df_emp = pd.read_csv(emp_lfp_csv)
    emp_markers = get_emp_markers(emp_marker_csv, dt=dt)

    # --------------- Get time window -----------------
    t_start, t_end = paradigm_stand[paradigm_type]['timerange']

    # --------------- Trim empirical -----------------
    emp_first = emp_markers[0]
    emp_mask = (df_emp['time'] >= emp_first + t_start) & (df_emp['time'] <= emp_first + t_end)
    df_emp_trimmed = df_emp.loc[emp_mask].copy()

    # Adjust empirical time to start at 0
    df_emp_trimmed['time'] = df_emp_trimmed['time'] - emp_first

    # --------------- Trim model -----------------
    mod_first = model_markers[0]
    mod_mask = (df_model['time'] >= mod_first + t_start) & (df_model['time'] <= mod_first + t_end)
    df_model_trimmed = df_model.loc[mod_mask].copy()

    # Adjust model time to start at 0
    df_model_trimmed['time'] = df_model_trimmed['time'] - mod_first

    ref_markers = emp_markers

    # Adjust so first marker is zero
    ref_markers_adj = ref_markers - ref_markers[0]

    # --------------- Save all -----------------

    df_emp_trimmed.to_csv(
        invivo_csv_savepath ,
        index=False
    )
    df_model_trimmed.to_csv(
        model_csv_savepath,
        index=False
    )
    pd.DataFrame({'time': ref_markers_adj}).to_csv(
        markers_csv_savepath,
        index=False
    )

    # --------------- Sanity check -----------------
    diff = np.max(np.abs((emp_markers-emp_markers[0]) - (model_markers-model_markers[0])))
    print(f"Max absolute difference between original marker sets: {diff:.3f} ms")

    if diff > 1.0:
        print("Warning: Empirical and model markers differ noticeably!")

    return df_emp_trimmed, df_model_trimmed, ref_markers_adj


def set_signal_baseline(lfp_csv, baseline_window,
                        savepath):
    df = pd.read_csv(lfp_csv)

    # Select rows in baseline window
    baseline_mask = (df['time'] >= baseline_window[0]) & (df['time'] <= baseline_window[1])
    baseline_data = df.loc[baseline_mask]

    # Compute mean baseline for each channel (exclude 'time')
    baseline_means = baseline_data.drop(columns='time').mean()

    # Subtract baseline from all signal columns
    for col in baseline_means.index:
        df[col] = df[col] - baseline_means[col]
    df.to_csv(savepath, index=False)

def _make_output_name(basename: str) -> str:
    name, ext = os.path.splitext(basename)
    repl = f"03_filt"
    if "02_based" in name:
        name = name.replace("02_based", repl)
    else:
        name = f"{name}_{repl}"
    return name + ext

def bandpass_filter_signal(lfp_csv, lp_cutoff: float, hp_cutoff: float | None,
                         order: int, fs, savepath) -> pd.DataFrame:
    df = pd.read_csv(lfp_csv)
    if 'time' not in df.columns:
        raise ValueError("CSV must contain a 'time' column (in milliseconds).")
    nyq = 0.5 * fs

    # Prepare normalized frequencies safely
    if hp_cutoff is None:
        # Low-pass
        high = min(lp_cutoff / nyq, 0.99)
        if not (0 < high < 1):
            raise ValueError(f"Low-pass cutoff invalid vs Nyquist: lp={lp_cutoff} Hz, nyq={nyq:.3f} Hz")
        sos = butter(order, high, btype='low', output='sos')
    else:
        # Band-pass
        low = max(hp_cutoff / nyq, 1e-6)
        high = min(lp_cutoff / nyq, 0.99)
        if not (0 < low < high < 1):
            raise ValueError(
                f"Band-pass cutoffs invalid: hp={hp_cutoff} Hz, lp={lp_cutoff} Hz, nyq={nyq:.3f} Hz "
                f"(normalized low={low:.4f}, high={high:.4f})"
            )
        sos = butter(order, [low, high], btype='band', output='sos')

    # Apply filter with SOS (stable) channel-wise
    filtered_df = df.copy()
    signal_cols = [c for c in df.columns if c != 'time']
    n_samples = len(df)
    if n_samples < 3 * (order + 1):
        print("⚠ Warning: very short signal for chosen order; edge effects may be significant.")

    for col in signal_cols:
        x = df[col].to_numpy(dtype=float)
        # Handle NaNs by simple interpolation to avoid sosfiltfilt failure
        if np.isnan(x).any():
            nans = np.isnan(x)
            not_nans = ~nans
            x[nans] = np.interp(np.flatnonzero(nans), np.flatnonzero(not_nans), x[not_nans])
        filtered_df[col] = sosfiltfilt(sos, x)

    # Save with your naming rule
    filtered_df.to_csv(savepath, index=False)
    print(f"Filtered signal saved to: {savepath}  |  fs={fs:.3f} Hz, nyq={nyq:.3f} Hz")

    return filtered_df

def manual_lfp_signal_processing(
    lfp_csv: str,
    save_path: str,
    trim_range: tuple[float, float],      # (t_start, t_end) in ms
    baseline_window: tuple[float, float], # (t0, t1) in ms
    lp_cutoff: float,
    hp_cutoff: float | None,
    order: int,
    fs: float
) -> pd.DataFrame:
    """
    Trim, baseline-correct, and bandpass filter an LFP CSV signal.

    Args:
        lfp_csv: Path to input CSV with 'time' column.
        save_path: Path to save processed CSV.
        trim_range: Tuple of (start_time, end_time) in ms.
        baseline_window: Time window to compute baseline mean (before trimming if desired).
        lp_cutoff: Low-pass cutoff frequency (Hz).
        hp_cutoff: High-pass cutoff frequency (Hz). None = low-pass only.
        order: Filter order.
        fs: Sampling frequency (Hz).
    Returns:
        Processed DataFrame.
    """
    # --- Load CSV ---
    df = pd.read_csv(lfp_csv)
    if 'time' not in df.columns:
        raise ValueError("CSV must contain a 'time' column.")

    # --- Trim signal ---
    t_start, t_end = trim_range
    mask = (df['time'] >= t_start) & (df['time'] <= t_end)
    df_trimmed = df.loc[mask].copy()
    df_trimmed['time'] = df_trimmed['time'] - t_start

    # --- Baseline correction ---
    baseline_mask = (df_trimmed['time'] >= baseline_window[0]) & (df_trimmed['time'] <= baseline_window[1])
    baseline_data = df_trimmed.loc[baseline_mask]
    baseline_means = baseline_data.drop(columns='time').mean()
    for col in baseline_means.index:
        df_trimmed[col] = df_trimmed[col] - baseline_means[col]

    # --- Bandpass / low-pass filter ---
    nyq = 0.5 * fs
    if hp_cutoff is None:
        high = min(lp_cutoff / nyq, 0.99)
        sos = butter(order, high, btype='low', output='sos')
    else:
        low = max(hp_cutoff / nyq, 1e-6)
        high = min(lp_cutoff / nyq, 0.99)
        sos = butter(order, [low, high], btype='band', output='sos')

    signal_cols = [c for c in df_trimmed.columns if c != 'time']
    for col in signal_cols:
        x = df_trimmed[col].to_numpy(dtype=float)
        if np.isnan(x).any():
            nans = np.isnan(x)
            not_nans = ~nans
            x[nans] = np.interp(np.flatnonzero(nans), np.flatnonzero(not_nans), x[not_nans])
        df_trimmed[col] = sosfiltfilt(sos, x)

    # --- Save processed CSV ---
    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    df_trimmed.to_csv(save_path, index=False)
    print(f"Processed LFP saved to: {save_path}")

    return df_trimmed

def plot_comparison_lfp(lfp_signals_csv, title, offset=0.010, save_path=None, make_dashed=True):
    """
    Plots LFP traces from CSV with vertical offsets and inverted site order
    (earlier sites appear at the top). Adds vertical dashed lines every 100 ms
    and a solid red line at 0 ms.

    Args:
        lfp_signals_csv: Path to CSV file with LFP data (columns: time, site_0, site_1, ...).
        offset (float): Vertical offset between traces.
        save_path: If provided, save figure instead of showing.
    """
    # Load CSV
    df = pd.read_csv(lfp_signals_csv)

    # First column is time
    time = df.iloc[:, 0]
    site_columns = df.columns[1:]

    # Sort site IDs numerically if possible (e.g. c_1, c_2, ..., c_10)
    def site_sort_key(s):
        m = re.search(r'(\d+)$', s)
        return int(m.group(1)) if m else s

    site_columns = sorted(site_columns, key=site_sort_key, reverse=True)  # reverse: smallest on top

    # Plot
    plt.figure(figsize=(12, 7))
    for i, col in enumerate(site_columns):
        trace = df[col] + i * offset
        plt.plot(time, trace, label=col)

    # --- Vertical dashed lines every 100 ms ---
    if make_dashed:
        t_min, t_max = time.iloc[0], time.iloc[-1]
        for tline in range(int(t_min // 100) * 100, int(t_max) + 100, 100):
            plt.axvline(tline, linestyle='--', color='gray', linewidth=0.5)

        # --- Solid reddish line at 0 ms ---
        plt.axvline(0, linestyle='-', color='#ff5555', linewidth=1.5)

    # Formatting
    plt.xlabel("Time [ms]")
    plt.ylabel("LFP (offset traces)")
    plt.title(title)
    plt.legend(loc="upper right", fontsize=8)
    plt.yticks([])  # hide y-axis ticks
    plt.tight_layout()

    # Save or show
    if save_path:
        plt.savefig(save_path, dpi=300)
        plt.close()
        print(f"\t \t l02: LFP plotted and saved to {save_path}")
    else:
        plt.show()

