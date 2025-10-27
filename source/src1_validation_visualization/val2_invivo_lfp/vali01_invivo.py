import h5py
import pandas as pd
import matplotlib.pyplot as plt
from source.src2_utils.ut0_random_manager import np

from source.src0_core.cr2_lfp_reconstruction.l01_electrode_setup.Electrode import Electrode

def explore_hdf5(file_path):
    """Recursively explore and print the structure of an HDF5 file."""

    def print_attrs(name, obj):
        indent = '  ' * (name.count('/') - 1)
        if isinstance(obj, h5py.Group):
            print(f"{indent} Group: {name}")
        elif isinstance(obj, h5py.Dataset):
            print(f"{indent} Dataset: {name} | shape={obj.shape} | dtype={obj.dtype}")
        else:
            print(f"{indent}❓ Unknown object: {name}")

        # Print attributes if any
        if obj.attrs:
            for key, val in obj.attrs.items():
                print(f"{indent}   └─ Attr: {key} = {val}")

    with h5py.File(file_path, 'r') as f:
        f.visititems(print_attrs)

def load_hdf5_data(file_path, paradigm, subtype):
    """
    Load 'epoch', 'marker_times', and 'std' datasets from a given paradigm/subtype.

    Parameters
    ----------
    file_path : str
        Path to the HDF5 file.
    paradigm : str
        Top-level group (e.g., 'train', 'prepulse', 'single').
    subtype : str
        Subgroup under paradigm (e.g., '40t', '100p', 'weak').

    Returns
    -------
    tuple of np.ndarray
        (epoch, marker_times, std)
    """
    with h5py.File(file_path, "r") as f:
        base_path = f"{paradigm}/{subtype}/average"

        epoch = np.array(f[f"{base_path}/epoch"])
        marker_times = np.array(f[f"{base_path}/marker_times"])
        std = np.array(f[f"{base_path}/std"])

    return epoch, marker_times, std

def get_significant_el_sites(av_sig: np.ndarray,
                             save_path,
                             el_topo_variant: str,
                             dvs_map: dict,
                             dt:float,
                             electrode_offset:float):
    """
    Extract signals from av_sig for the electrode sites defined in MasmanidisEl.site_topology.

    Parameters
    ----------
    av_sig : np.ndarray
        Array of shape (n_sites, n_samples) where each row is a site signal.
    save_path
        Path to save the resulting CSV.
    el_topo_variant : str
        Variant name for MasmanidisEl initialization.
    dvs_map : dict
        Dictionary mapping region letter to list of channel indices (1-based).
    fs : int
        Sampling frequency in Hz (default 4000).

    Returns
    -------
    pd.DataFrame
        DataFrame where first column is 'sample', second is 'time' (seconds),
        and remaining columns are site signals.
    """
    el = Electrode(z_offset=electrode_offset, topo_variant=el_topo_variant)

    # Ensure av_sig is (n_sites, n_samples) → transpose if needed
    if av_sig.shape[0] < av_sig.shape[1]:
        pass
    else:
        av_sig = av_sig.T

    n_samples = av_sig.shape[1]
    samples = np.arange(n_samples)
    times = samples * dt  # in ms

    data_dict = {
        #"sample": samples,
        "time": times
    }

    for site in el.site_topology:
        dvs_letter, dvs_index = site.split("_")
        dvs_index = int(dvs_index)
        channel_idx = dvs_map[dvs_letter][dvs_index] - 1
        data_dict[site] = av_sig[channel_idx, :]

    df = pd.DataFrame(data_dict)
    df.to_csv(save_path, index=False)
    return df


def plot_invivo_lfp(sig_df:pd.DataFrame, save_path):
    # Make sure save_dir exists if it's a directory
    plt.figure(figsize=(12, 6))

    for col in sig_df.columns:
        if col != "sample" and col != "time":
            plt.plot(sig_df["time"], sig_df[col], label=col, linewidth=0.8)

    plt.xlabel("Time [ms]")
    plt.ylabel("LFP amplitude")
    plt.legend(loc="upper right", fontsize="small", ncol=2)
    plt.tight_layout()

    plt.savefig(save_path, dpi=300)
    plt.close()

    print(f"In-vivo LFP plot saved to {save_path}")