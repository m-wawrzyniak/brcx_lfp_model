import csv
from neuron import h

from source.src2_utils.ut0_random_manager import np

from source.src0_core.cr1_simulation_run.s02_cells_init.CxCell import CxCell
from source.src0_core.cr1_simulation_run.s02_cells_init.VPMCell import VPMCell

def save_synapses_currents_csv(recordings:dict[str, dict], time_vector:h.Vector, save_filepath:str):
    """
    Saves recorded synaptic currents to a CSV file.

    Args:
        recordings (dict): Output from record_synapses_currents().
        time_vector (h.Vector): Used during simulation as time vector.
        save_filepath (str): Path where resulting *.csv file should be saved.

    """
    syn_ids = sorted(recordings.keys())
    num_points = len(time_vector)

    with open(save_filepath, mode='w', newline='') as file:
        writer = csv.writer(file)

        # Header
        header = ["time"] + syn_ids
        writer.writerow(header)

        # Rows
        for i in range(num_points):
            row = [time_vector[i]]
            for syn_id in syn_ids:
                row.append(recordings[syn_id]['vec'][i])
            writer.writerow(row)

    print(f"\t\t [INFO] Synaptic currents saved to {save_filepath}")

def save_cell_v_csv(cell_v_records:dict, time_vector:h.Vector, save_filepath:str):
    """
    Save all recorded voltages from `cell_v_records` to CSV.
    Rows are time points (from `time_vector`), columns are cell IDs.

    Args:
        cell_v_records (dict): <cell_id: h.Vector()>
        time_vector (h.Vector): Recording h._ref_t
        save_filepath (str): Output CSV filename
    """
    if not cell_v_records:
        print("No voltages recorded.")
        return

    if time_vector is None or len(time_vector) == 0:
        print("No time vector provided.")
        return

    times = np.array(time_vector)
    all_data = {cell_id: np.array(vec) for cell_id, vec in cell_v_records.items()}
    max_len = len(times)

    # Pad shorter vectors (shouldn't be needed if all recorded properly)
    for cell_id in all_data:
        if len(all_data[cell_id]) < max_len:
            all_data[cell_id] = np.pad(all_data[cell_id], (0, max_len - len(all_data[cell_id])), constant_values=np.nan)

    # Write to CSV
    with open(save_filepath, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['time(ms)'] + list(all_data.keys()))
        for i in range(max_len):
            row = [times[i]] + [all_data[cell_id][i] for cell_id in all_data]
            writer.writerow(row)

def save_cell_spikes_csv(cells:dict[str, CxCell|VPMCell], save_filepath:str):
    """
    Fetches spike times of each Cell provided in the dictionary and saves them in an individual *.csv file.
    Args:
        cells (dict): <cell_id, Cell>
        save_filepath (str): Where the *.csv should be saved.
    """
    with open(save_filepath, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["cell_id", "spikes"])  # header

        for cell_id, cell_obj in cells.items():
            times = cell_obj.get_spike_times()
            time_str = ' '.join(f"{t:.3f}" for t in times)
            writer.writerow([cell_id, time_str])
