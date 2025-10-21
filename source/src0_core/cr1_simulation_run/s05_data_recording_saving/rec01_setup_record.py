from neuron import h

from source.src0_core.cr1_simulation_run.s03_conn_init.CxSynapse import CxSynapse
from source.src0_core.cr1_simulation_run.s03_conn_init.TcCxSynapse import TcCxSynapse
from source.src0_core.cr1_simulation_run.s03_conn_init.PrvTcSynapse import PrvTcSynapse

from source.src0_core.cr1_simulation_run.s02_cells_init.CxCell import CxCell
from source.src0_core.cr1_simulation_run.s02_cells_init.VPMCell import VPMCell

def record_synapses_currents(synapses:dict[str, CxSynapse|TcCxSynapse|PrvTcSynapse]) -> dict[str, dict]:
    """
    Sets up recording of synaptic currents for a dictionary of synapse objects (CxSynapse or VPMSynapse).

    Args:
        synapses (dict): Dictionary of synapses meant to record their global _ref_i.

    Returns:
        dict : <syn_id: {'vec': h.Vector, 'ref': _ref_i}>
    """
    recordings = {}

    for syn_id, syn_obj in synapses.items():
        v = h.Vector()
        try:
            v.record(syn_obj.hsyn._ref_i)
            recordings[syn_id] = {'vec': v, 'ref': syn_obj.hsyn._ref_i}
        except AttributeError:
            print(f"\t\t [WARN] Synapse {syn_id} has no attribute '_ref_i'. Skipping.")
            continue

    return recordings

def record_soma_v(cells:dict[str, CxCell | VPMCell]) -> dict:
    """
    For each cell in the `cells` list or dict, record membrane voltage at soma(0.5).
    Handles both `soma(0.5)` and `soma[0](0.5)` formats.

    Args:
        cells (dict): Cortical cells which should record their membrane potential during simulation.

    Returns:
        dict : <cell_id: h.Vector>. Cell ids and their respective membrane potential vectors.
    """
    cell_v_records = {}
    for cell_id, cell in cells.items():
        v_vec = h.Vector()

        try:
            # Try the list-like access
            v_vec.record(cell.h_cell.soma[0](0.5)._ref_v)
        except:
            # Fallback for single-section soma
            v_vec.record(cell.soma(0.5)._ref_v)

        cell_v_records[cell_id] = v_vec

    return cell_v_records

def record_cell_spikes(cells:dict[str, CxCell|VPMCell]):
    """
    Sets internal spike detector of each CxCell to recording.

    Args:
          (dict): <cell_id, CxCell()>
    """
    for cx_c_id, cx_c_obj in cells.items():
        cx_c_obj.setup_spike_detector()