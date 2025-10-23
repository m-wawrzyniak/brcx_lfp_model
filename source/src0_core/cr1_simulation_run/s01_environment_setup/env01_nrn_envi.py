import os, sys
from pathlib import Path
import contextlib
from neuron import h

ROOT = Path("/home/mateusz-wawrzyniak/PycharmProjects/som_sens_processing")
MECH_DIR = Path("/home/mateusz-wawrzyniak/PycharmProjects/brcx_lfp_model/config_templates/nrn_mechanisms")
CELL_TEMP_DIR = "/home/mateusz-wawrzyniak/PycharmProjects/brcx_lfp_model/config_templates/cx_cell_templates"

@contextlib.contextmanager
def suppress_stdout_stderr():
    """
    While initializing the cells, a NEURON has a lot of stdout. This suppresses it.
    """
    with open(os.devnull, 'w') as fnull:
        fd_stdout = sys.__stdout__.fileno()
        fd_stderr = sys.__stderr__.fileno()

        def _redirect_fd(fd, target_fd):
            os.dup2(target_fd.fileno(), fd)

        orig_stdout = os.dup(fd_stdout)
        orig_stderr = os.dup(fd_stderr)

        try:
            _redirect_fd(fd_stdout, fnull)
            _redirect_fd(fd_stderr, fnull)
            yield
        finally:
            os.dup2(orig_stdout, fd_stdout)
            os.dup2(orig_stderr, fd_stderr)
            os.close(orig_stdout)
            os.close(orig_stderr)

# --- Ensure compiled mechanisms are loaded ---
def load_mechanisms(mech_dir):
    if not (mech_dir / "x86_64").exists():
        print(f"\tNo compiled mechanisms found in {mech_dir}/x86_64")
        print("\tCompiling MOD files...")
        os.system(f"cd {mech_dir} && nrnivmodl")
    h.nrn_load_dll(str(mech_dir / "x86_64" / "libnrnmech.so"))
    print(f"\tMechanisms loaded from {mech_dir}/x86_64")

    return h


def _load_template(template_dir:str, cx_cell_desc:str):
    """
    Loads template.hoc for specific cx_cell. Used to load (Markram, 2015) templates.
    Returns the template method directly from nrn.h

    Args:
        template_dir (str): Directory in which "template.hoc" for specific cell type can be found.
        cx_cell_desc (str): Cx cell me-type e.g. 'L4_SS_cAD'

    Returns:
        h.<method> : Cell template loaded into nrn.h.
    """
    temp_path = os.path.join(template_dir, "template.hoc")
    with suppress_stdout_stderr():
        h.load_file(temp_path)
    for method in dir(h):
        if cx_cell_desc in method:
            print(f"\t OK. Template for {cx_cell_desc} loaded.")
            return getattr(h, method)
    raise TypeError(f'ERR. Template for {cx_cell_desc} loading failed')

def load_all_templates(cx_cell_templates_dir, just_mtype:bool=True) -> dict:
    """
    Loads all CX cells templates located in cx_cell_templates_dir. The template files are structured like the original
    (Markram, 2015) templates provided.

    Args:
        cx_cell_templates_dir : Directory with all cx cell templates.

    Returns:
        dict: <cell_mtype: h.template>
    """

    temps = {}
    for cell_full_name in os.listdir(cx_cell_templates_dir):
        temp_path = os.path.join(cx_cell_templates_dir, cell_full_name)
        cell_mtype = "_".join(cell_full_name.split('_')[:2])
        if just_mtype:
            temps[cell_mtype] = _load_template(temp_path, cx_cell_desc=cell_mtype)
        else:
            cell_metype = "_".join(cell_full_name.split('_')[:3])
            temps[cell_metype] = _load_template(temp_path, cx_cell_desc=cell_mtype)

    return temps

# --- Initialize NEURON environment ---
def setup_neuron_env(mech_dir = MECH_DIR):
    print("[env01] Creating NEURON environment.")

    load_mechanisms(mech_dir)
    h.load_file("nrngui.hoc")
    h.load_file("import3d.hoc")
    # TODO: Change cell paths later
    cell_temps = load_all_templates(cx_cell_templates_dir=CELL_TEMP_DIR)

    print("[env01] SUCCESS: NEURON environment created.")
    return h, cell_temps
