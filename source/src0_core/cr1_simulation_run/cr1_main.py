import os

from source.src0_core.cr1_simulation_run.s01_environment_setup import (env01_nrn_envi as env01, env02_nrn_params as env02)
from source.src0_core.cr1_simulation_run.s02_cells_init import (initc01_cx as initc01,
                                                                initc02_tc as initc02)
from source.src0_core.cr1_simulation_run.s03_conn_init import (inits01_cxcx as inits01,
                                                               inits02_tccx as inits02)

import source.src3_utils.ut1_path_config_parser as ut1
import config_templates.conf0_model_parameters as conf0

MODEL_ROOT = os.path.join(conf0.ROOT, "data", conf0.MODEL_NAME)
MODEL_CONFIG_PATH = os.path.join("/home/mateusz-wawrzyniak/PycharmProjects/brcx_lfp_model/config_templates/", f"{conf0.MODEL_NAME}_structure.json")
ut1.folder_to_json_with_files(
    root_path=MODEL_ROOT,
    json_file=MODEL_CONFIG_PATH
)
paths = ut1.load_tree(config_file=MODEL_CONFIG_PATH, root=MODEL_ROOT)

# env1 - loading nrn mechanisms and cell templates
h, cx_templates = env01.setup_neuron_env()
# env2 - setting global nrn parameters
h = env02.set_simulation_params(h)

# initc01 - initialize cx cells
CX_ME_SUMMARY = paths["setup"]["cells"]["cx"]["cx02"]["me_comp_summary.json"]
CX_CELLS = paths["setup"]["cells"]["cx"]["cx04"]["cx04_int_ext.csv"]
cx_cells = initc01.init_cx_cells(
    cx_cells_csv=CX_CELLS,
    temp_dict=cx_templates,
    me_comp_summary_path=CX_ME_SUMMARY
)
for cell in cx_cells.values():
    cell.populate_segs_meta_for_synapses()

# inits01 - initialize cxcx synapses
CXCX_SUMMARY_PATH = paths["setup"]["conn"]["cx_cx"]["cxcx02"]["prune03"]["cxcx_summary.json"]
CXCX_CONN_PATH = paths["setup"]["conn"]["cx_cx"]["cxcx02"]["prune03"]["prune03_cxcx.json"]
# TODO: GLOBAL PATH
BBP_PHYSIO_DATA = "/home/mateusz-wawrzyniak/PycharmProjects/brcx_lfp_model/source/src0_core/cr0_model_setup/m00_bbp_parameters/pathways_physiology_factsheets_simplified.json"
CXCX_SYN_PATH = os.path.join(paths["simulation"]["conn"]["cx_cx"], "cxcx_syns.csv")
cxcx_parameters = inits01.get_synapse_parameters(
    cx_syns_summary_csv=CXCX_SUMMARY_PATH,
    bbp_physio_file=BBP_PHYSIO_DATA
)
cxcx_conn = inits01.create_cxcx_synapses(
    cx_cells=cx_cells,
    syn_params=cxcx_parameters,
    cxcx_synapses_data=CXCX_CONN_PATH
)
inits01.save_synapses_to_csv(
    synapse_map=cxcx_conn,
    save_path=CXCX_SYN_PATH
)

# initc02 - initalize tc cells
TC_CELLS = paths["setup"]["cells"]["tc"]["tc01"]["tc01_pop.csv"]
tc_cells = initc02.init_vpm_cells(
    vpm_cells_path=TC_CELLS
)

# inits02 - initialize tccx synapses
TCCX_CONN_PATH = paths["setup"]["conn"]["tc_cx"]["tccx02"]["tccx02_synapses.json"]
TCCX_SYN_PATH = os.path.join(paths["simulation"]["conn"]["tc_cx"], "tccx_syns.csv")
tccx_syns = inits02.create_tccx_synapses(
    tc_cells=tc_cells,
    cx_cells=cx_cells,
    tccx_synapses_data=TCCX_CONN_PATH
)
inits02.save_synapses_to_csv(tccx_syns, save_path=TCCX_SYN_PATH)

# TODO: p02!!!