import os

from source.src0_core.cr1_simulation_run.s01_environment_setup import (env01_nrn_envi as env01,
                                                                       env02_nrn_params as env02)
from source.src0_core.cr1_simulation_run.s02_cells_init.VPMCell import VPMCell

from source.src1_validation_visualization.val0_sim_indep import (vali01_cells as vali01,
                                                                 vali02_conns as vali02)

import source.src2_utils.ut1_path_config_parser as ut1
import config_templates.conf0_model_parameters as conf0


# Fetching model config file and loading paths
MODEL_NAME = "model_run_w0.1_2025-10-23_18-16-17"

MODEL_ROOT = os.path.join(conf0.ROOT, "data", MODEL_NAME)
MODEL_CONFIG_PATH = os.path.join("/home/mateusz-wawrzyniak/PycharmProjects/brcx_lfp_model/config_templates/",
                                 f"{MODEL_NAME}_structure.json")
ut1.folder_to_json_with_files(
    root_path=MODEL_ROOT,
    json_file=MODEL_CONFIG_PATH
)
paths = ut1.load_tree(config_file=MODEL_CONFIG_PATH, root=MODEL_ROOT)


# env01 - nrn envi setup
h, cx_templates = env01.setup_neuron_env()
# env02 - nrn params setup
h = env02.set_simulation_params(h)

"""
# vali01 - cell, simulation config independent val/vis
# cx_iclamp
CX_ICLAMP_SAVEDIR = os.path.join(paths["visualizations"]["sim_indep"]["cells"]["cx"], "cx_iclamp")
vali01.run_iclamp_cx(
    save_dir=CX_ICLAMP_SAVEDIR,
    cell_templates = cx_templates
)
#tc_iclamp
TC_ICLAMP_SAVEDIR = os.path.join(paths["visualizations"]["sim_indep"]["cells"]["tc"], "tc_iclamp")
tc_cell = VPMCell(cell_id="VPM_1")
vali01.run_iclamp_tc(
    save_dir=TC_ICLAMP_SAVEDIR,
    tc_cell=tc_cell
)
# cx morpho distribution along z-axis
CX_Z_DIST_SAVEDIR = os.path.join(paths["visualizations"]["sim_indep"]["cells"]["cx"], "cx_morpho_distribution")
vali01.plot_dend_axon_z_dist(
    cell_templates_path=env01.CELL_TEMP_DIR,
    outdir=CX_Z_DIST_SAVEDIR
)
"""

# vali02 - conns, simulation config independent val/vis
CXCX_SYNAPSES_CSV = paths["simulation"]["conn"]["cx_cx"]["cxcx_syns.csv"]
CX_POP_CSV = paths["setup"]["cells"]["cx"]["cx04"]["cx04_pop_post3.csv"]
CXCX_SYNAPSES_EPSP_SAVEDIR = paths["visualizations"]["sim_indep"]["conn"]["cx_cx"]["me_pair_epsp"]
"""
vali02.cxcx_conn_epsp_check(
    synapses_csv_path= CXCX_SYNAPSES_CSV,
    cell_pop_csv= CX_POP_CSV,
    cell_templates= cx_templates,
    save_dir= CXCX_SYNAPSES_EPSP_SAVEDIR
)

# TODO: they cant go both, first this then rerun
vali02.cxcx_conn_epsp_plot(
    me_pair_epsp_dir=CXCX_SYNAPSES_EPSP_SAVEDIR
)
"""
TCCX_SYNAPSES_CSV = paths["simulation"]["conn"]["tc_cx"]["tccx_syns.csv"]
TCCX_SYNAPSES_EPSP_SAVEDIR = paths["visualizations"]["sim_indep"]["conn"]["tc_cx"]["me_pair_epsp"]

vali02.tccx_conn_epsp_check(
    tccx_synapses_csv_path=TCCX_SYNAPSES_CSV,
    cell_pop_csv=CX_POP_CSV,
    cell_templates=cx_templates,
    save_dir=TCCX_SYNAPSES_EPSP_SAVEDIR
)

vali02.conn_epsp_plot(
    me_pair_epsp_dir=TCCX_SYNAPSES_EPSP_SAVEDIR
)