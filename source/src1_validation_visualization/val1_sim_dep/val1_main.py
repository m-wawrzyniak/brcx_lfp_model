import os

from source.src1_validation_visualization.val1_sim_dep import (vald01_cells as vald01,
                                                               vald02_conns as vald02)

import source.src2_utils.ut1_path_config_parser as ut1
import config_templates.conf0_model_parameters as conf0


# Fetching model config file and loading paths
MODEL_NAME = "model_run_w0.1_2025-10-22_22-47-14"

MODEL_ROOT = os.path.join(conf0.ROOT, "data", MODEL_NAME)
MODEL_CONFIG_PATH = os.path.join("/home/mateusz-wawrzyniak/PycharmProjects/brcx_lfp_model/config_templates/",
                                 f"{MODEL_NAME}_structure.json")
ut1.folder_to_json_with_files(
    root_path=MODEL_ROOT,
    json_file=MODEL_CONFIG_PATH
)
paths = ut1.load_tree(config_file=MODEL_CONFIG_PATH, root=MODEL_ROOT)

# vald01 - cell, simulation config dependent
stimulation_paradigm = conf0.WHISKER_STIMULATION_PARADIGMS[conf0.STIM_PARADIGM_TYPE][conf0.STIM_PARADIGM_SUBTYPE]
CX_SPIKES_CSV = paths["recordings"]["cells"]["cx"]["cx_spikes.csv"]
TC_SPIKES_CSV = paths["recordings"]["cells"]["tc"]["tc_spikes.csv"]
PRV_SPIKES_CSV = paths["setup"]["conn"]["prv_tc"]["prvtc01"]["prvtc01_synapses.json"]
CX_CELLS_CSV = paths["setup"]["cells"]["cx"]["cx04"]["cx04_int_ext.csv"]
RASTER_SAVEPATH = ut1.new_file_in_dir(paths["visualizations"]["sim_dep"]["cells"]["cx"], "rasterplot.pdf")

vald01.plot_rasterplot(
    save_path = RASTER_SAVEPATH,
    paradigm = stimulation_paradigm,
    cx_spikes_csv = CX_SPIKES_CSV,
    vpm_spikes_csv = TC_SPIKES_CSV,
    prv_spikes_json = PRV_SPIKES_CSV,
    cell_pop_path = CX_CELLS_CSV
)

CX_ME_SUMMARY_PATH = paths["setup"]["cells"]["cx"]["cx02"]["me_comp_summary.json"]
CX_COMP_EMP_SAVEPATH = ut1.new_file_in_dir(paths["visualizations"]["sim_dep"]["cells"]["cx"], "cx_composition_emp.jpg")
CX_COMP_THEO_SAVEPATH = ut1.new_file_in_dir(paths["visualizations"]["sim_dep"]["cells"]["cx"], "cx_composition_theo.jpg")
vald01.cx_plot_me_composition_emp(
    json_path=CX_ME_SUMMARY_PATH,
    save_path=CX_COMP_EMP_SAVEPATH
)
vald01.cx_plot_me_composition_theo(
    layer_comp_params=conf0.LAYER_COMP_PARAMS, #TODO: im not sure, i think it gets it from a different module
    save_path=CX_COMP_THEO_SAVEPATH
)

CX_CELLS_CSV = paths["setup"]["cells"]["cx"]["cx04"]["cx04_int_ext.csv"]
CX_ME_Z_DIST_SAVEPATH = ut1.new_file_in_dir(paths["visualizations"]["sim_dep"]["cells"]["cx"], "cx_me_z_distribution.jpg")

vald01.cx_somata_dist_z(
    cx_cells_path = CX_CELLS_CSV,
    save_path = CX_ME_Z_DIST_SAVEPATH
)

