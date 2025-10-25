import os

from source.src0_core.cr1_simulation_run.s01_environment_setup import (env01_nrn_envi as env01, env02_nrn_params as env02)
from source.src0_core.cr1_simulation_run.s02_cells_init import (initc01_cx as initc01,
                                                                initc02_tc as initc02)
from source.src0_core.cr1_simulation_run.s03_conn_init import (inits01_cxcx as inits01,
                                                               inits02_tccx as inits02,
                                                               inits03_prvtc as inits03)

from source.src0_core.cr1_simulation_run.s04_simulation_handling import (sim01_run as sim01)

from source.src0_core.cr1_simulation_run.s05_data_recording_saving import (rec01_setup_record as rec01,
                                                                           rec02_record_save as rec02,
                                                                           rec03_lfp_reconstruction_data as rec03)

import source.src2_utils.ut1_path_config_parser as ut1
import config_templates.conf0_model_parameters as conf0

def run():
    MODEL_ROOT = os.path.join(conf0.ROOT, "data", conf0.MODEL_NAME)
    MODEL_CONFIG_PATH = os.path.join("/home/mateusz-wawrzyniak/PycharmProjects/brcx_lfp_model/config_templates/", f"{conf0.MODEL_NAME}_structure.json")
    ut1.folder_to_json_with_files(
        root_path=MODEL_ROOT,
        json_file=MODEL_CONFIG_PATH
    )
    paths = ut1.load_tree(
        config_file=MODEL_CONFIG_PATH,
        root=MODEL_ROOT
    )

    # env1 - loading nrn mechanisms and cell templates
    h, cx_templates = env01.setup_neuron_env()
    # env2 - setting global nrn parameters
    h = env02.set_simulation_params(h)

    # initc01 - initialize cx cells
    CX_ME_SUMMARY = paths["setup"]["cells"]["cx"]["cx02"]["me_comp_summary.json"]
    CX_CELLS = paths["setup"]["cells"]["cx"]["cx04"]["cx04_pop_post2.csv"]
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
    BBP_PHYSIO_DATA = os.path.join(conf0.ROOT, "source/src0_core/cr0_model_setup/m00_bbp_parameters/pathways_physiology_factsheets_simplified.json")
    CXCX_SYN_PATH = os.path.join(paths["simulation"]["conn"]["cx_cx"], "cxcx_syns.csv")
    cxcx_parameters = inits01.get_synapse_parameters(
        cx_syns_summary_csv=CXCX_SUMMARY_PATH,
        bbp_physio_file=BBP_PHYSIO_DATA
    )
    cxcx_syns = inits01.create_cxcx_synapses(
        cx_cells=cx_cells,
        syn_params=cxcx_parameters,
        cxcx_synapses_data=CXCX_CONN_PATH
    )
    inits01.save_synapses_to_csv(
        synapse_map=cxcx_syns,
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
    inits02.save_synapses_to_csv(
        synapse_map=tccx_syns,
        save_path=TCCX_SYN_PATH
    )

    # inits03 - initialize prvtc synapses
    PRVTC_CONN_PATH = paths["setup"]["conn"]["prv_tc"]["prvtc01"]["prvtc01_synapses.json"]
    prvtc_syns = inits03.create_prvtc_synapses(
        tc_cells=tc_cells,
        prv_json_path=PRVTC_CONN_PATH
    )

    # rec01 - recording setup
    cx_soma_v = rec01.record_soma_v(
        cell_label="cx",
        cells=cx_cells
    )
    tc_soma_v = rec01.record_soma_v(
        cell_label="cx",
        cells=tc_cells
    )

    cxcx_i = rec01.record_synapses_currents(
        syn_label="cxcx",
        synapses=cxcx_syns
    )
    tccx_i = rec01.record_synapses_currents(
        syn_label="tccx",
        synapses=tccx_syns
    )
    prvtc_i = rec01.record_synapses_currents(
        syn_label="prvcx",
        synapses=prvtc_syns
    )

    rec01.record_cell_spikes(
        cell_label="cx",
        cells=cx_cells
    )
    rec01.record_cell_spikes(
        cell_label="tc",
        cells=tc_cells
    )

    # rec03 - lfp data recording setup
    rec03.start_recording_all_imem(cells_dict=cx_cells)

    # RUN SIMULATION
    t_vec = sim01.run_simulation()

    # rec03 - save lfp data
    IMEM_SAVEPATH = os.path.join(paths["recordings"]["lfp"]["raw"], "cx_imem.hdf")
    rec03.save_all_imem_h5(
        cells_dict=cx_cells,
        filename=IMEM_SAVEPATH
    )

    # rec02 - save other data
    CX_V_SAVEPATH = os.path.join(paths["recordings"]["cells"]["cx"], "cx_v.csv")
    TC_V_SAVEPATH = os.path.join(paths["recordings"]["cells"]["tc"], "tc_v.csv")
    rec02.save_cell_v_csv(
        cell_label="cx",
        cell_v_records=cx_soma_v,
        time_vector=t_vec,
        save_filepath=CX_V_SAVEPATH
    )
    rec02.save_cell_v_csv(
        cell_label="tc",
        cell_v_records=tc_soma_v,
        time_vector=t_vec,
        save_filepath=TC_V_SAVEPATH
    )

    CX_SP_SAVEPATH = os.path.join(paths["recordings"]["cells"]["cx"], "cx_spikes.csv")
    TC_SP_SAVEPATH = os.path.join(paths["recordings"]["cells"]["tc"], "tc_spikes.csv")
    rec02.save_cell_spikes_csv(
        cell_label="cx",
        cells=cx_cells,
        save_filepath=CX_SP_SAVEPATH
    )
    rec02.save_cell_spikes_csv(
        cell_label="tc",
        cells=tc_cells,
        save_filepath=TC_SP_SAVEPATH
    )

    CXCX_I_SAVEPATH = os.path.join(paths["recordings"]["conn"]["cx_cx"], "cxcx_i.csv")
    TCCX_I_SAVEPATH = os.path.join(paths["recordings"]["conn"]["tc_cx"], "tccx_i.csv")
    PRVTC_I_SAVEPATH = os.path.join(paths["recordings"]["conn"]["prv_tc"], "prvtc_i.csv")
    rec02.save_synapses_currents_csv(
        syn_label="cxcx",
        recordings=cxcx_i,
        time_vector=t_vec,
        save_filepath=CXCX_I_SAVEPATH
    )
    rec02.save_synapses_currents_csv(
        syn_label="tccx",
        recordings=tccx_i,
        time_vector=t_vec,
        save_filepath=TCCX_I_SAVEPATH
    )
    rec02.save_synapses_currents_csv(
        syn_label="prvtc",
        recordings=prvtc_i,
        time_vector=t_vec,
        save_filepath=PRVTC_I_SAVEPATH
    )
