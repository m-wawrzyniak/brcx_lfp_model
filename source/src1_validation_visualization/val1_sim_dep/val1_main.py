import os

from source.src1_validation_visualization.val1_sim_dep import (vald01_cells as vald01,
                                                               vald02_conns as vald02)

import source.src2_utils.ut1_path_config_parser as ut1
import config_templates.conf0_model_parameters as conf0

def run(model_name):
    MODEL_ROOT = os.path.join(conf0.ROOT, "data", model_name)
    MODEL_CONFIG_PATH = os.path.join("/home/mateusz-wawrzyniak/PycharmProjects/brcx_lfp_model/config_templates/",
                                     f"{model_name}_structure.json")
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
    CX_CELLS_CSV = paths["setup"]["cells"]["cx"]["cx04"]["cx04_pop_post2.csv"]
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
    CX_COMP_EMP_SAVEPATH = ut1.new_file_in_dir(paths["visualizations"]["sim_dep"]["cells"]["cx"], "cx_composition_emp.png")
    CX_COMP_THEO_SAVEPATH = ut1.new_file_in_dir(paths["visualizations"]["sim_dep"]["cells"]["cx"], "cx_composition_theo.png")
    vald01.cx_plot_me_composition_emp(
        json_path=CX_ME_SUMMARY_PATH,
        save_path=CX_COMP_EMP_SAVEPATH
    )
    vald01.cx_plot_me_composition_theo(
        layer_comp_params=conf0.LAYER_COMP_PARAMS,
        save_path=CX_COMP_THEO_SAVEPATH
    )

    CX_CELLS_CSV = paths["setup"]["cells"]["cx"]["cx04"]["cx04_pop_post2.csv"]
    CX_ME_Z_DIST_SAVEPATH = ut1.new_file_in_dir(paths["visualizations"]["sim_dep"]["cells"]["cx"], "cx_me_z_distribution.jpg")

    vald01.cx_somata_dist_z(
        cx_cells_path = CX_CELLS_CSV,
        save_path = CX_ME_Z_DIST_SAVEPATH
    )


    # vald02 - connections, simulation config dependent
    CXCX_CSV = paths["simulation"]["conn"]["cx_cx"]["cxcx_syns.csv"]
    CXCX_SYNMATRIX_EMP_SAVEPATH = ut1.new_file_in_dir(paths["visualizations"]["sim_dep"]["conn"]["cx_cx"], "syn_matrix_emp.jpg")
    vald02.conn_matrix_emp(
        csv_file= CXCX_CSV,
        save_path= CXCX_SYNMATRIX_EMP_SAVEPATH,
        layer_comp_params= conf0.LAYER_COMP_PARAMS,
    )
    CXCX_CONNMATRIX_EMP_SAVEPATH = ut1.new_file_in_dir(paths["visualizations"]["sim_dep"]["conn"]["cx_cx"], "conn_matrix_emp.jpg")
    vald02.conn_matrix_emp(
        csv_file= CXCX_CSV,
        save_path= CXCX_CONNMATRIX_EMP_SAVEPATH,
        layer_comp_params= conf0.LAYER_COMP_PARAMS,
        unique_connections=True
    )
    BBP_PATHWAY_JSON = os.path.join(conf0.ROOT, "source/src0_core/cr0_model_setup/m00_bbp_parameters/pathways_anatomy_factsheets_simplified.json")
    CXCX_CONNMATRIX_THEO_SAVEPATH = ut1.new_file_in_dir(paths["visualizations"]["sim_dep"]["conn"]["cx_cx"], "conn_matrix_theo.jpg")
    vald02.conn_matrix_theo(
        json_file= BBP_PATHWAY_JSON,
        save_path= CXCX_CONNMATRIX_THEO_SAVEPATH,
        layer_comp_params= conf0.LAYER_COMP_PARAMS
    )

    # tccx conn matricies
    TCCX_CSV = paths["simulation"]["conn"]["tc_cx"]["tccx_syns.csv"]
    TCCX_SYNMATRIX_EMP_SAVEPATH = ut1.new_file_in_dir(paths["visualizations"]["sim_dep"]["conn"]["tc_cx"], "tccx_syn_matrix_emp.jpg")
    vald02.conn_matrix_emp_tc(
        tccx_syns= TCCX_CSV,
        layer_comp_params= conf0.LAYER_COMP_PARAMS,
        savefig_path= TCCX_SYNMATRIX_EMP_SAVEPATH,
        unique_connections= False
    )
    TCCX_CONNMATRIX_EMP_SAVEPATH = ut1.new_file_in_dir(paths["visualizations"]["sim_dep"]["conn"]["tc_cx"], "tccx_conn_matrix_emp.jpg")
    vald02.conn_matrix_emp_tc(
        tccx_syns= TCCX_CSV,
        layer_comp_params= conf0.LAYER_COMP_PARAMS,
        savefig_path= TCCX_CONNMATRIX_EMP_SAVEPATH,
        unique_connections= True
    )

    # vald02- pruning distributions

    # interbouton distance
    CX_POP_PREPRUNE = paths["setup"]["cells"]["cx"]["cx04"]["cx04_pop_pre.csv"]
    HIST_PRE_SAVEPATH = ut1.new_file_in_dir(paths["visualizations"]["sim_dep"]["cells"]["cx"]["bouton_int"], "intb_dist_pre.jpg")
    CX_POP_POST1 = paths["setup"]["cells"]["cx"]["cx04"]["cx04_pop_post1.csv"]
    HIST_POST1_SAVEPATH = ut1.new_file_in_dir(paths["visualizations"]["sim_dep"]["cells"]["cx"]["bouton_int"], "intb_dist_post1.jpg")
    CX_POP_POST2 = paths["setup"]["cells"]["cx"]["cx04"]["cx04_pop_post2.csv"]
    HIST_POST2_SAVEPATH = ut1.new_file_in_dir(paths["visualizations"]["sim_dep"]["cells"]["cx"]["bouton_int"], "intb_dist_post2.jpg")
    CX_POP_POST3 = paths["setup"]["cells"]["cx"]["cx04"]["cx04_pop_post3.csv"]
    HIST_POST3_SAVEPATH = ut1.new_file_in_dir(paths["visualizations"]["sim_dep"]["cells"]["cx"]["bouton_int"], "intb_dist_post3.jpg")

    vald02.interbouton_int_histogram(
        cx_pop_csv= CX_POP_PREPRUNE,
        save_path= HIST_PRE_SAVEPATH
    )
    vald02.interbouton_int_histogram(
        cx_pop_csv= CX_POP_POST1,
        save_path= HIST_POST1_SAVEPATH
    )
    vald02.interbouton_int_histogram(
        cx_pop_csv= CX_POP_POST2,
        save_path= HIST_POST2_SAVEPATH
    )
    vald02.interbouton_int_histogram(
        cx_pop_csv= CX_POP_POST3,
        save_path= HIST_POST3_SAVEPATH
    )

    # synapses per connection
    SYNS_PREPRUNE_JSON = paths["setup"]["conn"]["cx_cx"]["cxcx01"]["cxcx01_appositions.json"]
    HIST_PRE_SAVEPATH = ut1.new_file_in_dir(paths["visualizations"]["sim_dep"]["conn"]["cx_cx"]["syns_per_conn"], "syns_per_conn_pre.jpg")
    SYNS_POST1_JSON = paths["setup"]["conn"]["cx_cx"]["cxcx02"]["prune01"]["prune01_cxcx.json"]
    HIST_POST1_SAVEPATH = ut1.new_file_in_dir(paths["visualizations"]["sim_dep"]["conn"]["cx_cx"]["syns_per_conn"], "syns_per_conn_post1.jpg")
    SYNS_POST2_JSON = paths["setup"]["conn"]["cx_cx"]["cxcx02"]["prune02"]["prune02_cxcx.json"]
    HIST_POST2_SAVEPATH = ut1.new_file_in_dir(paths["visualizations"]["sim_dep"]["conn"]["cx_cx"]["syns_per_conn"], "syns_per_conn_post2.jpg")
    SYNS_POST3_JSON = paths["setup"]["conn"]["cx_cx"]["cxcx02"]["prune03"]["prune03_cxcx.json"]
    HIST_POST3_SAVEPATH = ut1.new_file_in_dir(paths["visualizations"]["sim_dep"]["conn"]["cx_cx"]["syns_per_conn"], "syns_per_conn_post3.jpg")

    vald02.syn_per_conn_histogram(
        synapses_json= SYNS_PREPRUNE_JSON,
        savepath=HIST_PRE_SAVEPATH
    )
    vald02.syn_per_conn_histogram(
        synapses_json= SYNS_POST1_JSON,
        savepath=HIST_POST1_SAVEPATH
    )
    vald02.syn_per_conn_histogram(
        synapses_json= SYNS_POST2_JSON,
        savepath=HIST_POST2_SAVEPATH
    )
    vald02.syn_per_conn_histogram(
        synapses_json= SYNS_POST3_JSON,
        savepath=HIST_POST3_SAVEPATH
    )

    # connection prob. with distance
    PROB_CONN_PRE_SAVEPATH = ut1.new_file_in_dir(paths["visualizations"]["sim_dep"]["conn"]["cx_cx"]["conn_prob"], "conn_prob_pre.jpg")
    PROB_CONN_POST1_SAVEPATH = ut1.new_file_in_dir(paths["visualizations"]["sim_dep"]["conn"]["cx_cx"]["conn_prob"], "conn_prob_post1.jpg")
    PROB_CONN_POST2_SAVEPATH = ut1.new_file_in_dir(paths["visualizations"]["sim_dep"]["conn"]["cx_cx"]["conn_prob"], "conn_prob_post2.jpg")
    PROB_CONN_POST3_SAVEPATH = ut1.new_file_in_dir(paths["visualizations"]["sim_dep"]["conn"]["cx_cx"]["conn_prob"], "conn_prob_post3.jpg")

    vald02.plot_conn_prob_with_distance(
        cx_cells_csv=CX_POP_POST3,
        cxcx_synapses_json=SYNS_PREPRUNE_JSON,
        bin_size=50,
        max_dist=None,
        savepath=PROB_CONN_PRE_SAVEPATH
    )
    vald02.plot_conn_prob_with_distance(
        cx_cells_csv=CX_POP_POST3,
        cxcx_synapses_json=SYNS_POST1_JSON,
        bin_size=50,
        max_dist=None,
        savepath=PROB_CONN_POST1_SAVEPATH
    )
    vald02.plot_conn_prob_with_distance(
        cx_cells_csv=CX_POP_POST3,
        cxcx_synapses_json=SYNS_POST2_JSON,
        bin_size=50,
        max_dist=None,
        savepath=PROB_CONN_POST2_SAVEPATH
    )
    vald02.plot_conn_prob_with_distance(
        cx_cells_csv=CX_POP_POST3,
        cxcx_synapses_json=SYNS_POST3_JSON,
        bin_size=50,
        max_dist=None,
        savepath=PROB_CONN_POST3_SAVEPATH
    )
