import os
import json

import config_templates.conf0_model_parameters as conf0
import source.src2_utils.ut1_path_config_parser as ut1

from source.src0_core.cr0_model_setup.m01_cx_cells import (cx01_tissue_structure as cx01,
                                                           cx02_population_composition as cx02,
                                                           cx03_instantiate_morphologies as cx03,
                                                           cx04_summary as cx04)

from source.src0_core.cr0_model_setup.m02_cx_cx_conn import (cxcx01_appositions_search as cxcx01,
                                                             cxcx02_appositions_prune as cxcx02,
                                                             cxcx03_summary as cxcx03)

from source.src0_core.cr0_model_setup.m03_tc_cells import (tc01_population_count as tc01)

from source.src0_core.cr0_model_setup.m04_tc_cx_conn import (tccx01_appositions_search as tccx01, tccx02_assign_cells_to_conn as tccx02)

from source.src0_core.cr0_model_setup.m06_prv_tc_conn import (prvtc01_spike_trains as prvtc01)

def run():
    # create paths structure
    paths = ut1.__main__()

    # cx01 - cx somata positions
    CX_SOMA_POS_PATH = paths["data"][conf0.MODEL_NAME]["setup"]["cells"]["cx"]["cx01"]
    for l in conf0.LAYER_COMP_TARGET.keys():
        l_pos = cx01.nrn_positions_cyl(
            layer_name=l,
            show_info=True)
        cx01.save_nrn_positions(
            nrn_pos=l_pos,
            save_dir=CX_SOMA_POS_PATH,
            file_name=f"{l}_soma_pos.csv")

    # cx02 - cx me-types
    CX_ME_ASSIGNED_PATH = paths["data"][conf0.MODEL_NAME]["setup"]["cells"]["cx"]["cx02"]
    CX_ME_SUMMARY = os.path.join(CX_ME_ASSIGNED_PATH, "me_comp_summary.json")
    cx02.assign_me_type(
        input_dir=CX_SOMA_POS_PATH,
        output_dir=CX_ME_ASSIGNED_PATH)

    # cx03 - cx morphology instantiation
    # TODO: PATH!
    CELL_TEMPLATES_PATH = "/home/mateusz-wawrzyniak/PycharmProjects/brcx_lfp_model/config_templates/cx_cell_templates"
    CX_FULL_POP_DIR = paths["data"][conf0.MODEL_NAME]["setup"]["cells"]["cx"]["cx03"]
    CX_FULL_POP_CSV = os.path.join(CX_FULL_POP_DIR, "cx03_full_pop.csv")

    cx_cells = {}
    for pop_file in sorted(os.listdir(CX_ME_ASSIGNED_PATH)):
        if not pop_file.startswith("cx02_"):
            continue
        path = CX_ME_ASSIGNED_PATH / pop_file
        pop = cx03.load_population(path)
        cx_cells.update(pop)
    all_cells = cx03.instantiate_morphologies(cx_cells, CELL_TEMPLATES_PATH)
    cx03.save_all_cells(all_cells, CX_FULL_POP_CSV)

    # cx04 - summary
    cx04.get_population_me_composition(
        pop_files_dir=CX_ME_ASSIGNED_PATH,
        save_path=CX_ME_SUMMARY
    )

    # cxcx01 - appositons search
    CXCX_SYNAPSE_SAVE_DIR = paths["data"][conf0.MODEL_NAME]["setup"]["conn"]["cx_cx"]["cxcx01"]
    CXCX_SYNAPSE_SAVE_PATH = os.path.join(CXCX_SYNAPSE_SAVE_DIR, "cxcx01_appositions.json")
    all_synapses, syn_cnt = cxcx01.find_pop_appositions(
        max_dist=conf0.MAX_APP_DIST,
        all_cells=all_cells)
    with open(CXCX_SYNAPSE_SAVE_PATH, "w") as f:
        json.dump(all_synapses, f, indent=2)


    # tccx01 - appositions search
    TCCX01_DIR = paths["data"][conf0.MODEL_NAME]["setup"]["conn"]["tc_cx"]["tccx01"]
    TCCX01_PATH = os.path.join(TCCX01_DIR, "tccx01_synapses.json")
    TCCX_BD_DIST = ut1.new_file_in_dir(paths["data"][conf0.MODEL_NAME]["visualizations"]["sim_dep"]["conn"]["tc_cx"], "tc_bouton_dist.jpg")
    z_bin_dict = tccx01.get_segments_binning(
        all_cells=all_cells,
        h_start=conf0.GLOBAL_Z_RANGE[1],
        h_stop=conf0.GLOBAL_Z_RANGE[0])
    bd_emp = tccx01.calc_tc_bd_dist(
        h_start=conf0.GLOBAL_Z_RANGE[1],
        h_stop=conf0.GLOBAL_Z_RANGE[0],
        save_path=TCCX_BD_DIST)
    tccx01.tc_synapse_sampling(
        z_bin_dict=z_bin_dict,
        bd_emp=bd_emp,
        save_path=TCCX01_PATH)

    del all_cells

    # cxcx02 - pruning
    CX_POP_PRE = os.path.join(paths["data"][conf0.MODEL_NAME]["setup"]["cells"]["cx"]["cx04"], "cx04_pop_pre.csv")
    PRUNE00_DIR = paths["data"][conf0.MODEL_NAME]["setup"]["conn"]["cx_cx"]["cxcx02"]["prune00"]

    PRUNE01_DIR = paths["data"][conf0.MODEL_NAME]["setup"]["conn"]["cx_cx"]["cxcx02"]["prune01"]
    PRUNED01_CXCX_PATH = os.path.join(PRUNE01_DIR, "prune01_cxcx.json")
    CX_POP_POST1 = os.path.join(paths["data"][conf0.MODEL_NAME]["setup"]["cells"]["cx"]["cx04"], "cx04_pop_post1.csv")

    PRUNE02_DIR = paths["data"][conf0.MODEL_NAME]["setup"]["conn"]["cx_cx"]["cxcx02"]["prune02"]
    PRUNED02_CXCX_PATH = os.path.join(PRUNE02_DIR, "prune02_cxcx.json")
    CX_POP_POST2 = os.path.join(paths["data"][conf0.MODEL_NAME]["setup"]["cells"]["cx"]["cx04"], "cx04_pop_post2.csv")

    PRUNE03_DIR = paths["data"][conf0.MODEL_NAME]["setup"]["conn"]["cx_cx"]["cxcx02"]["prune03"]
    PRUNED03_CXCX_PATH = os.path.join(PRUNE03_DIR, "prune03_cxcx.json")
    CX_POP_POST3 = os.path.join(paths["data"][conf0.MODEL_NAME]["setup"]["cells"]["cx"]["cx04"], "cx04_pop_post3.csv")

    # Calculating sm for m:m type prior to pruning
    cxcx02.calculate_sm(
        all_cells_csv=CX_FULL_POP_CSV,
        synapses_json=CXCX_SYNAPSE_SAVE_PATH,
        save_dir=PRUNE00_DIR
    )
    cxcx02.calc_b_int(
        full_population_csv=CX_FULL_POP_CSV,
        synapse_json=CXCX_SYNAPSE_SAVE_PATH,
        savepath=CX_POP_PRE
    )

    # Calculating retention rate for each m:m used in 01 general pruning
    # TODO: PATH!!
    BBP_PATHWAY_INFO = "/home/mateusz-wawrzyniak/PycharmProjects/brcx_lfp_model/source/src0_core/cr0_model_setup/m00_bbp_parameters/pathways_anatomy_factsheets_simplified.json"
    cxcx02.setup_01_pruning(
        model_type_sm_csv=os.path.join(PRUNE00_DIR, 'types_sm.csv'),
        goal_type_sm_json=BBP_PATHWAY_INFO,
        save_dir=PRUNE01_DIR
    )

    # Applying 01 general pruning
    cxcx02.apply_general_pruning_01(
        prune_01_params_csv=os.path.join(PRUNE01_DIR, "prune_01_params.csv"),
        init_synapse_json=CXCX_SYNAPSE_SAVE_PATH,
        output_json=PRUNED01_CXCX_PATH)

    # Calculating sm for m:m type post 01 general pruning
    cxcx02.calculate_sm(
        all_cells_csv=CX_FULL_POP_CSV,
        synapses_json=PRUNED01_CXCX_PATH,
        save_dir=PRUNE01_DIR
    )
    cxcx02.calc_b_int(
        full_population_csv=CX_FULL_POP_CSV,
        synapse_json=PRUNED01_CXCX_PATH,
        savepath=CX_POP_POST1
    )

    # Calculating cutoff threshold for each m:m used in 02 multisyn pruning
    cxcx02.setup_02_pruning(
        model_type_sm_csv=os.path.join(PRUNE01_DIR, 'types_sm.csv'),
        goal_type_sm_json=BBP_PATHWAY_INFO,
        save_dir=PRUNE02_DIR
    )

    # Applying 02 multisyn pruning
    cxcx02.apply_multisyn_pruning_02(
        prune_02_params_csv=os.path.join(PRUNE02_DIR, "prune_02_params.csv"),
        all_cells_csv=CX_FULL_POP_CSV,
        synapses_json=PRUNED01_CXCX_PATH,
        save_path=PRUNED02_CXCX_PATH
    )

    # Calculate interbouton interval for each cell as b_int = len_axon / n_syn
    cxcx02.calc_b_int(
        full_population_csv=CX_FULL_POP_CSV,
        synapse_json=PRUNED02_CXCX_PATH,
        savepath=CX_POP_POST2
    )

    # Calculating retention ratio for 03 plasticity-reserve pruning
    cxcx02.setup_03_pruning(
        examined_pop_csv=CX_POP_POST2,
        save_dir=PRUNE03_DIR
    )

    # Applying 03 plasticity reserve pruning
    cxcx02.apply_plastres_pruning_03(
        post02_synapse_json=PRUNED02_CXCX_PATH,
        prune_03_params_json=os.path.join(PRUNE03_DIR, "prune_03_params.json"),
        save_path=PRUNED03_CXCX_PATH
    )
    cxcx02.calc_b_int(
        full_population_csv=CX_FULL_POP_CSV,
        synapse_json=PRUNED03_CXCX_PATH,
        savepath=CX_POP_POST3
    )

    #cxcx03 - creating synapse summary
    CXCX_SUMMARY_PATH = os.path.join(PRUNE03_DIR, "cxcx_summary.json")
    cxcx03.create_synapse_summary(synapses_json=PRUNED03_CXCX_PATH,
                                  summary_path=CXCX_SUMMARY_PATH)

    # tc01 - deciding on tc cells count
    TC_CELLS_INIT_DIR = paths["data"][conf0.MODEL_NAME]["setup"]["cells"]["tc"]["tc01"]
    TC01_PATH = os.path.join(TC_CELLS_INIT_DIR, 'tc01_pop.csv')
    tc01.setup_vpm_pop(
        tccx_synapses_path = TCCX01_PATH,
        tc_cells_path = TC01_PATH
    )

    # tccx02 - assign tc cells to tccx synapses
    TCCX02_DIR = paths["data"][conf0.MODEL_NAME]["setup"]["conn"]["tc_cx"]["tccx02"]
    TCCX02_PATH = os.path.join(TCCX02_DIR, "tccx02_synapses.json")
    tccx02.assign_tc_synapses(
        preassignment_synapses_path = TCCX01_PATH,
        assigned_synapses_path = TCCX02_PATH,
        tc_pop_path = TC01_PATH)

    # prvtc01 - assigning prv spike trains and tc cells
    PRVTC01_DIR = paths["data"][conf0.MODEL_NAME]["setup"]["conn"]["prv_tc"]["prvtc01"]
    PRVTC01_PATH = os.path.join(PRVTC01_DIR, "prvtc01_synapses.json")

    prvtc01.generate_prv_spikes(
        stim_paradigm_type = conf0.STIM_PARADIGM_TYPE,
        stim_paradigm_subtype = conf0.STIM_PARADIGM_SUBTYPE,
        tc_cells_path = TC01_PATH,
        prvtc_save_path = PRVTC01_PATH)
