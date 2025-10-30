import os

from source.src0_core.cr2_lfp_reconstruction.l01_electrode_setup.Electrode import Electrode
from source.src0_core.cr2_lfp_reconstruction.l02_lfp_computation import (lfp01_compute as lfp01,
                                                                         lfp02_visualize as lfp02)

import source.src2_utils.ut1_path_config_parser as ut1
import config_templates.conf0_model_parameters as conf0
import config_templates.conf02_lfp_parameters as conf02

def run():
    MODEL_ROOT = os.path.join(conf0.ROOT, "data", conf0.MODEL_NAME)
    MODEL_CONFIG_PATH = os.path.join("/home/mateusz-wawrzyniak/PycharmProjects/brcx_lfp_model/config_templates/", f"{conf0.MODEL_NAME}_structure.json")
    ut1.folder_to_json_with_files(
        root_path=MODEL_ROOT,
        json_file=MODEL_CONFIG_PATH
    )
    paths = ut1.load_tree(config_file=MODEL_CONFIG_PATH, root=MODEL_ROOT)

    # Electrode - creating electrode topology
    el = Electrode(topo_variant=conf02.ELECTRODE_VARIANT, z_offset=conf02.ELECTRODE_Z_OFFSET, cross_species_scale=conf02.CROSS_SPECIES_SCALE)

    # lfp01 - computing lfp
    IMEM_PATH = paths["recordings"]["lfp"]["raw"]["cx_imem.hdf"]
    RECONSTRUCTED_LFP_PATH = os.path.join(paths["recordings"]["lfp"]["reconstructed"], "component_lfp.hdf")
    lfp01.reconstruct_lfp_from_hdf(
        imem_hdf_path=IMEM_PATH,
        electrode=el,
        out_path=RECONSTRUCTED_LFP_PATH
    )
    lfp01.compute_and_save_net_lfp(
        lfp_hdf_path=RECONSTRUCTED_LFP_PATH
    )

    # lfp02 - rudimentary visualization
    LFP_COMP_VIS_SAVEDIR = paths["visualizations"]["lfp"]["component_lfp"]
    LFP_NET_VIS_SAVEPATH = os.path.join(paths["visualizations"]["lfp"]["net_lfp"], "net_lfp.jpg")
    lfp02.export_all_cells_lfp_plots(
        lfp_hdf_path=RECONSTRUCTED_LFP_PATH,
        imem_hdf_path=IMEM_PATH,
        electrode=el,
        out_dir=LFP_COMP_VIS_SAVEDIR,
    )
    stimulation_paradigm = conf0.WHISKER_STIMULATION_PARADIGMS[conf0.STIM_PARADIGM_TYPE][conf0.STIM_PARADIGM_SUBTYPE]
    lfp02.plot_net_lfp(
        lfp_hdf_path=RECONSTRUCTED_LFP_PATH,
        save_path=LFP_NET_VIS_SAVEPATH,
        stim_paradigm=stimulation_paradigm,
        offset=0.01)
