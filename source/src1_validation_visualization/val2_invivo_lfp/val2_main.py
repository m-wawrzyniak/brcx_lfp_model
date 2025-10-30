import os

from source.src2_utils.ut0_random_manager import np

from source.src1_validation_visualization.val2_invivo_lfp import (vali01_invivo as vali01,
                                                                  vali02_compare as vali02)

import source.src2_utils.ut1_path_config_parser as ut1

import config_templates.conf0_model_parameters as conf0
import config_templates.conf01_simulation_parameters as conf01
import config_templates.conf02_lfp_parameters as conf02


def run(model_name):

    MODEL_ROOT = os.path.join(conf0.ROOT, "data", model_name)
    MODEL_CONFIG_PATH = os.path.join("/home/mateusz-wawrzyniak/PycharmProjects/brcx_lfp_model/config_templates/",
                                     f"{model_name}_structure.json")
    ut1.folder_to_json_with_files(
        root_path=str(MODEL_ROOT),
        json_file=MODEL_CONFIG_PATH
    )
    paths = ut1.load_tree(config_file=MODEL_CONFIG_PATH, root=MODEL_ROOT)

    # BTBR
    av_sig, marker_times, std_sig = vali01.load_hdf5_data(
        file_path="/home/mateusz-wawrzyniak/PycharmProjects/brcx_lfp_model/source/src1_validation_visualization/val2_invivo_lfp/invivo_data/BTBRv30_mean_clean.hdf",
        paradigm=conf0.STIM_PARADIGM_TYPE,
        subtype=conf0.STIM_PARADIGM_SUBTYPE_INVIVO
    )
    BTBR_DF = ut1.new_file_in_dir(paths["recordings"]["invivo_lfp"]["invivo_exported"], "btbr_df.csv")
    BTBR_MARKERS = ut1.new_file_in_dir(paths["recordings"]["invivo_lfp"]["invivo_exported"], "btbr_markers.csv")
    btbr_sig_df = vali01.get_significant_el_sites(
        av_sig=av_sig,
        save_path=BTBR_DF,
        el_topo_variant = conf02.ELECTRODE_VARIANT,
        dvs_map = conf02.DVS_MAP,
        dt = conf01.DT,
        electrode_offset=conf02.ELECTRODE_Z_OFFSET
    )
    np.savetxt(BTBR_MARKERS, marker_times, delimiter=",")

    #C57
    C57_DF = ut1.new_file_in_dir(paths["recordings"]["invivo_lfp"]["invivo_exported"], "c57_df.csv")
    C57_MARKERS = ut1.new_file_in_dir(paths["recordings"]["invivo_lfp"]["invivo_exported"], "c57_markers.csv")
    av_sig, marker_times, std_sig = vali01.load_hdf5_data(
        file_path="/home/mateusz-wawrzyniak/PycharmProjects/brcx_lfp_model/source/src1_validation_visualization/val2_invivo_lfp/invivo_data/C57v31_mean_clean.hdf",
        paradigm=conf0.STIM_PARADIGM_TYPE,
        subtype=conf0.STIM_PARADIGM_SUBTYPE_INVIVO
    )
    c57_sig_df = vali01.get_significant_el_sites(
        av_sig=av_sig,
        save_path=C57_DF,
        el_topo_variant = conf02.ELECTRODE_VARIANT,
        dvs_map = conf02.DVS_MAP,
        dt = conf01.DT,
        electrode_offset=conf02.ELECTRODE_Z_OFFSET
    )
    np.savetxt(C57_MARKERS, marker_times, delimiter=",")

    # Plotting the signals
    BTBR_LFP = ut1.new_file_in_dir(paths["visualizations"]["invivo_lfp"]["invivo_exported"], "net_invivo_LFP_btbr.jpg")
    C57_LFP = ut1.new_file_in_dir(paths["visualizations"]["invivo_lfp"]["invivo_exported"], "net_invivo_LFP_c57.jpg")
    vali01.plot_invivo_lfp(btbr_sig_df, save_path=BTBR_LFP)
    vali01.plot_invivo_lfp(c57_sig_df, save_path=C57_LFP)

################
    RECONSTRUCTED_LFP = paths["recordings"]["lfp"]["reconstructed"]["component_lfp.hdf"]
    MODEL_LFP_CSV = ut1.new_file_in_dir(paths["recordings"]["invivo_lfp"]["invivo_exported"], "model_lfp.csv")
    vali02.export_net_lfp_to_csv(
        lfp_hdf_path=RECONSTRUCTED_LFP,
        csv_path=MODEL_LFP_CSV)

    # TRIMMING
    TRIMMED_MODEL = ut1.new_file_in_dir(paths["recordings"]["invivo_lfp"]["trimmed"],"model_lfp_01_trimmed.csv")
    TRIMMED_INVIVO = ut1.new_file_in_dir(paths["recordings"]["invivo_lfp"]["trimmed"], "invivo_lfp_01_trimmed.csv")
    TRIMMED_MARKERS = ut1.new_file_in_dir(paths["recordings"]["invivo_lfp"]["trimmed"],"markers.csv")
    vali02.trim_signals(
        model_lfp_csv=MODEL_LFP_CSV,
        emp_lfp_csv=C57_DF,
        emp_marker_csv=C57_MARKERS,
        paradigm_type=conf0.STIM_PARADIGM_TYPE,
        paradigm_subtype=conf0.STIM_PARADIGM_SUBTYPE,
        model_csv_savepath=TRIMMED_MODEL,
        invivo_csv_savepath=TRIMMED_INVIVO,
        markers_csv_savepath=TRIMMED_MARKERS,
        whisker_stim_paradigms=conf0.WHISKER_STIMULATION_PARADIGMS,
        paradigm_stand=conf02.PARADIGMS_STANDARIZATION,
        dt=conf01.DT
    )
    paradigm_stand_dict = conf02.PARADIGMS_STANDARIZATION[conf0.STIM_PARADIGM_TYPE]

    #BASELINE
    BASELINED_MODEL = ut1.new_file_in_dir(paths["recordings"]["invivo_lfp"]["baselined"],"model_lfp_02_baselined.csv")
    BASELINED_INVIVO = ut1.new_file_in_dir(paths["recordings"]["invivo_lfp"]["baselined"],"invivo_lfp_02_baselined.csv")
    vali02.set_signal_baseline(
        lfp_csv=TRIMMED_INVIVO,
        baseline_window=paradigm_stand_dict['baseline_win'],
        savepath=BASELINED_INVIVO)
    vali02.set_signal_baseline(
        lfp_csv=TRIMMED_MODEL,
        baseline_window=paradigm_stand_dict['baseline_win'],
        savepath=BASELINED_MODEL)

    # FILTERING
    FILT_MODEL = ut1.new_file_in_dir(paths["recordings"]["invivo_lfp"]["filtered"],"model_lfp_03_filtered.csv")
    FILT_INVIVO = ut1.new_file_in_dir(paths["recordings"]["invivo_lfp"]["filtered"],"invivo_lfp_03_filtered.csv")
    vali02.bandpass_filter_signal(
        lfp_csv=BASELINED_INVIVO,
        lp_cutoff=paradigm_stand_dict['lp_freq'],
        hp_cutoff=paradigm_stand_dict['hp_freq'],
        order=4,
        fs = 1000/conf01.DT,
        savepath=FILT_INVIVO
    )
    vali02.bandpass_filter_signal(
        lfp_csv=BASELINED_MODEL,
        lp_cutoff=paradigm_stand_dict['lp_freq'],
        hp_cutoff=paradigm_stand_dict['hp_freq'],
        order=4,
        fs=1000/conf01.DT,
        savepath=FILT_MODEL
    )

    # PLOTS
    chosen_plot_func = vali02.plot_comparison_lfp
    invivo_offset = 40
    model_offset = 0.008
    title_invivo = "LFP registered in-vivo:"
    title_model = "LFP reconstructed:"

    INVIVO_LFP_PLOT_TRIMMED = ut1.new_file_in_dir(paths["visualizations"]["invivo_lfp"]["trimmed"], "emp_lfp.jpg")
    MODEL_LFP_PLOT_TRIMMED = ut1.new_file_in_dir(paths["visualizations"]["invivo_lfp"]["trimmed"], "model_lfp.jpg")
    chosen_plot_func(
        lfp_signals_csv = TRIMMED_INVIVO,
        title = f"{title_invivo} peristimulus trimming applied",
        save_path=INVIVO_LFP_PLOT_TRIMMED,
        offset=invivo_offset
    )
    chosen_plot_func(
        lfp_signals_csv=TRIMMED_MODEL,
        title=f"{title_model} peristimulus trimming applied",
        save_path=MODEL_LFP_PLOT_TRIMMED,
        offset=model_offset
    )

    INVIVO_LFP_PLOT_BASED = ut1.new_file_in_dir(paths["visualizations"]["invivo_lfp"]["baselined"], "emp_lfp.jpg")
    MODEL_LFP_PLOT_BASED = ut1.new_file_in_dir(paths["visualizations"]["invivo_lfp"]["baselined"], "model_lfp.jpg")
    chosen_plot_func(
        lfp_signals_csv = BASELINED_INVIVO,
        title=f"{title_invivo} constant offset corrected",
        save_path=INVIVO_LFP_PLOT_BASED,
        offset=invivo_offset
    )
    chosen_plot_func(
        lfp_signals_csv=BASELINED_MODEL,
        title=f"{title_model} constant offset corrected",
        save_path=MODEL_LFP_PLOT_BASED,
        offset=model_offset
    )

    INVIVO_LFP_PLOT_FILT = ut1.new_file_in_dir(paths["visualizations"]["invivo_lfp"]["filtered"], "emp_lfp.jpg")
    MODEL_LFP_PLOT_FILT = ut1.new_file_in_dir(paths["visualizations"]["invivo_lfp"]["filtered"], "model_lfp.jpg")
    chosen_plot_func(
        lfp_signals_csv = FILT_INVIVO,
        title=f"{title_invivo} band-pass filtered",
        save_path=INVIVO_LFP_PLOT_FILT,
        offset=invivo_offset
    )
    chosen_plot_func(
        lfp_signals_csv=FILT_MODEL,
        title=f"{title_model} band-pass filtered",
        save_path=MODEL_LFP_PLOT_FILT,
        offset=model_offset
    )

    # ADDITIONAL INTERESTING RESPONSE:
    INIT_LFP_RESPONSE = ut1.new_file_in_dir(paths["visualizations"]["invivo_lfp"], "init_lfp.csv")
    INIT_LFP_RESPONSE_PLOT = ut1.new_file_in_dir(paths["visualizations"]["invivo_lfp"], "init_lfp.jpg")
    vali02.manual_lfp_signal_processing(
        lfp_csv=MODEL_LFP_CSV,
        save_path=INIT_LFP_RESPONSE,
        trim_range=(10, 160),
        baseline_window=(10, 50),
        lp_cutoff=paradigm_stand_dict['lp_freq'],
        hp_cutoff=paradigm_stand_dict['hp_freq'],
        order=4,
        fs=1000/conf01.DT
    )

    chosen_plot_func(
        lfp_signals_csv=INIT_LFP_RESPONSE,
        title=f"LFP during the initialization of the cortical cells, at the beginning of the simulation",
        save_path=INIT_LFP_RESPONSE_PLOT,
        offset=model_offset,
        make_dashed=False
    )

