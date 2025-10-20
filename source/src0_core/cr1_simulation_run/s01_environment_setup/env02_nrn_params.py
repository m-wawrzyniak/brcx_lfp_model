import config_templates.conf01_simulation_parameters as conf1


def set_simulation_params(h_envi):
    h_envi.celsius = conf1.TEMP
    return h_envi