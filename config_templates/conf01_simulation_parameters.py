# Globals
GLOBAL_SYN_CNT = 0

# Simulation params
TEMP = 21
DT = 0.25
SIM_TIME = 770

RESTING_POTENTIALS = {
    'L23_LBC_cAC': -72.6,
    'L23_LBC_cNAC': -72.6,
    'L23_PC_cAD': -74.5,
    'L4_LBC_cAC': -72,
    'L4_LBC_cNAC': -72,
    'L4_SS_cAD': -76,
    'L5_LBC_cNAC': -65,
    'L5_MC_cAC': -72.6,
    'L5_TTPC1_cAD': -73.6,
    'L5_TTPC2_cAD': -72.6
}

# Weights etc
WEIGHT_MATRIX = [
    #>PRV   >VPM    >L23_i  >L23_i  >L23_e  >L4_i   >L4_i   >L4_e   >L5_i   >L5_i   >L5_e   >L5_e
    [0.0000, 0.1000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000],  # PRV ->
    [0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.003, 0.0000, 0.0000, 0.003, 0.003],  # VPM ->
    [0.0000, 0.0000, 0.0005, 0.0005, 0.0050, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000],  # L23_LBC_cAC (inh) ->
    [0.0000, 0.0000, 0.0005, 0.0005, 0.0050, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000],  # L23_LBC_cNAC (inh) ->
    [0.0000, 0.0000, 0.0050, 0.0050, 0.0000, 0.0000, 0.0000, 0.0000, 0.0050, 0.0050, 0.0009, 0.0009],  # L23_PC_cAD (pyr) ->
    [0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0005, 0.0005, 0.0050, 0.0000, 0.0000, 0.0000, 0.0000],  # L4_LBC_cAC (inh) ->
    [0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0005, 0.0005, 0.0050, 0.0000, 0.0000, 0.0000, 0.0000],  # L4_LBC_cNAC (inh) ->
    [0.0000, 0.0000, 0.0050, 0.0050, 0.0009, 0.0050, 0.0050, 0.0009, 0.0000, 0.0000, 0.0000, 0.0000],  # L4_SS_cAD (pyr) - >
    [0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0005, 0.0005, 0.0050, 0.0050],  # L5_LBC_cNAC (inh) ->
    [0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0005, 0.0005, 0.0050, 0.0050],  # L5_MC_cAC (inh) ->
    [0.0000, 0.0000, 0.0050, 0.0050, 0.0009, 0.0050, 0.0050, 0.0009, 0.0050, 0.0050, 0.0009, 0.0009],  # L5_TTPC1_cAD (pyr) ->
    [0.0000, 0.0000, 0.0050, 0.0050, 0.0009, 0.0050, 0.0050, 0.0009, 0.0050, 0.0050, 0.0009, 0.0009],  # L5_TTPC2_cAD (pyr) ->
]
CELL_ORDER = ['PRV', 'VPM', 'L23_LBC_cAC', 'L23_LBC_cNAC', 'L23_PC_cAD', 'L4_LBC_cAC', 'L4_LBC_cNAC', 'L4_SS_cAD', 'L5_LBC_cNAC', 'L5_MC_cAC', 'L5_TTPC1_cAD', 'L5_TTPC2_cAD']
CELL_ORDER_MAP = {name: i for i, name in enumerate(CELL_ORDER)}

# TCCX synapses
TCCX_SYNAPSE_PARAMS = {
    "gsyn_mean": 0.3,
    "epsp_mean": 2.6,
    "risetime_std": 0.63,
    "f_std": 2.7,
    "gsyn_std": 0.11,
    "u_std": 0.048,
    "decay_mean": 61,
    "latency_mean": 1.3,
    "failures_mean": 0,
    "u_mean": 0.86,
    "d_std": 9.1,
    "synapse_type": "Excitatory, depressing",
    "space_clamp_correction_factor": 2.5,
    "latency_std": 0.32,
    "decay_std": 1.9,
    "cv_psp_amplitude_std": 0.036,
    "risetime_mean": 1.9,
    "cv_psp_amplitude_mean": 0.15,
    "epsp_std": 1.1,
    "d_mean": 670,
    "f_mean": 17,
    "failures_std": 0
}
TCCX_NMDA_COMPONENT = 0.1
TCCX_NC_THRESH = -20