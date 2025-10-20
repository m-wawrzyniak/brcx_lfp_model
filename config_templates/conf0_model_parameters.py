### Global params
GLOBAL_SEED = 121       # Randomization seed
GLOBAL_CX_CELL_CNT = 0  # Count of predefined Cx cells.
GLOBAL_PRV_CELL_CNT = 0
MODEL_NAME = 'test'
ROOT = '/home/mateusz-wawrzyniak/PycharmProjects/brcx_lfp_model'

### m01_cx_cells params

# Target layers and me-types included

# TODO: This can be stored in one single metadata structure
LAYER_COMP_TARGET = {
    "L23": ["cAC", "cNAC", "cADpyr"],
    "L4": ["cAC", "cNAC", "cADpyr"],
    "L5": ["cAC", "cNAC", "cADpyr"]
}

METYPE_Z_RANGES={
    'L23_LBC_cAC': (250, -600),
    'L23_LBC_cNAC': (200, -600),
    'L23_PC_cAD' : (250, -1100), #(500, -1300),
    'L4_LBC_cAC' : (500, -750),
    'L4_LBC_cNAC': (600, -950),
    'L4_SS_cAD': (600, -750),
    'L5_LBC_cNAC': (500, -400),
    'L5_MC_cAC': (1100, -200), #(1300, -400),
    'L5_TTPC1_cAD': (850, -800), #(1050, -1000),
    'L5_TTPC2_cAD': (1100, -500), #(1300, -750)
}

LAYER_COMP_PARAMS = {
    "L23":{
        "LBC_cAC": 0.072,
        "LBC_cNAC": 0.088,
        "PC_cAD": 0.84
    },
    "L4":{
        "LBC_cAC":0.063,
        "LBC_cNAC":0.037,
        "SS_cAD": 0.90
    },
    "L5":{
        "LBC_cNAC":0.050,
        "MC_cAC":0.12,
        "TTPC1_cAD": 0.415,
        "TTPC2_cAD": 0.415
    }
}


GOAL_N_CELLS = 30#200
SCALE_FACTOR = 0.1#0.5           # Scaling factor of the BrCx size. It affects only the number of cells predefined.
GLOBAL_M = GOAL_N_CELLS//3 #50                  # Number of minicolumns at each layer
MIN_SOMA_DISTANCE = 10              # Minimal um distance accepted during somata positioning.
Z_STD = 10/SCALE_FACTOR             # STD applied to z-coord during somata placement.
X_Y_JITTER_STD = 20
GLOBAL_Z_RANGE = (-2000.0, 0)

RADIUS = 120
TISSUE_PARAMS = {
    'L23':{
        'start_H':-165,
        'R': RADIUS ,
        'H': 502 ,
        'D': 10.07,
        'M': GLOBAL_M
    },
    'L4': {
        'start_H':-667,
        'R': RADIUS ,
        'H': 190 ,
        'D': 16.86,
        'M': GLOBAL_M
    },
    'L5': {
        'start_H': -857,
        'R': RADIUS ,
        'H': 525 ,
        'D': 7.62,
        'M': GLOBAL_M
    }
}

# TC

Z_BIN_SIZE = 25.0
TC_TARGET_MTYPES = {'L4_SS', 'L5_TTPC1', 'L5_TTPC2'}
# TODO Calculate so to get proper TC input per CX cell
SYN_PER_TC_CELL = 45  #TODO: Decide on SYN_PER_TC_CELL.
TCCX_SYNAPSE_SCALE = 0.002 #TODO: this has to be scale somehow

TC_BD_DIST = { # mean_bd [10^7 / mm^3], center [um from pia]
    'L3_4 cluster': {
        'mean_bd': 4,
        'std_bd': 400/2.355, # Calculated ad hoc
        'center': -600,
        'FWHM': 400,
    },
    'L5 cluster': {
        'mean_bd': 2.75,
        'std_bd': 125/2.355,
        'center': -1200,
        'FWHM': 125
    }
}

# CXCX params
NORM_BIAS_01 = 0.5 # pruning 1
LAMBDA_02, ALPHA_02 = 1, 1 # pruning 2
# pruning 3
BD_EMPIRICAL_DICT = {
    "L1_DAC": 0.18,
    "L1_DLAC": 0.21,
    "L1_HAC": 0.21,
    "L1_NGC_DA": 0.22,
    "L1_NGC_SA": 0.18,
    "L1_SLAC": 0.19,
    "L23_BP": 0.21,
    "L23_BTC": 0.20,
    "L23_DBC": 0.27,
    "L23_LBC": 0.24,
    "L23_MC": 0.21,
    "L23_NBC": 0.24,
    "L23_SBC": 0.20,
    "L4_BP": 0.23,
    "L4_BTC": 0.20,
    "L4_DBC": 0.19,
    "L4_LBC": 0.21,
    "L4_NBC": 0.20,
    "L4_SBC": 0.19,
    "L4_SP": 0.22,
    "L4_SS": 0.22,
    "L5_BP": 0.20,
    "L5_BTC": 0.21,
    "L5_DBC": 0.21,
    "L5_LBC": 0.17,
    "L5_MC": 0.19,
    "L5_NBC": 0.21,
    "L5_SBC": 0.20,
    "L5_TTPC1": 0.21,
    "L5_TTPC2": 0.15,
    "L5_UTPC": 0.19,
    "L6_BPC": 0.15,
    "L6_DBC": 0.20,
    "L6_LBC": 0.14,
    "L6_MC": 0.18,
    "L6_NBC": 0.19,
    "L6_NGC": 0.19,
    "L6_SBC": 0.19,
    "L6_UTPC": 0.20
}
REP_BD_EMP = 0.2

# PRV params

PRV_TC_MEAN_DELAY = 0
PRV_TC_MEAN_TAU = 1

PRV_FR_STD_FACTOR = 0.2
PRV_FR_FACTOR = 1
PRV_WEAK_FACTOR = 0.7

# PRVTC params
PRV_PER_VPM_CELL = 1

# STIMULATION params
STIM_PARADIGM_TYPE = 'single'
STIM_PARADIGM_SUBTYPE = 'str200'

WHISKER_STIMULATION_PARADIGMS = {
    'rest':{
        'r0':{'rest':('r', 500)}
    },
    'single':{
        'wk200':{'prestim':('r', 200), 'single':('wk', 10), 'poststim':('r', 500-210)},
        'str200':{'prestim':('r', 260), 'single':('str', 10), 'poststim':('r', 770-210)}
    },
    'ppi':{
        's76':{'prestim':('r', 200), 'prepulse':('str', 10), 'interval':('r', 76-10), 'pulse':('wk', 10), 'poststim':('r', 500-296)},
        's109':{'prestim':(), 'prepulse':(), 'interval':(), 'pulse':(), 'poststim':()},
        's209':{}
    },
    'fdd':{
        'f0.5':{},
        'f1':{},
        'f2':{},
        'f4.8':{},
        'f9.2':{},
        'f13.2':{},
        'f16.9':{},
        'f23.6':{},
        'f29.4':{}
    }
}