# Electrode
ELECTRODE_VARIANT = 'center_col_sparse' #'full'
ELECTRODE_Z_OFFSET = -5
CROSS_SPECIES_SCALE = 2

# invivo registration
DVS_MAP = {
    'l' : [47, 46, 45, 44, 43, 42, 41, 40, 39, 38, 37, 36, 35, 34, 33, 48, 49, 50, 51, 52, 53],
    'c' : [1, 64, 2, 63, 3, 62, 4, 61, 5, 60, 6, 59, 7, 58, 8, 57, 9, 56, 10, 55, 11, 54],
    'r' : [17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 16, 15, 14, 13, 12]
}
PARADIGMS_STANDARIZATION = {
    'single': {
        'timerange':(-200, 200),
        'baseline_win':(-150, 0),
        'lp_freq': 1900,
        'hp_freq': 5
    }
}