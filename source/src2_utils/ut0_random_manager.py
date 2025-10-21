import random
import numpy as np
from numpy.random import default_rng

from config_templates.conf0_model_parameters import GLOBAL_SEED

rng = default_rng(GLOBAL_SEED)

def set_global_seed(seed=GLOBAL_SEED):
    global rng
    random.seed(seed)
    np.random.seed(seed)
    rng = default_rng(seed)
    return seed


set_global_seed()