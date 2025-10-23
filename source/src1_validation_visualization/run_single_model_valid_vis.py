
from source.src1_validation_visualization.val0_sim_indep import val0_main as val0
from source.src1_validation_visualization.val1_sim_dep import val1_main as val1

MODEL_NAME = ""

def run():
    val0.run(model_name=MODEL_NAME)
    val1.run(model_name=MODEL_NAME)