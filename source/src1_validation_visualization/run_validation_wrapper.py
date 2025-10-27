import argparse
import importlib.util
import sys

parser = argparse.ArgumentParser()
parser.add_argument("--conf0", required=True)
parser.add_argument("--conf1", required=True)
parser.add_argument("--conf2", required=True)
parser.add_argument("--vis_dep", required=True)
parser.add_argument("--vis_indep", required=True)
args = parser.parse_args()


def load_conf_module(conf_path: str, module_name: str):
    spec = importlib.util.spec_from_file_location(module_name, conf_path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    sys.modules[module_name] = mod
    return mod


conf0 = load_conf_module(args.conf0, "config_templates.conf0_model_parameters")
conf01 = load_conf_module(args.conf1, "config_templates.conf01_simulation_parameters")
conf02 = load_conf_module(args.conf2, "config_templates.conf02_lfp_parameters")

from source.src1_validation_visualization.val0_sim_indep import val0_main as val0
from source.src1_validation_visualization.val1_sim_dep import val1_main as val1
from source.src1_validation_visualization.val2_invivo_lfp import val2_main as val2

def run():
    run_dep = bool(int(args.vis_dep))
    run_indep = bool(int(args.vis_indep))

    print(f"[WRAPPER] vis_indep={run_indep}, vis_dep={run_dep}")

    if run_indep:
        val0.run(model_name=conf0.MODEL_NAME)

    if run_dep:
        val1.run(model_name=conf0.MODEL_NAME)

    val2.run(model_name=conf0.MODEL_NAME)


if __name__ == "__main__":
    run()
