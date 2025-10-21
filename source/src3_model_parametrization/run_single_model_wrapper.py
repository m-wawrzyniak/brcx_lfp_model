import argparse
import importlib.util
import sys

parser = argparse.ArgumentParser()
parser.add_argument("--conf0", required=True)
parser.add_argument("--conf1", required=True)
parser.add_argument("--conf2", required=True)
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

# import your pipeline
from source.src0_core.cr0_model_setup import cr0_main as cr0
from source.src0_core.cr1_simulation_run import cr1_main as cr1
from source.src0_core.cr2_lfp_reconstruction import cr2_main as cr2

cr0.run()
cr1.run()
cr2.run()
