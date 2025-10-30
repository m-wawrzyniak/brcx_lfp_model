import os
import sys
import shutil
import importlib.util
import ast
from datetime import datetime
import subprocess
import astor

from source.src3_model_parametrization import BATCH_PARAMETERS
from config_templates import conf0_model_parameters as default_conf

PROJECT_ROOT = default_conf.ROOT
TEMPLATE_CONFIG_DIR = os.path.join(PROJECT_ROOT, "config_templates")
DATA_ROOT = os.path.join(PROJECT_ROOT, "data")

def create_run_config(base_config_dir: str, model_dir: str, overrides: dict):
    """Copy config templates into per-model folder and apply overrides safely."""
    os.makedirs(model_dir, exist_ok=True)

    for fname in os.listdir(base_config_dir):
        if not fname.endswith(".py"):
            continue
        src = os.path.join(base_config_dir, fname)
        dst = os.path.join(model_dir, fname)
        shutil.copy(src, dst)

        if not overrides:
            continue

        # Safely apply overrides
        with open(dst, "r") as f:
            tree = ast.parse(f.read(), filename=dst)

        class ConfigTransformer(ast.NodeTransformer):
            def visit_Assign(self, node):
                if isinstance(node.targets[0], ast.Name):
                    key = node.targets[0].id
                    if key in overrides:
                        # Replace value safely
                        node.value = ast.parse(repr(overrides[key])).body[0].value
                return node

        tree = ConfigTransformer().visit(tree)
        ast.fix_missing_locations(tree)

        with open(dst, "w") as f:
            f.write(astor.to_source(tree))

def load_conf_module(conf_path: str, module_name: str):
    """
    Dynamically load a config module and inject into sys.modules.
    All downstream imports see this version.
    """
    if module_name in sys.modules:
        del sys.modules[module_name]

    spec = importlib.util.spec_from_file_location(module_name, conf_path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)

    sys.modules[module_name] = mod
    return mod



def run_single_model_subprocess(conf0_p, conf01_p, conf02_p,
                                run_cr0=True, run_cr1=True, run_cr2=True):
    """Run a single simulation in a separate Python process to isolate NEURON."""
    cmd = [
        "python",  # or sys.executable to be safe
        os.path.join(PROJECT_ROOT, "source/src3_model_parametrization/run_single_model_wrapper.py"),
        "--conf0", conf0_p,
        "--conf1", conf01_p,
        "--conf2", conf02_p,
        "--run_cr0", str(int(run_cr0)),
        "--run_cr1", str(int(run_cr1)),
        "--run_cr2", str(int(run_cr2))
    ]
    subprocess.run(cmd, check=True)

def main():
    timestamp = datetime.now().strftime("%m-%d_%H-%M")
    print(f"[BATCH SIMULATION] Running batch simulation at {timestamp}")

    # Define parameter sweep
    PARAM_SETS = BATCH_PARAMETERS.PARAMETER_SET_B

    for params in PARAM_SETS:
        model_dir = os.path.join(DATA_ROOT, params["MODEL_NAME"], "config")
        create_run_config(TEMPLATE_CONFIG_DIR, model_dir, params)

        conf0_p = os.path.join(model_dir, "conf0_model_parameters.py")
        conf01_p = os.path.join(model_dir, "conf01_simulation_parameters.py")
        conf02_p = os.path.join(model_dir, "conf02_lfp_parameters.py")

        run_single_model_subprocess(conf0_p, conf01_p, conf02_p, run_cr0=False, run_cr1=False)

if __name__ == "__main__":
    main()