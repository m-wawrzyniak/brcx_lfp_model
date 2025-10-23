import os
import sys
import shutil
import importlib.util
import re
from datetime import datetime
import subprocess

from source.src3_model_parametrization.BATCH_PARAMETERS import PARAMETER_SET_A

PROJECT_ROOT = "/home/mateusz-wawrzyniak/PycharmProjects/brcx_lfp_model"
TEMPLATE_CONFIG_DIR = os.path.join(PROJECT_ROOT, "config_templates")
DATA_ROOT = os.path.join(PROJECT_ROOT, "data")

def create_run_config(base_config_dir: str, model_dir: str, overrides: dict):
    """Copy config templates into per-model folder and apply overrides."""
    os.makedirs(model_dir, exist_ok=True)
    for fname in os.listdir(base_config_dir):
        if not fname.endswith(".py"):
            continue
        src = os.path.join(base_config_dir, fname)
        dst = os.path.join(model_dir, fname)
        shutil.copy(src, dst)
        # Apply overrides (simple key=value replacement)
        with open(dst, "r+") as f:
            text = f.read()
            for key, val in overrides.items():
                text = re.sub(rf"^{key}\s*=.*$", f"{key} = {repr(val)}", text, flags=re.MULTILINE)
            f.seek(0)
            f.write(text)
            f.truncate()

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



def run_single_model_subprocess(conf0_p, conf01_p, conf02_p):
    """Run a single simulation in a separate Python process to isolate NEURON."""
    cmd = [
        "python",  # or sys.executable to be safe
        os.path.join(PROJECT_ROOT, "source/src3_model_parametrization/run_single_model_wrapper.py"),
        "--conf0", conf0_p,
        "--conf1", conf01_p,
        "--conf2", conf02_p
    ]
    subprocess.run(cmd, check=True)

def main():
    timestamp = datetime.now().strftime("%m-%d_%H-%M")
    print(f"[BATCH SIMULATION] Running batch simulation at {timestamp}")

    # Define parameter sweep
    PARAM_SETS = PARAMETER_SET_A

    for params in PARAM_SETS:
        model_dir = os.path.join(DATA_ROOT, params["MODEL_NAME"], "config")
        create_run_config(TEMPLATE_CONFIG_DIR, model_dir, params)

        conf0_p = os.path.join(model_dir, "conf0_model_parameters.py")
        conf01_p = os.path.join(model_dir, "conf01_simulation_parameters.py")
        conf02_p = os.path.join(model_dir, "conf02_lfp_parameters.py")

        run_single_model_subprocess(conf0_p, conf01_p, conf02_p)

if __name__ == "__main__":
    main()