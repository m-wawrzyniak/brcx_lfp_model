import os
import sys
import shutil
import re
import subprocess
from datetime import datetime

from source.src3_model_parametrization import BATCH_PARAMETERS
from config_templates import conf0_model_parameters as default_conf

PROJECT_ROOT = default_conf.ROOT
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
        with open(dst, "r+", encoding="utf-8") as f:
            text = f.read()
            for key, val in overrides.items():
                text = re.sub(rf"^{key}\s*=.*$", f"{key} = {repr(val)}", text, flags=re.MULTILINE)
            f.seek(0)
            f.write(text)
            f.truncate()


def run_validation_subprocess(conf0_p, conf01_p, conf02_p,
                              run_dep, run_indep):
    """Run val0 + val1 in a separate Python process with specified configs."""
    cmd = [
        sys.executable,
        os.path.join(PROJECT_ROOT, "source/src1_validation_visualization/run_validation_wrapper.py"),
        "--conf0", conf0_p,
        "--conf1", conf01_p,
        "--conf2", conf02_p,
        "--vis_dep", str(int(run_dep)),
        "--vis_indep", str(int(run_indep)),
    ]
    subprocess.run(cmd, check=True)


def main(run_sim_dep=True, run_sim_indep=True):
    timestamp = datetime.now().strftime("%m-%d_%H-%M")
    print(f"[BATCH VALIDATION] Running batch validation at {timestamp}")

    PARAM_SETS = BATCH_PARAMETERS.PARAMETER_SET_B

    for params in PARAM_SETS:
        model_name = params["MODEL_NAME"]
        model_dir = os.path.join(DATA_ROOT, model_name, "config")

        # Copy and override config templates
        create_run_config(TEMPLATE_CONFIG_DIR, model_dir, params)
        conf0_p = os.path.join(model_dir, "conf0_model_parameters.py")
        conf01_p = os.path.join(model_dir, "conf01_simulation_parameters.py")
        conf02_p = os.path.join(model_dir, "conf02_lfp_parameters.py")

        run_validation_subprocess(conf0_p, conf01_p, conf02_p, run_sim_dep, run_sim_indep)


if __name__ == "__main__":
    main(run_sim_indep=False)
