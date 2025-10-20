from datetime import datetime
from pathlib import Path
import json

from config_templates.conf0_model_parameters import MODEL_NAME, ROOT

def create_structure(structure, root=None, **placeholders):
    if root is None:
        root = ROOT

    timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    placeholders = {"timestamp": timestamp, "model_name": MODEL_NAME, **placeholders}

    for name, subtree in structure.items():
        name = name.format(**placeholders)
        path = Path(root) / name.rstrip("/")
        if name.endswith("/"):
            path.mkdir(parents=True, exist_ok=True)
            if subtree:
                create_structure(subtree, path, **placeholders)
        else:
            path.touch(exist_ok=True)

def load_structure(json_path):
    with open(json_path, 'r') as f:
        return json.load(f)

def load_tree(config_file, root=None, **vars):
    if root is None:
        root = ROOT

    with open(config_file) as f:
        tree = json.load(f)

    def resolve(node, prefix=Path(root)):
        out = {}
        for raw_key, val in node.items():
            key = raw_key.format(**vars).strip("/")
            path = prefix / key
            if isinstance(val, dict):
                out[key] = resolve(val, path)
            else:
                out[key] = path
        return out

    return resolve(tree)


def __main__():
    model_structure_path = "/home/mateusz-wawrzyniak/PycharmProjects/brcx_lfp_model/config_templates/model_structure_template.json"
    structure = load_structure(model_structure_path)
    create_structure(structure)

    return load_tree(config_file=model_structure_path, model_name=MODEL_NAME)
