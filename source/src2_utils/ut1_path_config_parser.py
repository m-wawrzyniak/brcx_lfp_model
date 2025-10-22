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


def folder_to_json_with_files(root_path: str, json_file: str):
    """
    Recursively scans a folder and saves a JSON structure compatible with `create_structure`,
    including both directories and files.

    Directories -> nested dictionaries (or null if empty)
    Files -> their relative paths as strings
    """
    root = Path(root_path)
    if not root.exists() or not root.is_dir():
        raise ValueError(f"{root} is not a valid directory")

    def scan(path: Path, prefix=""):
        tree = {}
        for child in sorted(path.iterdir()):
            if child.is_dir():
                key = f"{child.name}/"
                subtree = scan(child, prefix=prefix + child.name + "/")
                # if directory is empty, store null instead of {}
                tree[key] = subtree if subtree else None
            else:
                # store file as relative path from root
                key = child.name
                tree[key] = str(Path(prefix) / child.name)
        return tree

    structure = scan(root)

    with open(json_file, "w") as f:
        json.dump(structure, f, indent=2)

    print(f"JSON structure with files saved to {json_file}")

# TODO: change all os.path.join with this!
def new_file_in_dir(dir_dict: dict, filename: str) -> Path:
    """Return a new path in the same folder as any of the given paths in dict."""
    example_path = next(iter(dir_dict.values()))
    return example_path.parent / filename


def __main__():
    model_structure_path = "/home/mateusz-wawrzyniak/PycharmProjects/brcx_lfp_model/config_templates/model_structure_template.json"
    structure = load_structure(model_structure_path)
    create_structure(structure)

    return load_tree(config_file=model_structure_path, model_name=MODEL_NAME)
