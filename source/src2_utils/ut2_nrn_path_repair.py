import re
from pathlib import Path

import config_templates.conf0_model_parameters as conf0

def fix_hoc_paths(project_root: str, subdir: str = "config_templates/cx_cell_templates"):
    """
    Fixes absolute paths inside all .hoc files in the cell template directories.

    Replaces:
      load_file("...oldpath...")
      nl.input("...oldpath...")

    With absolute paths relative to each cell type folder.
    """
    root = Path(project_root).resolve()
    base_dir = root / subdir

    # Two regex patterns
    load_file_pattern = re.compile(r'load_file\("([^"]+)"\)')
    nl_input_pattern = re.compile(r'nl\.input\("([^"]+)"\)')

    for cell_dir in base_dir.iterdir():
        if not cell_dir.is_dir():
            continue

        hoc_files = list(cell_dir.glob("*.hoc"))
        if not hoc_files:
            print(f"No .hoc files in {cell_dir}")
            continue

        for hoc_file in hoc_files:
            with open(hoc_file, "r") as f:
                content = f.read()

            def replace_load_file(match):
                old_path = match.group(1)
                filename = Path(old_path).name
                new_path = str((cell_dir / filename).resolve())
                return f'load_file("{new_path}")'

            def replace_nl_input(match):
                old_path = match.group(1)
                filename = Path(old_path).name
                new_path = str((cell_dir / "morphology" / filename).resolve())
                return f'nl.input("{new_path}")'

            new_content = load_file_pattern.sub(replace_load_file, content)
            new_content = nl_input_pattern.sub(replace_nl_input, new_content)

            if new_content != content:
                with open(hoc_file, "w") as f:
                    f.write(new_content)
                print(f"Fixed paths in: {hoc_file}")
            else:
                print(f" No changes needed for: {hoc_file}")

if __name__ == "__main__":
    PROJECT_ROOT = conf0.ROOT
    fix_hoc_paths(PROJECT_ROOT)
