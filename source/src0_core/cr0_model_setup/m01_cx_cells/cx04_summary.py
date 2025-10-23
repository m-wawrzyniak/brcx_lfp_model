import json
import pandas as pd
import os
from glob import glob

def get_population_me_composition(pop_files_dir: str, save_path:str):
    print("[cx04] Creating me-type composition summary.")

    summary = {}
    csv_files = glob(os.path.join(pop_files_dir, "*.csv"))
    if not csv_files:
        raise FileNotFoundError(f"No CSV files found in {pop_files_dir}")

    dfs = [pd.read_csv(f) for f in csv_files]
    all_cells = pd.concat(dfs, ignore_index=True)

    for layer, layer_df in all_cells.groupby("layer"):
        tot_n_cells = len(layer_df)
        me_type_counts = layer_df['desc'].value_counts(normalize=True).round(4).to_dict()

        summary[layer] = {
            'tot_n_cells': tot_n_cells,
            'me_types': me_type_counts
        }

    total_n = len(all_cells)
    global_me_types = all_cells['desc'].value_counts(normalize=True).round(4).to_dict()
    summary['population'] = {
        'tot_n_cells': total_n,
        'me_types': global_me_types
    }

    with open(save_path, "w") as f:
        json.dump(summary, f, indent=4)

    print(f"[cx04] SUCCESS: Creating me-type composition summary at {save_path}.")