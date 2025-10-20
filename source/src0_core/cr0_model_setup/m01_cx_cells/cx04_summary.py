import json
import pandas as pd
import os
from glob import glob

def get_population_me_composition(pop_files_dir: str, save_path:str):
    """
    Based on all *.csv files in specific pop_files_dir, creates a summary me_summary.json file, where each
    key is a <layer name> and within the layer are keys:
        - 'tot_n_cells'
        - 'inh_fraction'
        - 'me_types': {
            - <me_type>: <me_type_fraction_in_population>
            }

    Args:
        pop_files_dir (str): Directory with all *csv concerning somata positions.
    """
    summary = {}

    csv_files = glob(os.path.join(pop_files_dir, "*.csv"))
    if not csv_files:
        raise FileNotFoundError(f"No CSV files found in {pop_files_dir}")

    dfs = [pd.read_csv(f) for f in csv_files]
    all_cells = pd.concat(dfs, ignore_index=True)

    # Per-layer
    for layer, layer_df in all_cells.groupby("layer"):
        tot_n_cells = len(layer_df)
        me_type_counts = layer_df['desc'].value_counts(normalize=True).round(4).to_dict()

        summary[layer] = {
            'tot_n_cells': tot_n_cells,
            'me_types': me_type_counts
        }

    # Whole population
    total_n = len(all_cells)
    global_me_types = all_cells['desc'].value_counts(normalize=True).round(4).to_dict()
    summary['population'] = {
        'tot_n_cells': total_n,
        'me_types': global_me_types
    }

    # Save JSON
    with open(save_path, "w") as f:
        json.dump(summary, f, indent=4)
