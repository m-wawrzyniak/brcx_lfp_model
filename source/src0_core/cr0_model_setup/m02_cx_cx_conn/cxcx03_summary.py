import json
from collections import defaultdict

def create_synapse_summary(synapses_json:str, summary_path:str):
    """
    Summarizes synapses.json as:
        - each key is "<pre_me_type>:<post_me_type>"
        - each value is another dict : {
            "n_synapses": number of synapses for this combination
            "n_connections": number of connections for this combination (connection is unique cell pair forming at least one synapse)
            }
    Args:
        synapses_json (str): Final synapses after pruning in json.
        summary_dir (str): Directoy in which the summary should be saved.

    """
    # Load JSON if path is given
    if isinstance(synapses_json, str):
        with open(synapses_json, "r") as f:
            syn_data = json.load(f)
    else:
        syn_data = synapses_json  # assume it's already a dict

    summary = defaultdict(lambda: {"n_synapses": 0, "n_connections": set()})

    for post_id_str, syn_list in syn_data.items():
        post_id = int(post_id_str)
        for syn in syn_list:
            pre_id = syn["pre_id"]
            pre_me = syn["pre_me_type"]
            post_me = syn["post_me_type"]
            pre_short = "_".join(pre_me.split('_')[:-1])
            post_short = "_".join(post_me.split('_')[:-1])
            key = f"{pre_short}:{post_short}"

            summary[key]["n_synapses"] += 1
            summary[key]["n_connections"].add((pre_id, post_id))

    # Convert sets to counts
    final_summary = {
        key: {
            "n_synapses": val["n_synapses"],
            "n_connections": len(val["n_connections"])
        }
        for key, val in summary.items()
    }

    # Save to file
    with open(summary_path, "w") as f:
        json.dump(final_summary, f, indent=4)