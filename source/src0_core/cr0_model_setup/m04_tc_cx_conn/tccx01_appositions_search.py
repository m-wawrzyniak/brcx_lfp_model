import os.path
import json
from collections import defaultdict
import matplotlib.pyplot as plt

import config_templates.conf0_model_parameters as conf0

from source.src2_utils.ut0_random_manager import np, random

def get_segments_binning(all_cells: dict, h_start: float, h_stop: float) -> dict:
    print("[tccx01] Segmenting cx layers with respect to z-axis.")

    z_bin_size = conf0.Z_BIN_SIZE

    if h_start > h_stop:
        n_bins = int((h_start - h_stop) / z_bin_size)
        bin_edges = [h_start - i * z_bin_size for i in range(n_bins)]
    else:
        n_bins = int((h_stop - h_start) / z_bin_size)
        bin_edges = [h_start + i * z_bin_size for i in range(n_bins)]

    bin_dict = {bin_z: [] for bin_z in bin_edges}

    for c_id, c_info in all_cells.items():
        cell_tag = f"{c_info['layer']}_{c_info['m_type']}"
        if cell_tag not in conf0.TC_TARGET_MTYPES:
            continue

        cell_center = c_info['pos']
        dend_segs = c_info['topo'].dendrite_pts_with_sec  # (x,y, z, sec_name)

        for x, y, z, sec_name in dend_segs:
            rel_seg = np.array([x, y, z]) - cell_center
            if (h_start >= z > h_stop) or (h_start <= z < h_stop):
                bin_idx = int((z - h_start) // z_bin_size) if h_stop > h_start else int((h_start - z) // z_bin_size)
                bin_z = h_start + bin_idx * z_bin_size if h_stop > h_start else h_start - bin_idx * z_bin_size
                if bin_z in bin_dict:
                    bin_dict[bin_z].append((c_id, cell_tag, sec_name, tuple(rel_seg)))
    print("[tccx01] SUCCESS: Segmenting cx layers with respect to z-axis.")

    return bin_dict

def calc_tc_bd_dist(h_start:float, h_stop:float, save_path) -> np.ndarray:
    print("[tccx01] Approximating empirical tc bouton density on cx cells with respect to z-axis")

    bin_size = conf0.Z_BIN_SIZE
    if h_start > h_stop:
        bin_centers = np.arange(h_start, h_stop, -bin_size)
    else:
        bin_centers = np.arange(h_start, h_stop, bin_size)

    bd_values = np.zeros_like(bin_centers, dtype=float)

    for cluster in conf0.TC_BD_DIST.values():
        A = cluster['mean_bd']
        mu = cluster['center']
        sigma = cluster['std_bd']
        bd_values += A * np.exp(-0.5 * ((bin_centers - mu) / sigma)**2)

    result = np.column_stack((bin_centers, bd_values))


    plt.figure(figsize=(6, 4))
    plt.plot(result[:, 0], result[:, 1], drawstyle='steps-mid')
    plt.xlabel("Depth from pia [um]")
    plt.ylabel("Bouton Density [$10^7 / mm^3$]")
    plt.title("Empirical TC-CX bouton density distribution with respect to cortical depth")
    plt.grid(True)
    plt.tight_layout()
    save_path = os.path.join(save_path)
    plt.savefig(save_path, dpi=200)

    print(f"[tccx01] SUCCESS: Approximating empirical tc bouton density saved at {save_path}")

    return result

def tc_synapse_sampling(z_bin_dict:dict, bd_emp:np.ndarray, save_path: str):
    print("[tccx01] Sampling cx segments for tccx based on emp.dist.")

    syn_dict = defaultdict(list)

    for bin_z, bd in bd_emp:
        segments = z_bin_dict.get(bin_z, [])
        if not segments or bd <= 0:
            continue

        bin_volume = np.pi * (conf0.RADIUS ** 2) * conf0.Z_BIN_SIZE / 1e9  # um3 >> mm3
        n_synapses = max(1, int(bd * 1e7 * bin_volume * conf0.TCCX_SYNAPSE_SCALE))
        sampled = random.choices(segments, k=n_synapses)

        for post_id, post_me_type, sec_name, rel_pos in sampled:
            syn_dict[str(post_id)].append({
                "pre_id": None,
                "pre_me_type": None,
                "post_me_type": post_me_type,
                "post_sec": sec_name,
                "post_loc": list(rel_pos)
            })

    with open(save_path, "w") as f:
        json.dump(syn_dict, f, indent=2)
    tot_syn = sum(len(v) for v in syn_dict.values())

    print(f"[tccx01] SUCCESS: Sampling cx segments for tccx based on emp.dist. Segments saved at {str(save_path)}")
