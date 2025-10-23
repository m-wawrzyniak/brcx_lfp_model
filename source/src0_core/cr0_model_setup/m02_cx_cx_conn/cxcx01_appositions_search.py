from scipy.spatial import cKDTree

def find_pop_appositions(max_dist, all_cells: dict) -> tuple[dict, int]:
    print(f"[cxcx01] Finding cxcx appositions with max_dist={max_dist}")

    syn_cnt = 0
    all_synapses = {}

    for post_id, post_cell in all_cells.items():
        syns = []
        post_pts = post_cell['dendrite_pts']
        post_pts_with_sec = post_cell['topo'].dendrite_pts_with_sec
        post_tree = cKDTree(post_pts)
        post_soma = post_cell['pos']

        for pre_id, pre_cell in all_cells.items():
            if pre_id == post_id:  # skip self
                continue

            contacts = post_tree.query_ball_point(pre_cell['axon_pts'], r=max_dist)
            for pre_idx, post_idxs in enumerate(contacts):
                for post_idx in post_idxs:
                    contact_global = post_pts[post_idx]
                    _, _, _, sec_name = post_pts_with_sec[post_idx]

                    contact_relative = contact_global - post_soma
                    syns.append({
                        "pre_id": pre_id,
                        "pre_me_type": f"{pre_cell['layer']}_{pre_cell['m_type']}_{pre_cell['e_type']}",
                        "post_me_type": f"{post_cell['layer']}_{post_cell['m_type']}_{post_cell['e_type']}",
                        "post_loc": contact_relative.tolist(),
                        "post_sec": sec_name
                    })
        syn_cnt += len(syns)
        all_synapses[post_id] = syns

    print(f"[cxcx01] SUCCESS: cxcx appositions found with max_dist={max_dist}")

    return all_synapses, syn_cnt