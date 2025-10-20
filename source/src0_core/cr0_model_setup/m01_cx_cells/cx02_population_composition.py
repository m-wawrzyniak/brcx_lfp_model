import pandas as pd
import os

from source.src3_utils.ut0_random_manager import np
import config_templates.conf0_model_parameters as conf0


def layer_soma_band_from_tissue(tissue_params: dict, layer: str):
    """(bottom, top) = (start_H - H, start_H)"""
    p = tissue_params[layer]
    top = float(p['start_H'])
    bottom = top - float(p['H'])
    return (bottom, top)

def feasible_soma_interval_for_mtype(layer: str,
                                     me_type_full: str,
                                     tissue_params: dict = conf0.TISSUE_PARAMS,
                                     global_z_range: tuple[float, float] = (-2000.0, 0.0)):
    """
    Uses your sign convention in METYPE_Z_RANGES: (above>0, below<0).
    Returns (lo, hi) feasible soma z interval, or None if empty/unknown.
    """
    if me_type_full not in conf0.METYPE_Z_RANGES:
        return None

    z_above, z_below = conf0.METYPE_Z_RANGES[me_type_full]    # e.g., (250, -600)
    zmin_global, zmax_global = global_z_range           # (-2000, 0)

    # Morph + global:
    # soma <= zmax - above   (no part above 0)
    # soma >= zmin - below   (no part below -2000; below is negative)
    lo = zmin_global - z_below      # -2000 - (-600) = -1400
    hi = zmax_global - z_above      # 0 - 250 = -250

    # Layer band
    lay_lo, lay_hi = layer_soma_band_from_tissue(tissue_params, layer)
    lo = max(lo, lay_lo)
    hi = min(hi, lay_hi)

    return (lo, hi) if lo <= hi else None

def mtypes_feasible_at_soma(layer: str, soma_z: float, candidates_full: list[str]):
    feas = []
    for mfull in candidates_full:
        iv = feasible_soma_interval_for_mtype(layer, mfull)
        if iv is None:
            continue
        lo, hi = iv
        if lo <= soma_z <= hi:
            feas.append(mfull)
    return feas

def nearest_feasible_mtype(layer: str, soma_z: float, candidates_full: list[str]):
    best, best_d = None, float("inf")
    for mfull in candidates_full:
        iv = feasible_soma_interval_for_mtype(layer, mfull)
        if iv is None:
            continue
        lo, hi = iv
        d = 0.0 if lo <= soma_z <= hi else min(abs(soma_z - lo), abs(soma_z - hi))
        if d < best_d:
            best, best_d = mfull, d
    return best, best_d


def assign_me_type(input_dir: str, output_dir: str,
                   layer_comp_params: dict = conf0.LAYER_COMP_PARAMS,
                   verbose: bool = True):
    """
    For each layer CSV in 01_soma_positions:
      - sample full m-type distribution defined in LAYER_COMP_PARAMS[layer]
      - enforce spatial feasibility (global limits, layer band, morph extents)
      - save to 02_me_types_assigned/<Layer>_me_pop.csv
    Assumes CSV has a 'z' column (soma depth, µm; negative down).
    """

    for filename in sorted(os.listdir(input_dir)):
        if not filename.endswith(".csv"):
            continue

        layer_name = filename.split(".")[0].split("_")[0]  # 'L23' from 'L23_pos.csv'
        if layer_name not in layer_comp_params:
            if verbose:
                print(f"[skip] Unrecognized layer in {filename}: {layer_name}")
            continue

        df = pd.read_csv(os.path.join(input_dir, filename))
        if 'z' not in df.columns:
            raise ValueError(f"{filename} must contain a 'z' column.")

        n = len(df)
        # SHORT keys and probs for this layer (e.g., 'LBC_cAC', ...)
        shorts = list(layer_comp_params[layer_name].keys())
        probs = np.array([layer_comp_params[layer_name][s] for s in shorts], dtype=float)
        probs = probs / probs.sum()

        # Map to FULL keys for feasibility (prefix with layer)
        short_to_full = {s: f"{layer_name}_{s}" for s in shorts}

        m_type_list, e_type_list, chosen_full = [], [], []
        rot_angles = np.random.uniform(0, 2*np.pi, size=n)

        # Precompute a preference order per row: we first try the sampled short type,
        # then the remaining types in descending probability to keep global proportions.
        prob_order = list(np.argsort(-probs))  # indices sorted by prob desc

        for i in range(n):
            z_soma = float(df.loc[i, "z"])

            # Draw one short type by target fractions
            s_pick = np.random.choice(shorts, p=probs)

            # Build ordered candidate list: [sampled] + [others by prob desc]
            # Keep uniques while preserving order
            order = [shorts.index(s_pick)] + [idx for idx in prob_order if shorts[idx] != s_pick]
            cands_short_order = [shorts[idx] for idx in order]
            cands_full_order = [short_to_full[s] for s in cands_short_order]

            # Enforce feasibility at this soma
            feas = mtypes_feasible_at_soma(layer_name, z_soma, cands_full_order)
            if feas:
                mfull = feas[0]               # first feasible in our biased order
                mshort = mfull.split("_", 1)[1]
            else:
                # Fallback: nearest feasible among all candidates
                mfull, dist = nearest_feasible_mtype(layer_name, z_soma, cands_full_order)
                if mfull is None:
                    # No interval configured at all; last resort: keep the sampled type
                    mfull = cands_full_order[0]
                    if verbose:
                        print(f"[warn] {layer_name} z={z_soma:.1f}: no feasible intervals; using {mfull}.")
                    mshort = mfull.split("_", 1)[1]
                else:
                    if verbose and dist > 0:
                        print(f"[warn] {layer_name} z={z_soma:.1f}: using nearest feasible {mfull} (Δ≈{dist:.1f} µm).")
                    mshort = mfull.split("_", 1)[1]

            m_type_list.append(mshort.split('_')[0])
            e_type_list.append(mshort.split('_')[1])
            chosen_full.append(mfull)


        # Write outputs
        df["layer"] = layer_name
        df["m_type"] = m_type_list                 # short (e.g., 'LBC_cAC')
        df["e_type"] = e_type_list
        df["desc"] = chosen_full              # full (e.g., 'L23_LBC_cAC')
        df["lay_m_type"] = df["layer"].astype(str) + "_" + df["m_type"].astype(str)
        df["rot_ang"] = rot_angles

        out_path = os.path.join(output_dir, f"cx02_{layer_name}_me_pop.csv")
        df.to_csv(out_path, index=False)
        if verbose:
            print(f"s02: Saved {layer_name}_me_pop.csv (m_type assigned with spatial constraints)")