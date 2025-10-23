import os
import pandas as pd
from scipy.stats import qmc

from source.src2_utils.ut0_random_manager import np, rng
import config_templates.conf0_model_parameters as conf0


def nrn_positions_cyl(layer_name: str, show_info: bool = False,
                      x_y_std:float = conf0.X_Y_JITTER_STD) -> np.ndarray:
    print(f"[cx01] Generating cx positions at {layer_name}")

    start_H = conf0.TISSUE_PARAMS[layer_name]['start_H']
    R = conf0.TISSUE_PARAMS[layer_name]['R']
    H = conf0.TISSUE_PARAMS[layer_name]['H']
    D = conf0.TISSUE_PARAMS[layer_name]['D'] / 1e6  # convert to nrn/um³
    M = max(1, conf0.TISSUE_PARAMS[layer_name]['M'])
    correction_factor = 4/np.pi  # correction bcs I discard calls outside cylinder, and Sobol creates rec.prism V_rec / V_cyl for R, H = const.

    v_cyl = np.pi * R**2 * H
    nrn_n = max(1, int(D * v_cyl * correction_factor * conf0.SCALE_FACTOR))
    nrn_per_col = max(1, int(nrn_n / M))

    # Sobol-distributed space-filling inside rectangular prism, then filtered by circle radius
    oversamp_fac = 4  # oversampling for discarding too close somatas
    sobol_seed = rng.integers(0, 2 ** 32 - 1)
    sampler = qmc.Sobol(d=2, scramble=True, seed=sobol_seed)
    raw_samp = sampler.random(n=M * oversamp_fac)
    col_x = raw_samp[:, 0] * 2 * R - R
    col_y = raw_samp[:, 1] * 2 * R - R
    col_xy = np.column_stack((col_x, col_y))

    # Making sure that chosen x,y are within the cylinder
    def _in_cylinder(x, y, R):
        return x**2 + y**2 <= R**2
    mask = np.array([_in_cylinder(xi, yi, R) for xi, yi in col_xy])
    col_xy = col_xy[mask][:M]

    # Deciding on z-coord for each chosen minicolumn
    nrn_pos = []
    layer_top = start_H
    layer_bottom = start_H - H  # since start_H is negative, this is more negative

    for x, y in col_xy:
        # Apply Gaussian jitter to x and y
        x_jittered = x + np.random.normal(0, x_y_std)
        y_jittered = y + np.random.normal(0, x_y_std)

        if nrn_per_col >= 5:  # distribute evenly and jitter z
            z_base = np.linspace(layer_top, layer_bottom, nrn_per_col)
            z_jitter = np.random.normal(0, conf0.Z_STD, size=nrn_per_col)
            z_vals = np.clip(z_base + z_jitter, layer_bottom, layer_top)
        else:  # few cells, uniform distribution in z
            #print("s01: <5 nrn per column: using uniform z-sampling")
            z_vals = np.random.uniform(layer_bottom, layer_top, size=nrn_per_col)
        z_vals.sort()

        # Ensure MIN_SOMA_DISTANCE separation in z
        accepted_z = []
        for zi in z_vals:
            if not accepted_z or abs(zi - accepted_z[-1]) >= conf0.MIN_SOMA_DISTANCE:
                accepted_z.append(zi)
                if len(accepted_z) == nrn_per_col:
                    break

        # Keep only cells within the cylinder based on jittered x, y
        for zi in accepted_z:
            if _in_cylinder(x_jittered, y_jittered, R):
                nrn_pos.append([x_jittered, y_jittered, zi])

    nrn_pos = np.array(nrn_pos)

    print(f"[cx01] SUCCESS: cx generated for {layer_name}")
    if show_info:
        print(f"\t Column volume: {v_cyl:.2f} µm³")
        print(f"\t Target total nrns: {nrn_n}")
        print(f"\t Target nrns per column: {nrn_per_col}")
        print(f"\t Obtained total nrns: {len(nrn_pos)}")

    return nrn_pos

def save_nrn_positions(nrn_pos: np.ndarray, file_name: str,
                       save_dir: str):
    print(f"[cx01] Saving cx positions for {file_name}")

    df = pd.DataFrame(nrn_pos, columns=["x", "y", "z"])
    df.insert(0, "cell_id", range(conf0.GLOBAL_CX_CELL_CNT, conf0.GLOBAL_CX_CELL_CNT + len(df)))
    conf0.GLOBAL_CX_CELL_CNT += len(df)
    full_path = os.path.join(save_dir, file_name)
    df.to_csv(full_path, index=False)

    print(f"[cx01] SUCCESS: Saved {file_name}")