import os
import pandas as pd
from scipy.stats import qmc

from source.src3_utils.ut0_random_manager import np, rng
import config_templates.conf0_model_parameters as conf0


def nrn_positions_cyl(layer_name: str, show_info: bool = False,
                      x_y_std:float = conf0.X_Y_JITTER_STD) -> np.ndarray:
    """
    Using Sobol-distributed space-filling algorithm, specifies the positions of cortical cells somata.
    Globals:
        - TISSUE_PARAMS
        - SCALE_FACTOR
        - MIN_SOMA_DISTANCE
        - GLOBAL_SEED
        - Z_STD

    Args:
        layer_name (str): Name of the cortical layer (e.g., 'L23', 'L4', 'L5').
        show_info (bool, optional): If True, prints debug information. Defaults to False.

    Returns:
        np.ndarray: Array of shape (N, 3), where each row is (x, y, z) position of a cell.
    """

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

    if show_info:
        print(f"s01: Specyfing soma positions for {layer_name}...")
        print(f"\t Column volume: {v_cyl:.2f} µm³")
        print(f"\t Target total nrns: {nrn_n}")
        print(f"\t Target nrns per column: {nrn_per_col}")
        print(f"\t Obtained total nrns: {len(nrn_pos)}")

    return nrn_pos

def save_nrn_positions(nrn_pos: np.ndarray, file_name: str,
                       save_dir: str):
    """
    Takes the <nrn_pos> - np.array with x,y,z done by nrn_positions(), adds the cell_id column as the df, and transforms the df into
    .csv file saving it at <path>

    Globals:
        GLOBAL_CX_CELL_ID

    Args:
        nrn_pos (np.array): Output from s01.nrn_positions_cyl(). Somata positions.
        file_name (str): Name of the saved *.csv file.
        save_dir (str): Directory where *csv should be saved.

    Returns:
        None
    """

    df = pd.DataFrame(nrn_pos, columns=["x", "y", "z"])
    df.insert(0, "cell_id", range(conf0.GLOBAL_CX_CELL_CNT, conf0.GLOBAL_CX_CELL_CNT + len(df)))
    conf0.GLOBAL_CX_CELL_CNT += len(df)
    full_path = os.path.join(save_dir, file_name)
    df.to_csv(full_path, index=False)
    print(f"s01: Saved {file_name} - somata positions")




## TODO: LEGACY
def __main__(save_dir, verbose=False):
    """
    1. Goes over each <l> in LAYER_COMP.
    2. Uses nrn_positions() for each <l> - specifying the somata position at each layer.
    3. Saves the somata positions using save_nrn_positions.
    4. Provides FINAL_OUTPUT

    Args:
        verbose (bool): Should the component functions plot/print all the info.
    """
    print(f"s01_population_topology:")

    for l in conf0.LAYER_COMP_TARGET.keys():
        l_pos = nrn_positions_cyl(l, show_info=True)
        save_nrn_positions(l_pos, save_dir=save_dir, file_name=f"cx01_{l}_soma_pos.csv")
    #TODO: create_nrn_pos_summary()

    #TODO: if verbose: show_layer_comp()

# TODO:
"""
if __name__ == "__main__":
    res = calc_goal_densities()
    pprint(res)
"""