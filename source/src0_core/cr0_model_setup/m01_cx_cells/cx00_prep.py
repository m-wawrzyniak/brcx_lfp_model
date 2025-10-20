
'''

# TODO: Used to limit cell me-type placement
def plot_dend_axon_z_dist(
    cell_templates_path: str = CELL_TEMPLATES_PATH,
    rotation=(np.pi/2, 0, 0),
    bins: int = 40,
    density: bool = True,
    mode: str = "stacked",
    outdir_name: str = "/home/mateusz-wawrzyniak/PycharmProjects/som_sens_processing/m01_brcx_setup/z_distributions",
    use_global_range: bool = True,
    verbose: bool = True
):
    """
    Goes over all cell template files:
        temp_file > 'morphology' > *.asc file.
    and for each creates a z-dist plot with title being 'temp_file' and saves it to CELL_TEMPLATES_PATH/z_distributions/*.jpg.
    """
    os.makedirs(os.path.join(cell_templates_path, outdir_name), exist_ok=True)

    # Collect candidate asc files (first per template folder)
    tasks = []
    for entry in sorted(os.listdir(cell_templates_path)):
        template_dir = os.path.join(cell_templates_path, entry)
        if not os.path.isdir(template_dir):
            continue
        morph_dir = os.path.join(template_dir, "morphology")
        if not os.path.isdir(morph_dir):
            if verbose:
                print(f"[skip] no morphology/ in {entry}")
            continue
        asc_files = sorted(glob.glob(os.path.join(morph_dir, "*.asc")))
        if not asc_files:
            if verbose:
                print(f"[skip] no .asc in {morph_dir}")
            continue
        tasks.append((entry, asc_files[0]))

    # Optional: compute a global z-range so all histograms share the same vertical extent
    bin_range = None
    if use_global_range and tasks:
        zmins, zmaxs = [], []
        for _, asc in tasks:
            cell = CellTopology(morph_path=asc, rotation=rotation)
            zs = []
            if len(cell.axon):     zs.append(np.asarray(cell.axon)[:,2])
            if len(cell.dendrite): zs.append(np.asarray(cell.dendrite)[:,2])
            if not zs:             continue
            allz = np.concatenate(zs)
            zmins.append(float(np.min(allz)))
            zmaxs.append(float(np.max(allz)))
        if zmins and zmaxs:
            bin_range = (min(zmins), max(zmaxs))

    saved = []
    outdir = os.path.join(cell_templates_path, outdir_name)

    for template_name, asc_path in tasks:
        try:
            cell = CellTopology(morph_path=asc_path, rotation=rotation)

            which = []
            if len(cell.axon): which.append("axon")
            if len(cell.dendrite): which.append("dendrite")
            if not which:
                if verbose:
                    print(f"[skip] {template_name}: no axon/dendrite points")
                continue

            # Use your existing plotting method
            fig, ax = cell.plot_z_distribution(
                which=tuple(which),
                bins=bins,
                bin_range=bin_range,
                density=density,
                mode=mode
            )
            ax.set_title(template_name)

            out_path = os.path.join(outdir, f"{template_name}.jpg")
            fig.savefig(out_path, dpi=200, bbox_inches="tight")
            plt.close(fig)

            saved.append((template_name, out_path))
            if verbose:
                print(f"[ok] saved {out_path}")

        except Exception as e:
            if verbose:
                print(f"[err] {template_name}: {e}")

    return saved

'''
"""
plot_dend_axon_z_dist(cell_templates_path=CELL_TEMPLATES_PATH,
        rotation=(np.pi / 2, 0, 0),
        bins=40,
        density=True,
        mode="stacked")
"""