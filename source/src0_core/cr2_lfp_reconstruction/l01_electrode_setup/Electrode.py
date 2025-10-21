class Electrode:
    """
    Electrode class, meant to reconstruct registration sites in Masmanidis electrode.
    It has topology and electrical parameters important during registration.



    Globals:
        NONE

    Attributes:
        site_topology (dict) : Each site represented by <site_id:(x, y, z>
        sigma (float): Extracellular conductivity [siemens/meter]. For brain tissue, it's about 0.3

    """
    def __init__(self, z_offset:float = 5, topo_variant:str = 'center_col_sparse', cross_species_scale=1):
        self.site_topology = self._setup_electrode_topology(variant=topo_variant, z_offset=z_offset, scale=cross_species_scale)
        self.sigma = 0.3

    def _setup_electrode_topology(self, variant: str, z_offset: float, scale: float):
        if variant == 'full':
            topo = {}
            # Left column (l_n)
            for i in range(21):
                x, y, z = -20, 0, z_offset - i * 50
                topo[f"l_{i}"] = (x * scale, y * scale, z * scale)

            # Central column (c_n)
            for i in range(21):
                x, y, z = 0, 0, z_offset - 25 - i * 50
                topo[f"c_{i}"] = (x * scale, y * scale, z * scale)

            # Right column (r_n)
            for i in range(21):
                x, y, z = 20, 0, z_offset - i * 50
                topo[f"r_{i}"] = (x * scale, y * scale, z * scale)

            return topo

        elif variant == 'center_col_sparse':
            raw_topo = {
                'c_0': (0, 0, -25 + z_offset),
                'c_2': (0, 0, -125 + z_offset),
                'c_4': (0, 0, -225 + z_offset),
                'c_6': (0, 0, -325 + z_offset),
                'c_8': (0, 0, -425 + z_offset),
                'c_10': (0, 0, -525 + z_offset),
                'c_12': (0, 0, -625 + z_offset),
                'c_14': (0, 0, -725 + z_offset),
                'c_16': (0, 0, -825 + z_offset),
                'c_18': (0, 0, -925 + z_offset),
                'c_20': (0, 0, -1025 + z_offset),
            }
            # Apply scaling
            topo = {k: (x * scale, y * scale, z * scale) for k, (x, y, z) in raw_topo.items()}
            return topo

        else:
            raise ValueError('No such variant in MasmanidisEl topology.')