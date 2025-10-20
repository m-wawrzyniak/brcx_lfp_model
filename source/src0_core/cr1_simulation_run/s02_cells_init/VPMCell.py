from neuron import h

from source.src3_utils.ut0_random_manager import np


class VPMCell:
    """
    VPMCell model, similar to (Destexhe, ????).

    Globals:
        NONE

    Attributes:
        cell_id (str): Cell id.
        cell_name (str): Cell short description: "VPM_standard_standard"
        variant (str): 'follower' or 'initiator'
        spike_times (h.Vector): Spiking times.
        spike_detector (h.NetCon): NetCon at soma(0.5), detecting spikes.
        self.bio_parameters (dict): All biophysics parameters assigned to this cell.
        self.tree (h.SectionList): All sections.
        self.soma (h.Section): Soma of h.cell
        self.dend (h.Section) Single dendrite section.
        self.kleak (h.PointProcess): Potassium leak at soma.
    """
    bio_parameters_dist = {
        "soma": {
            # Passive
            "e_pas": (-70, 0),
            "g_pas": (1e-5, 0),
            "ek": (-100, 0),
            "ena": (50, 0),

            # HH Na/K
            "vtraub_hh2": (-25, 0),
            "gnabar_hh2": (0.08, 0),  # ↓ Na (was 0.09)
            "gkbar_hh2": (0.05, 0),  # ↑ K (was 0.04)

            # Calcium IT
            "cai": (2.4e-4, 0),
            "cao": (2, 0),
            "eca": (120, 0),
            "gcabar_it": (0.0018, 0),  # slight ↓ T current (was 0.002)

            # Ih (IAR)
            "eh": (-40, 0),
            "nca_iar": (4, 0),
            "k2_iar": (0.0004, 0),
            "cac_iar": (0.002, 0),
            "nexp_iar": (1, 0),
            "k4_iar": (0.001, 0),
            "Pc_iar": (0.004, 0),  # ↓ binding affinity (was 0.005)
            "ginc_iar": (0.25, 0),  # further ↓ Ih gain (was 0.5)
            "ghbar_iar": (2.0e-5, 0),  # baseline Ih slightly ↓ (was 2.5e-5)

            # Calcium dynamics
            "depth_cad": (1, 0),
            "taur_cad": (2.5, 0),  # even faster clearance
            "cainf_cad": (2.4e-4, 0),
            "kt_cad": (0, 0),
        },

        "point_proc": {
            "k_leak": {
                "gmax": (0.005, 0.0)  # ↑ leak (was 0.004)
            }
        },

        "global": {
            "Erev_kleak": (-100.0, 0.0)
        }
    }


    def __init__(self, cell_id:str, soma_l=96, soma_diam=96, ra=100, cm=1,
                 ext_bio_params=None, rand_params=False):
        self.cell_id = cell_id
        self.cell_name = 'VPM_standard_standard'

        self.soma_len = soma_l
        self.soma_diam = soma_diam
        self.ra = ra
        self.cm = cm
        self.v_rest = -60.0

        self._set_morpho()
        self.soma.v = self.v_rest

        self.spike_times = h.Vector()
        self.spike_detector = None

        if ext_bio_params is not None:
            self.bio_parameters = ext_bio_params
        else:
            self.bio_parameters = self._set_bio_params(rand_params)

        self._set_biophys()

    def _set_morpho(self):
        self.soma = h.Section(name="soma", cell=self)
        self.soma.L = self.soma_len
        self.soma.diam = self.soma_diam
        self.tree = h.SectionList()
        self.tree.wholetree(self.soma)

    def _set_bio_params(self, rand=False):
        chosen_bio = {}
        for param_type, type_dict in self.bio_parameters_dist.items():
            chosen_bio[param_type] = {}
            if param_type == "point_proc":
                for pp_type, pp_dict in type_dict.items():
                    chosen_bio[param_type][pp_type] = {}
                    for pp_param, pp_dist in pp_dict.items():
                        val = np.random.normal(pp_dist[0], pp_dist[1]) if rand else pp_dist[0]
                        chosen_bio[param_type][pp_type][pp_param] = val
            else:
                for param, dist in type_dict.items():
                    val = np.random.normal(dist[0], dist[1]) if rand else dist[0]
                    chosen_bio[param_type][param] = val
        return chosen_bio

    def _set_biophys(self):
        for sec in self.tree:
            sec.Ra = self.ra
            sec.cm = self.cm

        self.soma.insert("pas")
        self.soma.insert("hh2")
        self.soma.insert("cad")
        self.soma.insert("iar")
        self.soma.insert("it")

        for param_name, val in self.bio_parameters["soma"].items():
            setattr(self.soma, param_name, val)

        self.k_leak = h.kleak(self.soma(0.5))
        for _, pp_dict in self.bio_parameters["point_proc"].items():
            for param_name, val in pp_dict.items():
                setattr(self.k_leak, param_name, val)

        for param_name, val in self.bio_parameters["global"].items():
            setattr(h, param_name, val)

    def setup_spike_detector(self, threshold:float = 0.0):
        """
        Sets up a NetCon to detect APs in the soma(0.5).

        Args:
            threshold (float): Potential value at which spike is registered.
        """
        soma_seg = self.soma(0.5)
        self.spike_detector = h.NetCon(soma_seg._ref_v, None, sec=self.soma)
        self.spike_detector.threshold = threshold
        self.spike_detector.record(self.spike_times)

    def get_spike_times(self) -> list:
        """
        Returns spike_times after the simulation.
        """
        return list(self.spike_times)

