import config_templates.conf01_simulation_parameters as conf01

class CxSynapse:
    """
    CxSynapse based on two synapse templates provided in (Markram, 2015) - ProbAMPANMDA and ProbGABAAB.

    Globals:
        GLOBAL_SYN_CNT (int) - Used for indexing synapses globally.

    Attributes:
        syn_id (str): Synapse id e.g. 'cxcx_0'
        hsyn (h.PointProcess): Synapse object of NEURON, connects pre_id and post_id cells.
        pre_id (str): Always cortical id.
        post_id (str): Always cortical id.
        pre_me_type (str): e.g. 'L4_SS_cAD'
        post_me_type (str): e.g. 'L4_SS_cAD'
        post_loc_rel (tuple): Synapse (x, y, z) position relative to post_cell.soma[0]
        syn_type (str): 'e' or 'i'
        u_val (float): ???
        d_val (float): ???
        f_val (float): ??/
        risetime_val (float): ???
        decay_val (float): ???
        gsyn_val (float): ???
        late_comp_ratio (float): ???
        activation_times (h.Vector): Saving synapse activation times during simulation.
        nc (h.NetCon): Used for activation detection.
    """
    def __init__(self, syn_obj, netcon, pre_id, post_id, pre_me_type, post_me_type,
                 post_loc_rel, xyz,
                 syn_type, u_val, d_val, f_val, risetime_val, decay_val,
                 gsyn_val, late_comp_ratio,
                 syn_id = None):

        # Synapse ID
        if syn_id is None:
            self.syn_id = f"cxcx_{conf01.GLOBAL_SYN_CNT}"
            conf01.GLOBAL_SYN_CNT += 1
        else:
            self.syn_id = syn_id

        # NEURON objects
        self.hsyn = syn_obj
        self.nc = netcon

        # Metadata
        self.pre_id = pre_id
        self.post_id = post_id
        self.pre_me_type = pre_me_type
        self.post_me_type = post_me_type
        self.post_loc_rel = tuple(post_loc_rel)
        self.xyz = xyz
        self.syn_type = syn_type  # 'e' or 'i'

        # Sampled parameters
        self.u_val = u_val
        self.d_val = d_val
        self.f_val = f_val
        self.risetime_val = risetime_val
        self.decay_val = decay_val
        self.gsyn_val = gsyn_val
        self.late_comp_ratio = late_comp_ratio

        # Optional spike recorder
        self.activation_times = None

    def __repr__(self):
        return f"CxSynapse({self.syn_id})"


    def get_params(self):
        """
        Returns all stored attributes. Used during data dumping after simulation for LFP signal reconstruction.

        Returns:
            dict : A lot.
        """
        return {
            "syn_id": self.syn_id,
            "pre_id": self.pre_id,
            "post_id": self.post_id,
            "pre_me_type": self.pre_me_type,
            "post_me_type": self.post_me_type,
            "post_loc_rel": self.post_loc_rel,
            "syn_type": self.syn_type,
            "u_val": self.u_val,
            "d_val": self.d_val,
            "f_val": self.f_val,
            "risetime_val": self.risetime_val,
            "decay_val": self.decay_val,
            "gsyn_val": self.gsyn_val,
            "late_comp_ratio": self.late_comp_ratio,
        }
