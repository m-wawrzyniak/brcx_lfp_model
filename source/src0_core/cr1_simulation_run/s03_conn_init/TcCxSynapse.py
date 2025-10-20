from neuron import h

import config_templates.conf01_simulation_parameters as conf01

class TcCxSynapse:
    """
    TC-CX synapse. Implemented as h.ProbAMPANMDA model.

    Globals:
        - GLOBAL_SYN_CNT
        -

    Attributes:
        hsyn (h.HocObject): NEURON synapse object.
        nc (h.NetCon): NetCon between pre_cell.soma(0.5) and hsyn.
        activation_times (h.Vector): Used to store when synapse has been activated.
        pre_id (str): e.g. 'v_0'
        post_id (str): e.g. '0'
        pre_me_type (str): This will always be 'VPM_standard_standard'
        post_me_type (str): e.g. 'L4_SS_cAD'
        post_loc_rel (tuple): Position of the synapse relative to post_cell.soma[0]
        syn_type (str): 'e' or 'i'
        u_val (float):
        d_val (float):
        f_val (float):
        risetime_val (float):
        decay_val (float):
        gsyn_val (float):
        late_comp_ratio (float):

    """
    def __init__(self, syn_obj, netcon:h.NetCon,
                 pre_id:str, post_id:str, pre_me_type:str, post_me_type:str,
                 post_loc_rel:tuple, syn_type:str,
                 u_val:float, d_val:float, f_val:float, risetime_val:float, decay_val:float, gsyn_val:float, late_comp_ratio:float):


        # Synapse ID
        self.syn_id = f"tccx_{conf01.GLOBAL_SYN_CNT}"
        conf01.GLOBAL_SYN_CNT += 1

        # NEURON objects
        self.hsyn = syn_obj
        self.nc = netcon

        # Metadata
        self.pre_id = pre_id
        self.post_id = post_id
        self.pre_me_type = pre_me_type
        self.post_me_type = post_me_type
        self.post_loc_rel = tuple(post_loc_rel)
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
        return f"{self.syn_id}"

    def record_activation(self):
        """
        Sets synapse activation recording.
        """
        self.activation_times = h.Vector()
        self.nc.record(self.activation_times)

    def get_activation_times(self):
        """
        Returns the activation times of the synapse.
        Returns:
            list(float): Times of synapse activation [ms]
        """
        return list(self.activation_times)

    def get_params(self):
        """
        Returns all parameters as a dictionary for easy logging.
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