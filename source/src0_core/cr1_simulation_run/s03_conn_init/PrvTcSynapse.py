from neuron import h
import config_templates.conf01_simulation_parameters as conf01

class PrvTcSynapse:
    """
    Minimalistic PrV -> VPM excitatory synapse
    Fast AMPA kinetics + short-term depression.
    """
    def __init__(self, syn_obj, nc, vestim, spike_vec, weight=0.1, u=0.7, tau_rec=800.0):
        """
        Parameters
        ----------
        weight : float
            Initial synaptic weight (uS for Exp2Syn)
        u : float
            Release probability
        tau_rec : float
            Recovery time from depression (ms)
        """

        # Synapse ID
        self.syn_id = f"prvtc_{conf01.GLOBAL_SYN_CNT}"
        conf01.GLOBAL_SYN_CNT += 1

        # Other
        self.nc = nc
        self.vecstim = vestim
        self.spike_vec = spike_vec

        # Synapse kinetics
        self.hsyn = syn_obj
        self.hsyn.tau1 = 0.2   # ms, rise time
        self.hsyn.tau2 = 1.5   # ms, decay time
        self.hsyn.e = 0        # mV, reversal potential for glutamate

        # STP parameters
        self.U = u
        self.tau_rec = tau_rec
        self.last_spike = -1e9
        self.R = 1.0  # fraction of resources available

        # Base weight
        self.base_weight = weight

    def _on_spike(self):
        """Callback called on presynaptic spike to apply depression."""
        t = h.t
        dt = t - self.last_spike
        # Recover resources
        self.R += (1.0 - self.R) * (1 - np.exp(-dt / self.tau_rec))
        # Resources used this spike
        used = self.U * self.R
        self.R -= used
        # Scale weight accordingly
        self.hsyn._ref_weight[0] = self.base_weight * used
        self.last_spike = t