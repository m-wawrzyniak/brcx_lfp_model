from neuron import h

from source.src3_utils.ut0_random_manager import rng, np

import config_templates.conf0_model_parameters as conf0

"""
Each PrV will have the same paradigm, either rest, PPI or FDD sharing onset, interval, etc.
For each PrV cell, the paradigm will be individual described by resting and onset states in terms of mean FR and STD.
For each repetition of resting and onset states, specific times of spikes will be created according to Poisson dist.
"""

class PrvCell:
    """
    Klasa PrV_Cell.
    Funkcjonalności:
        - łączenie się innymi komórkami (???) przez synapse typu (???)
        - AP w stanie spoczynku: pewien FR, który ma być z rozkładu normalnego o danej średniej i odchyleniu.
        - AP w stanie pojawienia się bodźca: pewien krótki okres wysokiego FR, wybranego z NORM o średniej i odchyleniu

    """
    gid_prv = 0
    fr_fact = conf0.PRV_FR_FACTOR
    fr_std = conf0.PRV_FR_STD_FACTOR
    param_dists = {
        'stim_on_fr': (2.65*50*fr_fact, fr_std*1.6*50*fr_fact),  # FRs were defined for 20ms window, so to get in Hz, has to be multiplied
        'rest_fr': (11.90*fr_fact, fr_std*12.7*fr_fact)
    }

    def __init__(self, rgen, cell_id=None):

        # Identification objects:
        self.cell_id = cell_id if cell_id is not None else f"PrV_{conf0.GLOBAL_PRV_CELL_CNT}"
        conf0.GLOBAL_PRV_CELL_CNT += 1

        # Simulation objects:
        self.out_nc = []  # Outgoing NetCons
        self._rgen = rgen  # Used random generator - will probably be shared by all PrV cells
        self._strong_fr, self._weak_fr, self._rest_fr = self.pick_frs()  # Picking FRs, based on empirical dist.

        # Registration objects:

    def __repr__(self):
        return str(self.cell_id)

    def set_spike_sequence(self, paradigm_type=None, paradigm_subtype=None):
        """
        Sketch:
            - take correct parameters of the paradigm from simulation_paradigms (dict)
            - for this, pick resting and onset firing rate
            - construct spike-event vector based on chosen paradigm and FRs

        :param paradigm:
        :return:
        """
        ppi_params = conf0.WHISKER_STIMULATION_PARADIGMS[paradigm_type][paradigm_subtype]

        curr_t = 0
        spike_seq = np.array([])
        for stage, params in ppi_params.items():
            if params[0] == 'r':
                rate = self._rest_fr
            elif params[0] == 'wk':
                rate = self._weak_fr
            elif params[0] == 'str':
                rate = self._strong_fr
            else:
                raise ValueError('Paradigm stage labels incorrect')
            subseq = create_sp_arr(period=params[1], start_t=curr_t, rate=rate)
            spike_seq = np.concatenate([spike_seq, subseq])

            curr_t += params[1]


        return spike_seq


    def pick_frs(self):
        wk_scale = conf0.PRV_WEAK_FACTOR
        str_fr = self._rgen.normal(self.param_dists['stim_on_fr'][0], self.param_dists['stim_on_fr'][1])
        str_fr = max(str_fr, 0)
        wk_fr = self._rgen.normal(self.param_dists['stim_on_fr'][0] * wk_scale, self.param_dists['stim_on_fr'][1] * wk_scale)
        wk_fr = max(wk_fr, 0)
        rest_fr = self._rgen.normal(self.param_dists['rest_fr'][0], self.param_dists['rest_fr'][1])
        rest_fr = max(rest_fr, 0)
        return str_fr, wk_fr, rest_fr


def create_sp_arr(period, start_t, rate):

    isi_arr = rng.exponential(1000/rate,  # mean ISI
                                    int(2 * rate * period / 1000))  # Expected number of spikes, times 2, in order not to run out.
    spike_times_arr = np.cumsum(isi_arr) + start_t  # Times of each spike, shifted by start time of activity
    spike_times_arr = spike_times_arr[spike_times_arr < (period + start_t)]

    return spike_times_arr