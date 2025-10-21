import config_templates.conf01_simulation_parameters as conf1

from neuron import h

def run_simulation():
    h.celsius = conf1.TEMP
    h.dt = conf1.DT
    tstop = conf1.SIM_TIME

    t = h.Vector()
    t.record(h._ref_t)
    h.finitialize()


    dt_chunk = 10  # ms
    while h.t < tstop:
        h.continuerun(h.t + dt_chunk)
        print(f"\t Sim. time: {h.t:.1f} / {tstop:.1f} ms")
    print("M02: Simulation finished successfully.")
    print('\n')

    return t