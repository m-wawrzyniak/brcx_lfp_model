COMMENT
IClamp which injects current with specific period
ENDCOMMENT

NEURON {
	POINT_PROCESS IClampFreq
	RANGE del, freq, amp, total_dur, duty_cycle, i
	ELECTRODE_CURRENT i
}
UNITS {
	(nA) = (nanoamp)
}

PARAMETER {
	del = 0 (ms)
	freq = 3 (Hz)
	amp = 0.0 (nA)
	total_dur = 100(ms)
	duty_cycle = 0.5
}
ASSIGNED { 
	i (nA) 
	period (ms)
	pulse_dur (ms)

}

INITIAL {
	i = 0
	if (freq > 0) {
		period = 1000 / freq
		pulse_dur = period * duty_cycle
	} else {
		period = 1e9
		pulse_dur = 0
	}
}

BREAKPOINT {

	if (t < del + total_dur && t >= del) {
		LOCAL local_time, mod_result

		local_time = t - del
		mod_result = local_time - period * floor(local_time / period)
		
		if (mod_result < pulse_dur) {
			i = amp
		} else {
			i = 0
		}
	} else {
		i = 0
	}
}