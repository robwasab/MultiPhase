import matplotlib.pyplot as plt
import numpy as np
import lpf

#   phases    s1           s2                s3
#INPUT----[-]---->[sin( )]---->[*G]----------------------+
#          ^                                             |
#          | phivco                                      |
#          |         s6        Loop Filter Section       |
#          +---[-1/s]<---[+]<----------------------------+
#                         ^                              |
#                         | s5                           |
#                         |          s4a      s4         |
#                         +----[1/s]-----[+]<----[*a]<---+
#                         |               |
#                         | s5            |
#                         |               |
#                         +----[*-a]------+


# Simulation Variables
def pseudo_transient(pll_fnat, pll_zeta, biq_qual, biq_fcut, npts, fs, phase_step):
	twofs = 2.0*fs
	twopi = 2.0*np.pi

	# Loop Variables
	G = 4.0*np.pi*pll_zeta*pll_fnat
	a = pll_fnat*np.pi/pll_zeta
	phivco = 0

	# Direct Form II Variables
	w2b = 0
	w2c = 0

	# Generate Arrays
	nstep = int(npts*0.1)
	phases = np.append(np.zeros(nstep), phase_step*np.ones(npts-nstep))
	errors = np.zeros(npts)
	fvco = np.zeros(npts)
	pvco = np.zeros(npts)
	filt = lpf.LowPass(biq_fcut, fs, biq_qual)

	for n in xrange(npts):
		pvco[n] = phivco
		s1 = phases[n] - phivco
		s2 = filt.work(np.sin(s1))
		s3 = G*s2
		s4 = a*s3

		w1b = s4+w2b
		w2b = s4+w1b

		s5 = w1b/twofs
		s6 = s3+s5

		w1c = s6+w2c
		w2c = s6+w1c
		phivco = w1c/twofs
		errors[n] = s1
		fvco[n] = s6/twopi
	return pvco

if __name__ == '__main__':
	pll_fnat = 1E3
	pll_zeta = 1.0/np.sqrt(2.0)
	biq_qual = 1.0/np.sqrt(2.0)
	biq_fcut = 100.0
	npts = 2000
	fs   = 44.1E3
	phase_step = np.pi
	phase_resp = pseudo_transient(pll_fnat, pll_zeta, biq_qual, biq_fcut, npts, fs, phase_step)
	plt.plot(phase_resp)
	plt.show()
