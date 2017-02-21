import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np
from LowPass import LowPass

class Integrator(object):
	def __init__(self, fs):
		self.lasty = 0
		self.twofs = 2.0*fs
	
	def work(self, value):
		y = value+self.lasty
		self.lasty = value+y
		return y/self.twofs
	
	def value(self):
		return self.lasty/self.twofs

pts = 4000
fs = 44.1E3
fc = 2E3

biq_qual = 1.0/np.sqrt(2.0)
biq_fcut = 550
pll_zeta = 1.0/np.sqrt(2.0)
pll_fnat = 200

G = 4.0*np.pi*pll_zeta*pll_fnat
a = pll_fnat*np.pi/pll_zeta

lpf1 = LowPass(biq_fcut, fs, biq_qual)
lpf2 = LowPass(biq_fcut, fs, biq_qual)
lpf3 = LowPass(biq_fcut, fs, biq_qual)
lpf4 = LowPass(biq_fcut, fs, biq_qual)

s0 = 0.0
s1 = np.pi/2.0
s2 = np.pi
s3 = 3.0*np.pi/2.0

nstep = int(pts/2)
times = np.arange(pts)/fs
step  = np.append(np.zeros(nstep), s2*np.ones(pts-nstep))
phase = 2.0*np.pi*fc*times + step


sig = np.cos(phase)
sig_pwr = np.average(np.power(sig, 2.0))
snr_db = 6 
snr = np.power(10.0, snr_db/10.0)
noise_pwr = sig_pwr/snr
print sig_pwr
print noise_pwr
sig += noise_pwr*np.random.randn(len(sig))

in_phase = np.zeros(pts)
qu_phase = np.zeros(pts)
raw = np.zeros(pts)
err = np.zeros(pts)
inc = 2.0*np.pi*fc
vco = Integrator(fs)
amp = Integrator(fs)

for n in xrange(pts):
	ph = vco.value()
	x1 = 2.0 * np.cos(ph)
	x2 = 2.0 * np.cos(ph + np.pi/4.0)
	x3 = 2.0 * np.cos(ph + np.pi/2.0)
	x4 = 2.0 * np.cos(ph + np.pi*3.0/4.0)

	o1 = lpf1.work(x1 * sig[n])
	o2 = lpf2.work(x2 * sig[n])
	o3 = lpf3.work(x3 * sig[n])
	o4 = lpf4.work(x4 * sig[n])
	s2 = o1 * o2 * o3 * o4
	s3 = G*s2
	s4 = a*s3
	s5 = amp.work(s4)
	s6 = s3+s5
	vco.work(inc+s6)
	in_phase[n] = o1
	qu_phase[n] = o3
	raw[n] = x1/2.0
	err[n] = s2

plt.figure(1)

plt.plot(sig)
plt.plot(raw)

plt.figure(2)
plt.plot(in_phase, qu_phase, 'x')
plt.gca().set_ylim((-1,1))
plt.gca().set_xlim((-1,1))

plt.figure(3)
plt.plot(err)

plt.show()
