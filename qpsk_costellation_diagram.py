import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from LowPass import LowPass
import numpy as np

def plot_const(inp, qup):
	plt.plot(inp, qup, 'x')
	plt.gca().set_xlim((-1,1))
	plt.gca().set_ylim((-1,1))

pi = np.pi
fs = 44.1E3
fc = 3E3
deg = 0
pts = 1000

gain = 2.0
x1 = gain * np.cos(2.0*pi*fc/fs*np.arange(pts) + 0.0)
x2 = gain * np.cos(2.0*pi*fc/fs*np.arange(pts) + pi/4.0)
x3 = gain * np.cos(2.0*pi*fc/fs*np.arange(pts) + pi/2.0)
x4 = gain * np.cos(2.0*pi*fc/fs*np.arange(pts) + 3.0*pi/4.0)

out = np.zeros(pts)
out1 = np.zeros(pts)
out2 = np.zeros(pts)
out3 = np.zeros(pts)
out4 = np.zeros(pts)

lpf1 = LowPass(500, fs)
lpf2 = LowPass(500, fs)
lpf3 = LowPass(500, fs)
lpf4 = LowPass(500, fs)

degs = [pi/4.0, 3.0/4.0*pi, 5.0/4.0*pi, 7.0/4.0*pi]
errs = np.zeros(len(degs))
in_phase = np.zeros((len(degs),pts))
qu_phase = np.zeros((len(degs),pts))

for n in xrange(len(degs)):
	deg = degs[n]
	sig = np.cos(2.0*pi*fc/fs*np.arange(pts) + deg)

	err1 = x1 * sig
	err2 = x2 * sig
	err3 = x3 * sig
	err4 = x4 * sig

	lpf1.reset()
	lpf2.reset()
	lpf3.reset()
	lpf4.reset()

	for k in xrange(pts):
		out1[k] = lpf1.work(err1[k]) 
		out2[k] = lpf2.work(err2[k]) 
		out3[k] = lpf3.work(err3[k]) 
		out4[k] = lpf4.work(err4[k])

	out = out1 * out2 * out3 * out4
	errs[n] = np.average(out)
	in_phase[n,:] = out1
	qu_phase[n,:] = out3

for n in xrange(len(degs)):
	plt.subplot(221 + n)
	plot_const(in_phase[n,:], qu_phase[n,:])
plt.show()
