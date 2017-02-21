import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from LowPass import LowPass
import numpy as np

pi = np.pi
fs = 44.1E3
fc = 3E3
deg = 0
pts = 1000

gain = 4.0
x1 = gain * np.cos(2.0*pi*fc/fs*np.arange(pts) + 0.0)
x2 = gain * np.cos(2.0*pi*fc/fs*np.arange(pts) + pi/4.0)
x3 = gain * np.cos(2.0*pi*fc/fs*np.arange(pts) + pi/2.0)
x4 = gain * np.cos(2.0*pi*fc/fs*np.arange(pts) + 3.0*pi/4.0)

out = np.zeros(pts)

lpf1 = LowPass(500, fs)
lpf2 = LowPass(500, fs)
lpf3 = LowPass(500, fs)
lpf4 = LowPass(500, fs)

degs = np.arange(100)/float(100) * 2.0 * pi
errs = np.zeros(len(degs))

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
		out[k] = lpf1.work(err1[k]) * lpf2.work(err2[k]) * lpf3.work(err3[k]) * lpf4.work(err4[k])

	errs[n] = np.average(out)

plt.plot(degs/pi, errs)
plt.show()
