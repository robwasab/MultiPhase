import numpy as np

class LowPass(object):
	def __init__(self, fc, fs, Q = 0.7071):
		#Biquad filter
		self.a0 = self.a1 = self.a2 = 0
		self.b0 = self.b1 = self.b2 = 0
		self.x1 = self.x2 = 0
		self.y1 = self.y2 = 0
		A = 1.0
		omega = 2.0*np.pi*fc/fs
		sn = np.sin(omega)
		cs = np.cos(omega)
		alpha = sn / (2.0*Q)
		beta  = np.sqrt(A+A)
		self.init_lp(A, omega, sn, cs, alpha, beta)
		self.b0 /= self.a0
		self.b1 /= self.a0
		self.b2 /= self.a0
		self.a1 /= self.a0
		self.a2 /= self.a0

	def init_lp(self, A, omega, sn, cs, alpha, beta):
		self.b0 = (1.0 - cs) / 2.0
		self.b1 =  1.0 - cs
		self.b2 = (1.0 - cs) / 2.0
		self.a0 =  1.0 + alpha
		self.a1 = -2.0 * cs
		self.a2 =  1.0 - alpha
	
	def work(self, x):
		y = self.b0 * x + self.b1 * self.x1 + self.b2 * self.x2 - self.a1 * self.y1 - self.a2 * self.y2
		self.x2, self.x1 = self.x1, x
		self.y2, self.y1 = self.y1, y
		return y

	def reset(self):
		self.x2 = 0
		self.x1 = 0
		self.y2 = 0
		self.y1 = 0

if __name__ == '__main__':
	import matplotlib
	matplotlib.use('TkAgg')
	import matplotlib.pyplot as plt
	from numpy.fft import fft

	fs = 44.1E3
	lpf = LowPass(500, fs)
	sig = np.zeros(4096)
	out = np.zeros(4096)
	sig[0] = 1.0

	for k in range(len(sig)):
		out[k] = lpf.work(sig[k])

	freqs = np.arange(len(sig)) / float(len(sig)) * fs
	sig_fft = fft(sig)
	out_fft = fft(out)
	resp = out_fft / sig_fft

	plt.subplot(121)
	plt.plot(freqs, 20.0*np.log10(np.abs(resp)))

	plt.subplot(122)
	plt.plot(freqs, np.angle(resp)*180.0/np.pi)

	plt.show()
