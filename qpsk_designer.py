import matplotlib
matplotlib.use('TkAgg')
from matplotlib.widgets import Slider, Button, RadioButtons
import matplotlib.pyplot as plt
import numpy as np
import step
import pdb

def biquad(s, q, fc):
	wc = 2.0*np.pi*fc
	return (wc**2.0)/(s**2.0+wc/q*s+wc**2.0)

# loop amp to achieve 
# perfect second order loop
def loop_amp(s, fnat, zeta):
	gain = 4.0*np.pi*zeta*fnat
	a = np.pi*fnat/zeta
	return gain*(1.0+a/s)

def pll(s, loop_fnat, loop_zeta, biquad_qual, biquad_fcut):
	lpf = 4.0*biquad(s, biquad_qual, biquad_fcut)
	amp = loop_amp(s, loop_fnat, loop_zeta)
	return lpf*amp/(s+lpf*amp)

def bilinear_transform(fs, freqs):
	inv_z = np.exp(-2j*np.pi*freqs/fs)
	s = 2.0*fs*(1.0-inv_z)/(1.0+inv_z)
	return s

def plot_mag(info, laplace_resp, digital_resp, txt):
	continu = 10.0*np.log10(np.abs(laplace_resp))
	digital = 10.0*np.log10(np.abs(digital_resp))
	if info['handles'] is None:
		handles = {}
		ax = info['ax']
		handles['laplace'], = ax.semilogx(freqs, continu, ':', label='S-Domain')
		handles['digital'], = ax.semilogx(freqs, digital, color='m', label='Bilinear')
		ax.set_title('%s Loop Magnitude Response'%txt)
		ax.set_xlabel('Frequency [Hz]')
		ax.set_ylabel('Magnitude [u]')
		ax.set_xlim((0,freqs[-1]))
		info['handles'] = handles
	else:
		info['handles']['laplace'].set_ydata(continu)
		info['handles']['digital'].set_ydata(digital)
		ax = info['ax']
		ax.relim()
		ax.autoscale_view()
	#ax.legend(handles=[laplace, digital], loc='upper left', bbox_to_anchor=(1,1))

def plot_ang(info, laplace_resp, digital_resp, txt):
	continu = np.angle(laplace_resp)/np.pi
	digital = np.angle(digital_resp)/np.pi
	if info['handles'] is None:
		handles = {}
		ax = info['ax']
		handles['laplace'], = ax.semilogx(freqs, continu, ':', label = 'S-Domain')
		handles['digital'], = ax.semilogx(freqs, digital, color='m', label = 'Bilinear')
		#ax.set_title('%s Loop Phase Response'%txt)
		ax.set_xlabel('Frequency [Hz]')
		ax.set_ylabel('Phase/$\pi$ [rad]')
		ax.set_xlim((0,freqs[-1]))
		info['handles'] = handles
	else:
		info['handles']['laplace'].set_ydata(continu)
		info['handles']['digital'].set_ydata(digital)
		ax = info['ax']
		ax.relim()
		ax.autoscale_view()
	#ax.legend(handles=[laplace, digital], loc='upper left', bbox_to_anchor=(1,1))

def plot_resp(pll_fnat, pll_zeta, biq_qual, biq_fcut):
	digital = {}
	laplace = {}

	digital['open'] = biquad(z2s, biq_qual, biq_fcut)*loop_amp(z2s, pll_fnat, pll_zeta)
	laplace['open'] = biquad(  s, biq_qual, biq_fcut)*loop_amp(  s, pll_fnat, pll_zeta)

	digital['closed'] = pll(z2s, pll_fnat, pll_zeta, biq_qual, biq_fcut)
	laplace['closed'] = pll(  s, pll_fnat, pll_zeta, biq_qual, biq_fcut)

	plot_mag(axes['mag-closed'], laplace['closed'], digital['closed'], 'Closed')
	plot_ang(axes['ang-closed'], laplace['closed'], digital['closed'], 'Closed')

	plot_mag(axes['mag-open'], laplace['open'], digital['open'], 'Open')
	plot_ang(axes['ang-open'], laplace['open'], digital['open'], 'Open')
	plot_transient(pll_fnat, pll_zeta, biq_qual, biq_fcut)

def setup_slider(pll_fnat, pll_zeta, biq_qual, biq_fcut, pll_fnat_cb, pll_zeta_cb, biq_qual_cb, biq_fcut_cb):
	plt.subplots_adjust(left=0.2, bottom=0.25)
	axcolor = 'lightgoldenrodyellow'
	ax1 = plt.axes([0.2, 0.00, 0.65, 0.03], axisbg=axcolor)
	ax2 = plt.axes([0.2, 0.05, 0.65, 0.03], axisbg=axcolor)
	ax3 = plt.axes([0.2, 0.10, 0.65, 0.03], axisbg=axcolor)
	ax4 = plt.axes([0.2, 0.15, 0.65, 0.03], axisbg=axcolor)

	sfnat = Slider(ax1, 'pll_fnat', 10, fs/2.0, valinit=pll_fnat)
	szeta = Slider(ax2, 'pll_zeta', 0.05, 10.0, valinit=pll_zeta)
	squal = Slider(ax3, 'biq_qual', 0.05, 10.0, valinit=biq_qual)
	sfcut = Slider(ax4, 'biq_fcut', 10, fs/2.0, valinit=biq_fcut)

	sfnat.on_changed(pll_fnat_cb)
	szeta.on_changed(pll_zeta_cb)
	squal.on_changed(biq_qual_cb)
	sfcut.on_changed(biq_fcut_cb)
	return (sfnat, szeta, squal, sfcut)

def plot_transient(pll_fnat, pll_zeta, biq_qual, biq_fcut):
	phivco = step.pseudo_transient(pll_fnat, pll_zeta, biq_qual, biq_fcut, 3000, fs, np.pi)
	phivco/= np.pi
	ax = axes['transient']['ax']
	handle = axes['transient']['handle']
	if handle is None:
		handle, = ax.plot(phivco)
		ax.set_title('Transient Response')
		ax.set_ylabel('VCO Phase/$\pi$')
		axes['transient']['handle'] = handle
		ax.set_ylim((0, 2))
	else:
		handle.set_ydata(phivco)
		ax.relim()
		ax.autoscale_view()

plt.figure(0, facecolor='white')
pts = 500
fs = 44.1E3
freqs = (1+np.arange(pts-1))/float(pts)*fs/2.0
s = 2j * np.pi * freqs
z2s = bilinear_transform(fs, freqs)
pll_fnat = 1E3
pll_zeta = 1.0/np.sqrt(2.0)
biq_qual = 1.0/np.sqrt(2.0)
biq_fcut = 100

axes = {}
axes['mag-closed'] = {'ax':plt.subplot(321), 'handles':None}
axes['ang-closed'] = {'ax':plt.subplot(323), 'handles':None}
axes['mag-open'] =   {'ax':plt.subplot(322), 'handles':None}
axes['ang-open'] =   {'ax':plt.subplot(324), 'handles':None}
axes['transient'] =  {'ax':plt.subplot(313), 'handle':None}

def pll_fnat_cb(val):
	global pll_fnat
	pll_fnat = val
	plot_resp(pll_fnat, pll_zeta, biq_qual, biq_fcut)

def pll_zeta_cb(val):
	global pll_zeta
	pll_zeta = val
	plot_resp(pll_fnat, pll_zeta, biq_qual, biq_fcut)

def biq_qual_cb(val):
	global biq_qual
	biq_qual = val
	plot_resp(pll_fnat, pll_zeta, biq_qual, biq_fcut)

def biq_fcut_cb(val):
	global biq_fcut
	biq_fcut = val
	plot_resp(pll_fnat, pll_zeta, biq_qual, biq_fcut)

sfnat, szeta, sfcut, squal = setup_slider(pll_fnat, pll_zeta, biq_qual, biq_fcut, pll_fnat_cb, pll_zeta_cb, biq_qual_cb, biq_fcut_cb)
plot_resp(pll_fnat, pll_zeta, biq_qual, biq_fcut)
plt.show()
