#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pyfits
import numpy as np
import matplotlib

matplotlib.use('WXAgg')
import matplotlib.pyplot as plt 

from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas

from pylab import *
from scipy.interpolate import interp1d
import scipy.signal


# hdulist = pyfits.open('/work1/jwe/IC4756/data/LRc06_E2_0155.fits')
hdulist = pyfits.open('/work2/jwe/CoRoT/LRc04/data/LRc04_E2_3913.fits')
hdr = hdulist[0].header
tbdata = hdulist[1].data
time = tbdata.field('datehel')
status = tbdata.field('status')
flux = tbdata.field('redflux')
hdulist.close()

i = np.where(status == 0)
time = time[i] #filter(lambda i: status[i] == 0, time)
flux = flux[i] #filter(lambda i: status[i] == 0, flux)

m = np.median(flux)
flux /= m

f = interp1d(time, flux)
# f2 = interp1d(x, y, kind='cubic')
time1 = np.arange(time[0], time[-1], 32./86400.)
flux1 = f(time1)

# f1 = np.fft.fft(flux1 - flux1.mean())
# q1 = fftfreq(len(time1),32./86400.)

# data = numpy.random.random(100)
# bins = linspace(0, 1, 10)
# digitized = numpy.digitize(data, bins)
# bin_means = [data[digitized == i].mean() for i in range(1, len(bins))]

data = time
bins = np.arange(time[0], time[-1], 0.1)
digitized = np.digitize(data, bins)
dig = list(set(digitized))
# time2 = [time[digitized == i].mean() for i in range(1, len(bins))]
time2 = bins[0:-1]
flux2 = [np.median(flux[digitized == i]) for i in dig]
# flux2 = [flux[digitized == i].mean() for i in range(1, len(bins))]

pf = np.polyfit(time2, flux2, 2)
p = np.poly1d(pf)
flux3 = flux2 - p(time2) + 1.0

# filter derivatives:
# df = np.diff(flux) / np.diff(time)


# flux2 = scipy.signal.medfilt(flux1, 15)

t = time2 - time2[0]

# pdf = PdfPages('detrend.pdf')

# plt.show()
# pdf.savefig() # here's another way - or you could do pdf.savefig(1)
# pdf.close()
matplotlib.rcParams.update({'font.size': 8})
fig = Figure(figsize=(6,4))
canvas = FigureCanvas(fig) # Create a canvas and add the figure to it.
ax = fig.add_subplot(111) # Create a subplot.
ax.set_title(hdr['COROTID'], fontsize=10) # Set the title.
ax.set_xlabel('days',fontsize=8) # Set the X Axis label.
ax.set_ylabel('norm. flux',fontsize=8) # Set the Y Axis label.
ax.grid(True,linestyle='-',color='0.75') # Display Grid.
ax.plot(t, flux2,'g')
ax.plot(t, p(time2),'b--')
ax.plot(t, flux3,'r')

# Save the generated Scatter Plot to a PNG file.
canvas.print_figure('detrend1.pdf',dpi=300)
