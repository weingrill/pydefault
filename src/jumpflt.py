#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pyfits
import numpy as np
import matplotlib

matplotlib.use('WXAgg')
import matplotlib.pyplot as plt 

# from matplotlib.backends.backend_pdf import PdfPages
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
time = time[i]
flux = flux[i]

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
df = np.diff(flux) #/ np.diff(time)
sdf = np.std(df)
# max = a if (a > b) else b;
#df1 = map(lambda x: x if (abs(x) < 3*sdf) else math.copysign(3*sdf,x), df)
df1 = (1-exp(-abs(df)/sdf))*np.copysign(sdf,df)
y = np.add.accumulate(df1)
y = np.insert(y, 0, 0)

plt.plot(flux-1,'r')
plt.plot(y)
plt.show()
# flux2 = scipy.signal.medfilt(flux1, 15)

t = time2 - time2[0]

pdf = PdfPages('jmpflt.pdf')
rc('font', size=16)
figure(figsize=(2*8/2.54,2*6/2.54))
axes([0.17,0.15,.8,.75])

# plt.plot(time,flux,',')
title(hdr['COROTID'], fontsize=20)
xlabel('days')
ylabel('norm. flux')
plt.plot(t, flux2,'g')
plt.plot(t, p(time2),'b--')
plt.plot(t, flux3,'r')

# plt.show()
pdf.savefig() # here's another way - or you could do pdf.savefig(1)
pdf.close()
