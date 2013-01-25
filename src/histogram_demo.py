#!/usr/bin/env python
# -*- coding: utf-8 -*-
# This is a demo of creating a pdf file with several pages.

import numpy as np
import matplotlib
from matplotlib.backends.backend_pdf import PdfPages
from pylab import *

# Create the PdfPages object to which we will save the pages:
pdf = PdfPages('multipage_pdf.pdf')

rc('text', usetex=True)
rc('font', size=8)
figure(figsize=(8/2.54,6/2.54))
x = np.arange(0,5,0.1)
axes([0.1,0.15,.8,.8])
plot(x, np.sin(x), 'b-')
plot(x, np.cos(x), 'r-')
xlabel('time')
#title('Page Two', fontsize=10)
pdf.savefig() # here's another way - or you could do pdf.savefig(1)
close()


# Remember to close the object - otherwise the file will not be usable
pdf.close()