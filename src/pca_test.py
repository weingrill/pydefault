#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on Jun 1, 2015

@author: Joerg Weingrill <jweingrill@aip.de>
'''
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import svd

stars = 1000
epochs = 100

# (observations, features) matrix
M = np.zeros([epochs, stars])

assert(epochs<stars)

t = np.random.rand(epochs)*30
t = np.sort(t)
#t = np.linspace(0, 30, epochs)

#mag0 = np.arange(stars)*0.0
mag0 = np.random.rand(stars)*8.0 + 7.0
amp = np.random.randn(stars)*0.1
period = np.random.rand(stars)*12 + 0.2
noise = np.random.randn(epochs)*0.2
phase = np.random.rand(stars)*2.0*np.pi

for i in range(stars):
    #if amp[i]>0.2: amp[i]=0.0
    data = amp[i]*np.sin(2*np.pi*t/period[i]+phase[i]) + mag0[i]
    photnoise = np.random.randn(epochs)*0.01
    M[:, i] = data + noise + photnoise


print t.shape #(100,)
print M.shape #(100, 262144)

U, s, Vt = svd(M, full_matrices=False)
V = Vt.T

# sort the PCs by descending order of the singular values (i.e. by the
# proportion of total variance they explain)
ind = np.argsort(s)[::-1]
U = U[:, ind]
s = s[ind]
V = V[:, ind]

# if we use all of the PCs we can reconstruct the noisy signal perfectly
S = np.diag(s)
Mhat = np.dot(U, np.dot(S, V.T))
print "Using all PCs, MSE = %.6G" %(np.mean((M - Mhat)**2))

# if we use only the first n PCs the reconstruction is less accurate
level = 1

k = 100

Mhat2 = np.dot(U[:, :level], np.dot(S[:level, :level], V[:,:level].T))
print "Using first %d PCs, MSE = %.6G" %(level, np.mean((M - Mhat2)**2))

fig, [ax1, ax2, ax3, ax4] = plt.subplots(4, 1)
ax1.scatter(t, M[:, k], edgecolor='none')
x = np.arange(30, step=0.1)
ax1.plot(x, amp[k]*np.sin(2*np.pi*x/period[k]+phase[k]) + mag0[k], 'r')
ax1.axis('tight')

ax2.scatter(t, Mhat[:, k]-noise, edgecolor='none')
ax2.axis('tight')

ax3.scatter(t, M[:, k]-Mhat2[:, k], color='g')
ax3.plot(x, amp[k]*np.sin(2*np.pi*x/period[k]+phase[k]),'r')
ax3.axis('tight')

ax4.semilogy(s[:10])
plt.show()

fig, [ax1, ax2] = plt.subplots(2, 1)
ax1.imshow(M)
ax2.imshow(M-Mhat2)
plt.show()

