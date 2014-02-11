'''
Created on Jan 9, 2014

@author: jwe
'''
import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt

fitfunc = lambda p, x: p[0]*np.cos(2*np.pi/p[1]*x+p[2])
# Distance to the target function
errfunc = lambda p, x, y: y - fitfunc(p, x)

func = lambda x, a, p, pha: a*np.cos(2*np.pi/p*x+pha)


n = 1000

p = [0.01, 4.0, 0.1]

x = np.linspace(0.0, 11.0, n)
#y = p[0]*np.cos(2*np.pi/p[1]*x+p[2])
y = func(x, p[0], p[1], p[2])
noise = 0.01*(np.random.rand(n)-0.5)
signal = y + noise

p0 = [0.1, 5.0, 0.5] # Initial guess for the parameters

ntests = 1000

periods = np.empty(ntests)

for i in range(ntests):
    p0[1]= np.random.rand(1)*11.0
#    p1, cov, _, _, success = optimize.leastsq(errfunc, 
#                                              p0[:], 
#                                              args=(x, signal), 
#                                              full_output=1)
    
    p1, cov = optimize.curve_fit(func, x, signal, p0=p0)
    
    e = np.sqrt(np.diag(cov))
    periods[i] = p1[1]  
print p
print p1
print e
print np.mean(periods), np.std(periods)
print e[1]/np.std(periods)

plt.subplot(2, 1, 1)
plt.hist(periods, bins=np.sqrt(ntests))
plt.subplot(2, 1, 2)
plt.scatter(x,signal)
plt.plot(x,y,'r')
plt.plot(x, fitfunc(p1,x),'g')
plt.show()