'''
Created on 28.04.2013

@author: jwe
'''

def kappa(t,rho):
    a = 5.6707e-8 # radiation constant
    c = 2.99792458e8 # speed of light
    
    return 4.0*a*c*t**3/(3.0*rho*lamb)

def deriv(y,t): # return derivatives of the array y
    a = -2.0
    b = -0.1
    return array([ y[1], a*y[0]+b*y[1] ])

if __name__ == '__main__':
    from scipy.integrate import odeint
    from pylab import linspace,array,figure,plot,xlabel,ylabel,show # for plotting commands
    time = linspace(0.0,10.0,1000)
    yinit = array([0.0005,0.2]) # initial values
    y = odeint(deriv,yinit,time)
    figure()
    plot(time,y[:,0]) # y[:,0] is the first column of y
    xlabel('t')
    ylabel('y')
    show()
    