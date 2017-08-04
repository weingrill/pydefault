'''
Created on Nov 27, 2013

@author: jwe
'''
def tau(bv, gilliland = True):
    """
    Table 1 from Barnes & Kim 2010 p.678:
    
    Global \tau_C (d)
    log Teff to B-V conversion by Mamajek
    grid extended with values from Gilliland
    """
    import numpy as np
    from scipy.interpolate import interp1d
    a = np.array([[1.640, 339.80000], [1.618, 349.08619], [1.595, 358.37237], 
                  [1.573, 367.65856], [1.550, 387.12493], [1.538, 406.86305], 
                  [1.529, 454.74676], [1.520, 470.89879], [1.501, 380.32005], 
                  [1.495, 297.41098], [1.489, 220.68335], [1.483, 198.54819], 
                  [1.475, 176.62500], [1.461, 164.12500], [1.447, 151.62500], 
                  [1.432, 142.09721], [1.416, 134.43794], [1.398, 126.77867], 
                  [1.378, 120.56730], [1.350, 114.67376], [1.323, 108.78023], 
                  [1.298, 104.14129], [1.274, 100.13722], [1.249, 96.13314], 
                  [1.224, 92.29436], [1.200, 89.00896], [1.175, 85.72356], 
                  [1.151, 82.43816], [1.131, 79.52713], [1.121, 76.85144], 
                  [1.111, 74.17575], [1.101, 71.50006], [1.075, 69.31314], 
                  [1.047, 67.27058], [1.020, 65.22801], [0.992, 63.46932], 
                  [0.969, 61.58586], [0.946, 59.74624], [0.923, 57.93177], 
                  [0.899, 56.11865], [0.869, 54.27821], [0.842, 52.50865], 
                  [0.827, 50.77602], [0.803, 49.02645], [0.766, 47.11458], 
                  [0.735, 44.91198], [0.722, 42.67782], [0.709, 40.34714], 
                  [0.694, 38.09877], [0.671, 35.79291], [0.654, 33.10995], 
                  [0.629, 30.34246], [0.601, 27.43598], [0.578, 24.41565], 
                  [0.558, 20.94642], [0.542, 17.18560], [0.528, 13.60996], 
                  [0.511, 10.06345], [0.493, 6.62785], [0.474, 4.37202], 
                  [0.455, 2.71518], [0.436, 1.67648], [0.420, 0.96018], 
                  [0.400, 0.56288], [0.382, 0.29435], [0.369, 0.15539], 
                  [0.355, 0.08260], [0.341, 0.04843], [0.324, 0.04190], 
                  [0.304, 0.05169], [0.287, 0.00717], [0.274, 0.00477], 
                  [0.261, 0.00432], [0.249, 0.00433], [0.238, 0.00547], 
                  [0.226, 0.01860], [0.215, 0.05191], [0.201, 0.10048], 
                  [0.182, 0.15311], [0.166, 0.22745], [0.155, 0.30827], 
                  [0.145, 0.42672], [0.131, 0.58051], [0.115, 0.79248], 
                  [0.098, 1.07463], [0.087, 1.46778], [0.081, 2.03562], 
                  [0.076, 2.77835], [0.068, 3.80599], [0.060, 5.01401], 
                  [0.052, 6.65116], [0.044, 9.07035], [0.034, 12.37053], 
                  [0.024, 17.18070], [0.014, 23.52624], [0.005, 30.60496]])
    
    tau_int = interp1d(a[:,0], a[:,1], kind='slinear',assume_sorted=False )
    try:
        return tau_int(bv)
    except ValueError:
        return np.nan


def gyroage(bv, P, P0=1.1):
    """
    calculates t from P derived from Barnes 2010 
    result in Myr
    """
    import numpy as np
    
    k_C = 0.646 # days Myr^-1 Barnes 2010 p.224
    k_I = 452.0 # Myr day^-1 Barnes 2010 p.224
    if P>P0:
        t = (tau(bv)/k_C)*np.log(P/P0) + (0.5*k_I/tau(bv))*(P**2 - P0**2) 
        return t
    else:
        return np.NaN

def gyroperiod(bv, age, P0=1.1, version=2010):
    """
    calculates the expected period according to the given age and color
    """
    from scipy.optimize import minimize
    import numpy as np

    if version == 2010:
        P = np.empty(len(bv))
        
        def min_func(x, bv, age, P0):
            return abs(gyroage(bv, x, P0) - age)
    
        x0 = 10.0
        for i in range(len(bv)):
            res = minimize(min_func, x0, method='nelder-mead', 
                           args=(bv[i], age, P0),
                           options={'xtol': 1e-3, 'disp': False})
            P[i] = res['x'][0]
        P[P<P0] = P0
        return P

    if version == 2007:
        """periods as a function of age and B-V from Barnes 2007"""
        P = 0.7725*np.power(bv-0.4,0.601)*np.power(age, 0.5189)
        
        return P
    
    if version == 2003:
        #eqns. 1 and 2 from Barnes 2003
        pi = np.sqrt(age)*(np.sqrt(bv - 0.5) - 0.15*(bv - 0.5) )
        pi[0] = 0.1
        #eqn. 15 from Barnes 2003
        pc=0.2*np.exp((age/100.0)/np.power(bv + 0.1 - (1.0/3000.)*age, 3.))
        pc[pi<pc] = pi
        
        return pi, pc
        


if __name__ == '__main__':
    #import matplotlib
    #matplotlib.use('WXAgg')
    import numpy as np
    from scipy.optimize import minimize
    import matplotlib.pyplot as plt
    
    bv = np.linspace(0.0, 1.6, num=100)
    #bv = np.linspace(0.473, 1.631, num=100)
    P = np.empty(len(bv))
    
    age = 100
    
    P11 = gyroperiod(bv, age)
    P34 = gyroperiod(bv, age, P0=3.4)
    P01 = gyroperiod(bv, age, P0=0.1)
    
    for k in zip(bv,P01,P11,P34):
        print '%.2f\t%.3f\t%.3f\t%.3f' % k
    
    from breakup import breakup
    bv_b, P_b = breakup()
    
    plt.plot(bv, P11, 'g--')
    plt.plot(bv, P01, 'b--')
    
    plt.plot(bv, P34, 'r')
    plt.plot(bv_b,P_b, 'k')
    plt.xlabel('B - V')
    plt.ylabel('P rot')
    plt.title('%d Myr' % age)
    plt.grid()
    plt.show()
    exit()
    def min_func(x, bv, age):
        return abs(gyroage(bv, x) - age)
    x0 = 10.0

    for age in [500.,800.,1000.,2000.,4000.,4570.]:
        print age,' Myr'    
        for i in range(len(bv)):
            res = minimize(min_func, x0, method='nelder-mead', args=(bv[i], age),
                       options={'xtol': 1e-3, 'disp': False})
            P[i] = res['x'][0]
            #print '%.3f %6.3f' % (bv[i], P[i])
        #plt.scatter(bv, P)
        plt.plot(bv, P, 'r')
    plt.scatter(0.65,26.09)
    plt.xlabel('B-V')
    plt.ylabel('period (days)')
    plt.grid()
    plt.show()
