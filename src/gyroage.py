'''
Created on Nov 27, 2013

@author: jwe
'''
def tau_vdb(bv):
    import numpy as np
    from scipy.interpolate import interp1d
    bv_tau = np.array([[1.696686, 296.359650],
                       [1.679887,  333.540527],
                       [1.641342,  368.980011],
                       [1.604276,  416.280823],
                       [1.574467,  529.159973],
                       [1.546652,  283.349030],
                       [1.527973,  199.954010],
                       [1.511048,  159.338181],
                       [1.489343,  135.041138],
                       [1.456988,  111.606415],
                       [1.411223,   96.741936],
                       [1.312827,   83.472977],
                       [1.183548,   74.010582],
                       [1.058591,   65.404373],
                       [0.946852,   58.010460],
                       [0.855742,   50.797001],
                       [0.778817,   43.986633],
                       [0.712589,   37.630024],
                       [0.654392,   31.698168],
                       [0.598389,   25.231077],
                       [0.549715,   19.091480],
                       [0.508855,   13.113885],
                       [0.469622,    7.338559],
                       [0.435101,    2.345757]])
    bv_vdb = bv_tau[::-1,0] 
    tau_C = bv_tau[::-1,1]
    tau_int = interp1d(bv_vdb, tau_C, kind='linear')
    return tau_int(bv)

def tau(bv):
    """
    Table 1 from Barnes & Kim 2010 p.678
    """
    import numpy as np
    from scipy.interpolate import interp1d
    B = np.array([15.056,14.000,13.316,12.761,12.308,11.873,11.402,10.895,10.343,9.789,
         9.230,8.696,8.120,7.543,7.040,6.605,6.204,5.835,5.499,5.190,4.898,
         4.623,4.363,4.110,3.873,3.648,3.432,3.227])
    V = np.array([13.425,12.464,11.788,11.247,10.809,10.398,9.949,9.472,8.958,8.455,
         7.959,7.499,7.022,6.558,6.157,5.809,5.480,5.173,4.888,4.620,4.363,
         4.119,3.890,3.667,3.457,3.260,3.076,2.896])
    tau_C = np.array([339.8,367.9,408.6,493.0,358.4,222.1,176.9,146.8,125.4,106.8,93.05,
             81.17,70.89,62.54,55.01,47.90,41.24,34.87,28.64,21.93,14.67,8.141,
             2.394,0.0,0.0,0.0,0.0,0.0])
    BV = (B-V)
    try:
        tau_int = interp1d(BV[::-1], tau_C[::-1], kind='linear')
    except ValueError:
        print bv
    return tau_int(bv)


def gyroage(bv, P, P0=1.1):
    """
    calculates t from P derived from Barnes 2010 
    result in Myr
    """
    import numpy as np
    
    k_C = 0.646 # days Myr^-1 Barnes 2010 p.224
    k_I = 452.0 # Myr day^-1 Barnes 2010 p.224
    
    t = (tau(bv)/k_C)*np.log(P/P0) + (0.5*k_I/tau(bv))*(P**2 - P0**2) 
    return t

def gyroperiod(bv, age, P0=1.1):
    """
    calculates the expected period according to the given age and color
    """
    from scipy.optimize import minimize
    import numpy as np

    P = np.empty(len(bv))
    
    def min_func(x, bv, age, P0):
        return abs(gyroage(bv, x, P0) - age)

    x0 = 10.0
    for i in range(len(bv)):
        res = minimize(min_func, x0, method='nelder-mead', 
                       args=(bv[i], age, P0),
                       options={'xtol': 1e-3, 'disp': False})
        P[i] = res['x'][0]
    
    return P

if __name__ == '__main__':
    import matplotlib
    matplotlib.use('WXAgg')
    import numpy as np
    
    bv = np.linspace(0.5, 1.6, num=10)
    #bv = np.linspace(0.473, 1.631, num=100)
    P = np.empty(len(bv))
    
    def min_func(x, bv, age):
        return abs(gyroage(bv, x) - age)
    x0 = 10.0
    from scipy.optimize import minimize
    import matplotlib.pyplot as plt

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
