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

def tau_gill(bv):
    """
    taken from Gilliland 1985
    """
    from scipy.interpolate import interp1d
    BV = [-0.044, -0.035, -0.022, -0.013, -0.011, -0.001, 0.005, 0.015, 0.027, 
          0.034, 0.042, 0.052, 0.054, 0.065, 0.067, 0.075, 0.078, 0.088, 0.091, 
          0.102, 0.107, 0.115, 0.121, 0.126, 0.132, 0.137, 0.137, 0.143, 0.149, 
          0.152, 0.157, 0.157, 0.163, 0.166, 0.171, 0.172, 0.174, 0.177, 0.180, 
          0.180, 0.183, 0.183, 0.185, 0.186, 0.191, 0.194, 0.209, 0.234, 0.256, 
          0.272, 0.285, 0.288, 0.289, 0.292, 0.292, 0.292, 0.295, 0.296, 0.302, 
          0.312, 0.318, 0.322, 0.329, 0.332, 0.339, 0.346, 0.349, 0.356, 0.359, 
          0.366, 0.370, 0.377, 0.380, 0.387, 0.390, 0.397, 0.401, 0.408, 0.415, 
          0.422, 0.429, 0.436, 0.447, 0.454, 0.465, 0.472, 0.479, 0.490, 0.501, 
          0.515]
    tau_C = [25.3569, 21.1591, 13.3237, 11.1177, 9.6580, 8.0589, 6.2051, 
             4.5892, 3.4631, 2.8898, 2.3633, 1.9721, 1.7132, 1.4295, 1.2418, 
             1.0362, 0.9002, 0.7511, 0.6395, 0.5336, 0.4364, 0.3642, 0.2804, 
             0.2340, 0.1992, 0.1415, 0.1629, 0.1181, 0.0856, 0.0966, 0.0584, 
             0.0714, 0.0415, 0.0478, 0.0283, 0.0346, 0.0205, 0.0152, 0.0071, 
             0.0115, 0.0053, 0.0088, 0.0042, 0.0045, 0.0037, 0.0033, 0.0034, 
             0.0039, 0.0045, 0.0053, 0.0067, 0.0092, 0.0124, 0.0158, 0.0182, 
             0.0301, 0.0222, 0.0272, 0.0272, 0.0272, 0.0333, 0.0407, 0.0497, 
             0.0608, 0.0744, 0.0909, 0.1112, 0.1332, 0.1629, 0.1992, 0.2436, 
             0.2978, 0.3642, 0.4453, 0.5445, 0.6658, 0.7978, 0.9755, 1.1929, 
             1.4586, 1.7834, 2.1373, 2.6133, 3.1955, 3.8294, 4.6825, 5.7254, 
             6.8615, 8.3897, 10.0542]
    tau_int = interp1d(BV, tau_C, kind='linear')
    return tau_int(bv)


def tau(bv, gilliland = True):
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
    
    if bv<=0.515 and gilliland:
        return tau_gill(bv)
    if bv<=0.46 and not gilliland:
        return 0.001
    
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
    
    bv = np.linspace(0.0, 1.6, num=50)
    #bv = np.linspace(0.473, 1.631, num=100)
    P = np.empty(len(bv))
    
    P11 = gyroperiod(bv, 400)
    P34 = gyroperiod(bv, 400, P0=3.4)
    P01 = gyroperiod(bv, 400, P0=0.1)
    
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
    plt.title('400 Myr')
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
