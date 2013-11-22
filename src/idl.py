'''
Created on Jan 28, 2013

@author: jwe

See /usr/local/itt/idl71/lib for details


'''
def asin(x):
    from numpy import arcsin
    return arcsin(x)

def shift(x, d):
    from numpy import roll
    return roll(x, d)

def deriv(x, y = []):
    """Perform numerical differentiation using 3-point, Lagrangian
       interpolation."""
    from numpy import roll
    n = len(x)
    if n < 3:
        print 'Parameters must have at least 3 points'
        return
    if len(y) > 0 and len(y) != len(x):
        print 'Vectors must have same size'
        return
    if len(y) == 0:
        d = (roll(x,-1) - roll(x,1))/2.
        d[0] = (-3.0*x[0] + 4.0*x[1] - x[2])/2.
        d[n-1] = (3.*x[n-1] - 4.*x[n-2] + x[n-3])/2.
        return d   
    x12 = x - roll(x,-1)   # x1 - x2
    x01 = roll(x,1) - x    # x0 - x1
    x02 = roll(x,1) - roll(x,-1) # x0 - x2
    
    # Middle points
    d = roll(y,1) * (x12 / (x01*x02)) + \
        y * (1./x12 - 1./x01) - \
        roll(y,-1) * (x01 / (x02 * x12))

    # Formulae for the first and last points:
    # First point
    d[0] = y[0] * (x01[1]+x02[1])/(x01[1]*x02[1]) - \
    y[1] * x02[1]/(x01[1]*x12[1]) + \
    y[2] * x01[1]/(x02[1]*x12[1])
    
    n2 = n-2
    # Last point
    d[n-1] = -y[n-3] * x12[n2]/(x01[n2]*x02[n2]) + \
    y[n-2] * x02[n2]/(x01[n2]*x12[n2]) - \
    y[n-1] * (x02[n2]+x12[n2]) / (x02[n2]*x12[n2])
    
    return d

def n_elements(x):
    return len(x)

def int_tabulated(x, f, sort=False):
    """This function integrates a tabulated set of data { x(i) , f(i) },
       on the closed interval [min(X) , max(X)]."""
    from numpy import sum, arange
    from scipy.interpolate import interp1d
    
    xsegments = len(x) - 1
    
    if sort:
        ii = x.sorted()
        x = x[ii]
        f = f[ii]
    
    # Uniform step size.
    while (xsegments % 4) != 0L:
        xsegments += 1

    xmin = min(x)
    xmax = max(x)

    h = (xmax+0.0 - xmin) / xsegments
    
    # Compute the interpolates at Xgrid.
    # x values of interpolates >> Xgrid = h * FINDGEN(Xsegments + 1L) + Xmin
    fi = interp1d(x, f, kind='cubic')
    x1 = h * arange(xsegments + 1) + xmin
    z = fi(x1)
    
    # Compute the integral using the 5-point Newton-Cotes formula.
    ii = (arange((len(z) - 1)/4)+1) * 4
    
    return sum(2. * h * (7. * (z[ii-4] + z[ii]) + \
            32 * (z[ii-3] + z[ii-1]) + 12. * z[ii-2]) / 45.)
    