# -*- coding: utf-8 -*-
'''
Created on Jan 15, 2015

@author: jwe
'''
def poly(x, c):
    """
    The POLY function evaluates a polynomial function of a variable. 
    © 2015 Exelis Visual Information Solutions
    
    The result is equal to:
    C[0] + C[1]*x + C[2]*x**2 + ...
    """
    from numpy import polyval
    return polyval(c[::-1], x)