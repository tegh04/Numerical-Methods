import numpy as np
from math import *

"""
Minimum bracketing

Inputs:

f: function
(a, b): interval
n: number of points
"""
def minbracket(f, a, b, n):
    
    xvalues = np.linspace(a, b, n)
    
    minimums = []
    h =  xvalues[1] - xvalues[0]
    
    for i in range(1, n):
        x = h*i + xvalues[0] # + a
        if f(x) < f(x-h) and f(x) < f(x+h):
            minimums.append((x-h, x+h))
            
    
    return minimums
    

"""
Golden-section search

Inputs:

f: function
(a, b): interval
tol: desired length of the interval
"""
def gss(f, a, b, tol):
   
    gr = (sqrt(5) + 1) / 2
    
    while abs(b - a) > tol:
        
        # We recompute both c and d here to avoid loss of precision which may lead to incorrect results or infinite loop
        c = b - (b - a) / gr
        d = a + (b - a) / gr

        
        if f(c) < f(d):  # f(c) > f(d) to find the maximum
            b = d
        else:
            a = c

       
    return (b + a) / 2
