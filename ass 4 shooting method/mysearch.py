# Standard Python Packages
import sys
import numpy as np
import matplotlib.pyplot as plt
from math import *


"""
Function Name: secant

Description: Finds the root of a given function i.e. f

Inputs:
- f: function
- x0, x1: estimates where the root is bounded i.e. [x0, x1]
- tol: how off a point is allowed to be from x-axis

Outputs:
- x: estimation of root

Author: Tegh 
"""
def secant(f, x0, x1, tol):

    F0 = f(x0)
    F1 = f(x1)
    iteration_counter = 0
    
    while abs(F1) > tol and iteration_counter < 100:
        
        try:
            denominator = (F1 - F0) / (x1 - x0)
            x = x1 - F1 / denominator
        
        except ZeroDivisionError:
            sys.exit(f'Error! - denominator zero for s = {x1}') 
        
        x0 = x1
        x1 = x
        F0 = F1
        F1 = f(x1)
        iteration_counter = iteration_counter + 1
           
    # here, either a solution is found or too many iterations
    if abs(f(x)) > tol:
        sys.exit("Unable to converge to the solution in the given number of iterations")
        
    return x




