# Standard Python Packages
import numpy as np
from math import *
from mysearch import secant


"""
Function Name: rk4

Description: A function for solving a system of DE using Runge-Katta 4

Inputs:
=======
system_of_DE : the right side fucntion of the system of ODE
x0, xmax: The x-values to apply Runge-Katta 4 across 
y0 : vector containing y-value (0th index), dydx (1st Index)
d2ydx2: function handle. Takes x, y, dydx as input and gets the second derivative as output

Returns:
========
(X, Y)

X: x-coordinates
Y: vector containing y-value (0th index), dydx (1st Index)

Author: Tegh
"""

def rk4(system_of_DE, x0, xmax, y0, d2ydx2):
    
    N = 100000 # number of steps
    h = xmax/N # step size
    
    X = [x0] # list for storing the grid points
    Y = [y0] # list for storing the numerical solution
    
    # Implementation of RK4 algorithm
    
    # when y0 get passed into f we need get first derivative from y0
    # Get the second derivative from d2ydx2
    for i in range(1, N + 1): 
        
        # Note d2ydx2 is a function, the x/y-coordinate of the function 
        # "system_of_DE" are taken as input in d2ydx2. 
        k1 = system_of_DE(x0, y0, d2ydx2)
        k2 = system_of_DE(x0+h/2, y0+h*k1/2, d2ydx2)
        k3 = system_of_DE(x0+h/2, y0+h*k2/2, d2ydx2)
        k4 = system_of_DE(x0+h, y0+h*k3, d2ydx2)
        y0 = y0 + h * (k1 + 2 * k2 + 2 * k3 + k4)/6
        x0 = x0 + h
        
        X.append(x0)
        Y.append(y0)
        
    return np.array(X), np.array(Y)


"""
Function Name: System of DE

Description: 

Right-handside of the first order system of differential equation.
System of DE are used for solving IVP. 

Inputs:

- x: x-coordinate
- y: contains y-value (0th index), and derivative(1st index) corresponding to x-value

Returns: 
- Right-handside of the first order system of differential equation

Author: Tegh
"""
def system_of_DE(x, y, d2ydx2):
    
    f = np.zeros(2) # Think of f as a vector function
    
    f[0] = y[1] 
    f[1] = d2ydx2(x, y[0], y[1])
    
    return f


"""
Function Name: bvpsolve

Description:
============
Solve the boundary value problem numerically using the shooting method.
Using the runge-katta 4 to solve the IVP from the shooting method and
the secant method to find the best derivative.


Inputs:
=======
d2ydx2 : function outputs the second derivative given x, y, dydx
L : x-value such that f(L) = b, where f(x) is the solution to the differential equation 
a : y-value such that f(0) = a
b: y-value such that f(L) = b
N : number of points uses


Outputs:
========
X: List of x-coordinates
Y: list of coordinates plus the derivative

Author: Tegh
"""

def bvpsolve(d2ydx2, L, a, b, tol):
    
    """
    Description: 
    This is the function we are trying from the root of for the shooting method.
    I.e. This function is used as input to secant method find best value of dydx
    which will be used to construct the solution of BVP differential equation.

    Inputs
    a, b: These represent the y-values at the boundary points
    dydx: derivatice

    Return:
    number indicating how above or below we are from the root 

    """

    def F(dydx):

        nonlocal a, b, L
        
        # initial conditions for IVP that will be solved by RK4
        initial_conditions = np.array([a, dydx])

        X, Y = rk4(system_of_DE, x0=0, xmax=L, y0=initial_conditions, d2ydx2 = d2ydx2)

        F = (Y[-1])[0] - b

        return F

    # We arbitrarily assume the gradient is between -10, 10 when applying shooting method
    # d: for derivative
    d0 = -10
    d1 = 10

    dydx = secant(F, d0, d1, tol)

    initial_conditions = np.array([a, dydx])

    X, Y = rk4(system_of_DE, x0=0, xmax=L, y0=initial_conditions, d2ydx2 = d2ydx2)


    return X, Y



