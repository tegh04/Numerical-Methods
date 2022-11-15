import numpy as np
import sys
from math import *
import mysearch_v1 as mys


"""
Function Name: calculate_derivative
Author: Tegh Bir Singh

Description: calculate the derivative i.e. grad f by using the forward step derivative.

Inputs:
- f: Function
- X: 1-D vector that is the argument for f  

Outputs: gradient vector
"""
def calculate_derivative(f, X):
    
    h = 1e-7
    gradient_vector = []
    
    for index in range(len(X)):
        
        Xh = X.copy()
        Xh[index] = X[index] + h
        grad_i = (f(Xh)-f(X))/h
        gradient_vector.append(grad_i)
        
    print(f'gradient_vector: {gradient_vector}')
    return np.array(gradient_vector)


"""
Function Name: gradsearch
Author: Tegh Bir Singh

Description: Uses gradient descent to find the minimum of some function.

Inputs: 
-         f: some mathematical function e.g. ???
-         X: 1D input vector for function f 
- tolerance: termination parameter (see note)

Outputs: estimated minimum of the function.

Technical Notes: 

-In place of computing the actual gradient vector the partial 
forward step derivative was used instead.

-We use golden section search i.e. mys.gss, to determine the stepsize.
"""
def gradsearch(f, X, tolerance):
    
    """
    We use the fline(t) function to convert a 3D function to a 1D function.
    That way it can optimised using golden search in the mysearch_v1 module
    
    Inputs:
    
    t: scalar variable, can be thought of as stepsize in a cerain direction
    """
    ##### Required functions #####
    def fline(t):
        nonlocal X
        nonlocal direction
        Xn = X + t*direction 
        
        fline = f(Xn)
        return fline
    ##############################
    
    
    # Gradient descent variables
    f_values = np.array([f(X)])
    X_list = np.array(X, ndmin = 2)
    X = np.array(X, dtype=float)
    num_iterations = 0
    
    while True:
        
        # Calculate direction
        gradient_vector = calculate_derivative(f, X)
        l2_norm = sqrt((gradient_vector**2).sum())
        direction = -gradient_vector/l2_norm
        
        
        # Find the stepsize i.e. tmin using gss i.e. golden section search
        tmin = mys.gss(fline, 0, 5, 1e-5)
        
        # Apply Gradient Descent equation 
        X = X + tmin*direction
        
            
        # Update X_list and f_values
        X_list = np.append(X_list, [X], axis=0)
        f_values = np.append(f_values, f(X))
        
        
        ##### Perform checks to determine whether to terminate the program #####
        
        # check if change < tolerance
        if len(X_list)>4:  # Due to the initial point more than 4 points are needed
            
            """ We will see how much change has occur current iteration vs the previous 3 
            iterations. i.e. by finding the distance between X and the points X_list[-4:-1]. 
            Then taking sum of thoses distances"""
        
            change = np.array([sqrt(((change-X)**2).sum()) for change in X_list[-4:-1]])
            change = change.sum()
            
            if (change < tolerance):
                
                print('change < tolerance. Returning best value of X')
                index = np.argmin(f_values)
                return X_list[index]
        
        # check if convergence is too slow
        num_iterations = num_iterations + 1
        if num_iterations == 1000:
            sys.exit('The algorithm has not converged for given number of iterations. Ending the program.')
        
        # check if the magnitude of the gradient vector is too small
        if sqrt((gradient_vector**2).sum()) < 0.0001:
            print("The gradient vector is too small. Returning the best value of X")
            index = np.argmin(f_values)
            
            return X_list[index]

