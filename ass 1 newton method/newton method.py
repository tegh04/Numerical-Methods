# Basic Python packages
import numpy as np
import sys


"""
Function Name: newtonzero
Author: Tegh Bir Singh

Description: Use quasi-Newton's method to calculate the roots of some function. 
In place of the actual derivative the symmetric derivative was used instead.

Inputs: 
-         f: some mathematical function e.g. x^3-x+2
-        x0: initial point 
- tolerance: maximum error that is acceptable

Outputs: estimate of the root of the function

Note: A division by zero error can occur if dfdx = 0
"""

def newtonzero(f, x0, tolerance):

    ### Variables ###
    
    h = 10**-6   # step-size

    # Flags to determine if program needs terminate
    
    num_iterations = 0 # Store number of time Newton's method runs
    error = abs(f(x0)) # how far y-value away from zero
    xvalues = np.array([x0]) # store xvalues to determine looping
    
    
    ##### Beginning of Newton's Method #####
    while(error > tolerance):

        # calculate the derivative
        dfdx = (f(x0+h)-f(x0-h))/(2*h)

        # use the equation for Newton's method
        x0 = x0 - f(x0)/dfdx  


        ####### Check for issues with Newton's method #######
        
        if x0 in xvalues:
            sys.exit("Looping has occured. Terminating program")
    
        
        num_iterations = num_iterations + 1
        if num_iterations == 50:
            sys.exit("Convergence is too slow. Terminating program.")
            
        # check for running off into +/-infinity
        if num_iterations>=6:
            
            # if the new-value is greater than or less than the all the last 5 values
            # from Newton's Method then we assume it is heading to +/-infinity and hence return None
            # Note: x0 is not in xvalues
            if (x0 > xvalues[-5:]).sum()==5:
                sys.exit("The program seems to be running to infinity. Terminating program.")
                
            elif (x0 < xvalues[-5:]).sum()==5:
                sys.exit("The program seems to be running to -infinity. Terminating program.")
                
        
        ####### All checks completed #######
        
        
        # Update variables
        xvalues = np.append(xvalues, x0)
        error = abs(f(x0))
            
    
    return x0

