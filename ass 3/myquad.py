import numpy as np
from scipy.stats import qmc


"""
Author:Tegh

Function: monty3d(f, region, a, b, num_points)

Description: Perform Monte Carlo integration in a 3 dimensional space.  

Inputs:
             f: function that takes as input x, y, z
        region: A function that returns a 1D boolean array. This allows a certain region
                to be integrated but not others.
          a, b: Some constants. This is the bound for integration for x, y, z e.g. a=0, b=5
    num_points: Controls the number of points to be used for integration. the numbe points used 
                is given by 2**num_points   

"""

def monty3d(f, region, a, b, num_points):
    
    # generate a column of quasi-random numbers for each variable (x, y, z)
    sampler = qmc.Sobol(d=3)
    p = sampler.random_base2(num_points)
    
    # Since p~U(0, 1) we need to apply a transformation to get variables in bounds [a, b]
    x = p[:, 0]*(b-a) + a
    y = p[:, 1]*(b-a) + a
    z = p[:, 2]*(b-a) + a
    
    # Get points which satify criteria set in the user-defined region function
    bool_mask = region(x, y, z)

    # Get the coordinates of x, y, z which satisfy the criteria using bool_mask
    x = p[bool_mask, 0]
    y = p[bool_mask, 1]
    z = p[bool_mask, 2]
    
    # Evaluate the function at each point
    f_values = np.array([f(xi, yi, zi) for xi, yi, zi in zip(x, y, z)])
    
    # Get the mean value across the region [a, b] then multiply by the volume i.e. (b-a)**3
    return f_values.sum()/len(p[:,0]) * (b-a)**3