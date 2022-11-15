import numpy
from myquad import monty3d

"""
Test case 2 (Volume of unit cube)

Answer: 1^3 = 1 
"""

def region(x,y,z):
    
    mask1 = x < 1
    mask2 = y < 1
    mask3 = z < 1
    
    bool_mask = mask1 & mask2 & mask3
    
    return bool_mask

def f(x,y,z):
    f = 1 
    return f

# Test case 1
a,b = 0,2
M = 20
unit_cube = monty3d(f,region,a,b,M)

print(f"The volume of a unit cube is: {unit_cube}")
