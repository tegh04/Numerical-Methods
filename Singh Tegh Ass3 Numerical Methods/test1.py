import numpy
from myquad import monty3d


"""
Test case 1 (Volume of unit sphere)

Answer: 4/3 *pi*(1^3) = 4.19 (approx)

"""

def region(x,y,z):
    arr = x**2+y**2+z**2
    bool_mask = arr <1
    
    return bool_mask


def f(x,y,z):
    f = 1 
    return f


a,b = -1,1
M = 20
volume_of_sphere = monty3d(f,region,a,b,M)

print(f"The volume of a unit sphere is approximately: {volume_of_sphere}")


