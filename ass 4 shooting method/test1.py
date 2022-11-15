#!/usr/bin/env python
# coding: utf-8

# In[3]:


from math import *
from myodes import *
import matplotlib.pyplot as plt


# In[ ]:





# In[2]:


# Test case 1
def d2ydx2(x, y, dydx):
    return -y

X, Y = bvpsolve(d2ydx2 = d2ydx2, L = pi/2, a=0, b = 2, tol=0.01)


# In[ ]:





# In[5]:


# Plotting of test case 1
plt.plot(X, Y[:, 0])


# In[ ]:




