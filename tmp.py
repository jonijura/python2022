# -*- coding: utf-8 -*-
"""
Created on Mon Nov  7 15:19:44 2022

@author: jonijura
"""

import numpy as np

x = 4
y = 2
t = 1

A = np.array([[-x,-x,x*t],[y,0,y*t],[0,t,-x*t]])

Ai = np.linalg.inv(A)

b = np.matmul(Ai,np.array([1,0,0]))
b2 = np.matmul(Ai,np.array([0,1,0]))
c = np.cross(b,b2)

# d11 = np.matmul(A,np.array([0,1,0]))
# d12 = np.matmul(A,np.array([0,0,1]))
# d1 = c[0]*np.cross(d11,d12)

# d21 = np.matmul(A,np.array([0,0,1]))
# d22 = np.matmul(A,np.array([1,0,0]))
# d2 = c[1]*np.cross(d21,d22)

# d31 = np.matmul(A,np.array([1,0,0]))
# d32 = np.matmul(A,np.array([0,1,0]))
# d3 = c[2]*np.cross(d31,d32)

# print(d1+d2+d3)

print(c[2]/2)
print(x*y/2)

# ans = c[2]*np.linalg.det(A)

# print(ans/2)

# we = (x*t+y*x+y*t)*c[2]/2
# print(we)

# area = x*y/np.sqrt(x*x*y*y+x*x*t*t+y*y*t*t)
# print(x*y/2)
