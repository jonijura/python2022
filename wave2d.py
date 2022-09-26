# -*- coding: utf-8 -*-
"""
Created on Mon Sep 26 09:31:52 2022

@author: jonijura

vibrating rectangular membrane
f = string displacement, dual vertices
lap f = f'' 
-> system
    d_t f = *dg
    d_t g=*df
leapfrog time discretization

boundary condition
    f=0, comes automatically from placing f on dual vertices
initial condition
    f(x,y,0)=sin(x)sin(2y)
    g=0
solution
    f(x,t)=cos(sqrt(5)t)sin(x)sin(2y)
"""

import numpy as np
from scipy.sparse import spdiags, diags
import matplotlib.pyplot as plt

N=10
s=0
e=np.pi
dx = (e-s)/N

ts = 0
te = 20
tN = 100
dt = (te-ts)/tN

lpd = np.linspace(s+dx/2, e-dx/2,N) #locations of dual vertices and circumcenters of primal edges
lp = np.linspace(s, e, N+1) #primal vertices and dual edge circumcenters
d = spdiags([-np.ones(N+1), np.ones(N+1)],[0,1], N,N+1) #discrete derivative
#primal edges to dual vertices
h1 = diags((1/dx)*np.ones(N))
#dual edges to primal vertices, inverse of h0
ih0 = diags(np.concatenate(([2/dx],(1/dx)*np.ones(N-1),[2/dx])))

#time stepping matrices from system
A=dt*h1*d
#discrete derivative for dual vertices: -d.T, deducted from dual edge orientation
B=dt*ih0*-d.T


f = np.sin(lpd)
#g at time dt/2, euler forward
g=0.5*B*f

hist = np.zeros(tN+1)
mp = int(N/2)
hist[0]=f[mp]
for i in range(tN):
    f+=A*g
    g+=B*f
    hist[i+1]=f[mp]

expected = np.cos(te)*np.sin(lpd)

plt.title("lopputulos")
plt.plot(lpd, f,'r')
plt.plot(lpd, expected,'b')
plt.show()

tlp = np.linspace(ts, te, tN+1)
plt.title("keskipisteen historia")
plt.plot(tlp, hist, 'r')
plt.plot(tlp, np.sin(lpd[mp])*np.cos(tlp))
plt.show()
