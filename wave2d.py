# -*- coding: utf-8 -*-
"""
Created on Mon Sep 26 09:31:52 2022

@author: jonijura

vibrating rectangular membrane
f = string displacement, primal vertices
lap f = f'' 
-> system
    d_t f = *dg
    d_t g=*df
leapfrog time discretization

boundary condition
    f=0, encorporated to h0i
initial condition
    f(x,y,0)=sin(x)sin(2y)
    g=0
solution
    f(x,t)=cos(sqrt(5)t)sin(x)sin(2y)
"""

import numpy as np
from scipy.sparse import diags
import matplotlib.pyplot as plt
from gmesh_recmesh import makerec
from pydec import simplicial_complex

Nxy=10
bxy = np.pi
dxy = bxy/Nxy

bt = 5
Nt = 50
dt = bt/Nt

V,E = makerec(bxy,bxy,dxy)
sc = simplicial_complex(V,E)

d0=sc[0].d
h1=sc[1].star
h0i=sc[0].star_inv
#some mesh quality statistics
print("h1 min: "+str(np.min(h1.diagonal())))
print("h1 max: "+str(np.max(h1.diagonal())))
print("h0i min: "+str(np.min(h0i.diagonal())))
print("h0i max: "+str(np.max(h0i.diagonal())))
#set boundary condition f=0
boundEdgs = sc.boundary()
boundVerts = set([boundVert for edge in boundEdgs for boundVert in edge.boundary()])
boundVertInds = [sc[0].simplex_to_index[v] for v in boundVerts]
h0ivals = h0i.diagonal()
h0ivals[boundVertInds]=0
h0i = diags(h0ivals)

g = np.zeros(sc[1].num_simplices)
f = np.zeros(sc[0].num_simplices)

A = dt*h1*d0
B = dt*h0i*(-d0.T)

#initial values f(0), g(dt/2)
f = np.sin(V[:,0])*np.sin(2*V[:,1])
g = 0.5*A*f

for i in range(Nt):
    f+=B*g
    g+=A*f

ax = plt.axes(projection='3d')
ax.set_zlim( [-1,1])
plt.title("Result")
ax.plot_trisurf(V[:,0], V[:,1], f,
                cmap='viridis', edgecolor='none')
plt.show()
    













