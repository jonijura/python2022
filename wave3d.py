# -*- coding: utf-8 -*-
"""
Created on Tue Sep 27 13:07:22 2022

@author: jonijura
"""

import numpy as np
from scipy.sparse import diags
import matplotlib.pyplot as plt
from gmesh_recmesh import makecube
from pydec import simplicial_complex
import time

Nxyz = 10
bxyz = np.pi
dxyz = bxyz/Nxyz

bt = 6
Nt = 400
dt = bt/Nt

prev = time.time()
print("making mesh")
V,E = makecube(bxyz,bxyz,bxyz,dxyz)
print("\tmesh complete in "+str(time.time()-prev))
prev = time.time()
print("making complex")
sc = simplicial_complex(V,E)
print("\tcomplex complete in "+str(time.time()-prev))

print("calculating operations")
prev = time.time()
d0=sc[0].d
print("\td0 calculated in "+str(time.time()-prev))
prev=time.time()
h0i=sc[0].star_inv
print("\th0i calculated in "+str(time.time()-prev))
prev=time.time()
h1=sc[1].star
print("\th1 calculated in "+str(time.time()-prev))
#some mesh quality statistics
print("h1 min: "+str(np.min(h1.diagonal())))
print("h1 max: "+str(np.max(h1.diagonal())))
print("h0i min: "+str(np.min(h0i.diagonal())))
print("h0i max: "+str(np.max(h0i.diagonal())))

#set boundary condition f=0
prev = time.time()
print("finding boundary nodes")
boundFcs = sc.boundary()
boundEdgs = set([boundEdg for face in boundFcs for boundEdg in face.boundary()])
boundVerts = set([boundVert for edge in boundEdgs for boundVert in edge.boundary()])
boundVertInds = [sc[0].simplex_to_index[v] for v in boundVerts]
print("found boundary nodes in "+str(time.time()-prev))

h0ivals = h0i.diagonal()
h0ivals[boundVertInds]=0
h0i = diags(h0ivals)

g = np.zeros(sc[1].num_simplices)
f = np.zeros(sc[0].num_simplices)

A = dt*h1*d0
B = dt*h0i*(-d0.T)

#initial values f(0), g(dt/2)
f = np.sin(V[:,0])*np.sin(2*V[:,1])*np.sin(V[:,2])
g = 0.5*A*f

refp = len(V)-200

hist = np.zeros(Nt+1)
hist[0]=f[refp]

prev=time.time()
print("timestepping")
for i in range(Nt):
    f+=B*g
    g+=A*f
    hist[i+1]=f[refp]
print("timestepping complete in "+str(time.time()-prev))

lp = np.linspace(0,bt, Nt+1)
plt.title("History of reference point")
plt.plot(lp, hist, 'r')
plt.plot(lp, np.sin(V[refp,0])*np.sin(2*V[refp,1])*np.sin(V[refp,2])*np.cos(np.sqrt(5)*lp),'b')
plt.legend(["simulation", "analytic"],loc=2)
plt.show()













