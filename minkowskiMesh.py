# -*- coding: utf-8 -*-
"""
Created on Thu Sep 29 15:53:20 2022

@author: jonijura
"""

import numpy as np
import matplotlib.pyplot as plt
from gmesh_recmesh import makerec
from pydec import simplicial_complex

Nxy = 10
bxy = 2*np.pi
dxy = bxy/Nxy

circ = 15

V,E = makerec(bxy,bxy,dxy)
sc = simplicial_complex(V,E)

plt.figure(figsize=(8, 8), dpi=80)

plt.triplot(V[:,0], V[:,1], E)
plt.plot(sc[2].circumcenter[:,0],sc[2].circumcenter[:,1],'.')

# r = np.linalg.norm(V[sc[0].simplices[sc[2].simplices[circ][1]][0]]-sc[2].circumcenter[circ,:])
# circle2 = plt.Circle((sc[2].circumcenter[circ,0],sc[2].circumcenter[circ,1]), r, color='b', fill=False)
# ax = plt.gca()
# ax.add_patch(circle2)

plt.show()

lap = sc[0].star_inv * sc[0].d.T * sc[1].star * sc[0].d
print("max: "+ str(np.max(lap)) + "  min: "+str(np.min(lap)))

xmin = np.min(V[:,0])
xmax = np.max(V[:,0])
tmin = np.min(V[:,1])
tmax = np.max(V[:,1])

eps = 1e-10
tbord = np.where((abs(V[:,0]-xmin)<eps) | (abs(V[:,0]-xmax)<eps))
xinit = np.where(abs(V[:,1]-tmin)<eps)
init = np.where((abs(V[:,0]-xmin)<eps) | (abs(V[:,0]-xmax)<eps) | (abs(V[:,1]-tmin)<eps))
initC = np.where( ~ ((abs(V[:,0]-xmin)<eps) | (abs(V[:,0]-xmax)<eps) | (abs(V[:,1]-tmin)<eps)))

A2 = lap.T[init].T
A1 = lap.T[initC].T

f = np.zeros(len(V))
f[xinit] = np.sin(V[xinit,0])

ax = plt.axes(projection='3d')
ax.set_zlim( [-1,1])
plt.title("Initial values")
ax.plot_trisurf(V[:,0], V[:,1], f,cmap='viridis', edgecolor='none')
plt.show()

b = -A2*f[init]
sol = np.linalg.lstsq(A1.toarray(),b,rcond=None)[0]
f[initC] = sol

ax = plt.axes(projection='3d')
ax.set_zlim( [-1,1])
plt.title("Solution")
ax.plot_trisurf(V[:,0], V[:,1], f,cmap='viridis', edgecolor='none')
plt.show()

print("solution error: "+str(np.linalg.norm(lap*f)))

exact = np.sin(V[:,0])*np.cos(V[:,1])
ax = plt.axes(projection='3d')
ax.set_zlim( [-1,1])
plt.title("Exact solution")
ax.plot_trisurf(V[:,0], V[:,1], exact,cmap='viridis', edgecolor='none')
plt.show()

print("exact solution error: "+str(np.linalg.norm(lap*exact)))

# =============================================================================
# pts = [[2.62729316, 1.53518664],[2.10085398, 1.22306732],[2.16674694, 1.70210904]]
# pts = [[-1,0],[1,0],[0,2]]
# pts = np.asarray(pts)
# A = pts[1:,:]-pts[0,:]
# A[:,1]=-A[:,1]
# b = (A[:,0]*A[:,0]-A[:,1]*A[:,1])/2
# x = np.linalg.solve(A,b) + pts[0,:]
# dst = np.sqrt((x[1]-pts[0][1])**2 - (x[0]-pts[0][0])**2)
# 
# def md(a,b):
#     return (a[:,1]-b[1])**2 - (a[:,0]-b[0])**2
# =============================================================================
