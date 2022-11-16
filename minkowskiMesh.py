# -*- coding: utf-8 -*-
"""
Created on Thu Sep 29 15:53:20 2022

@author: jonijura
same problem as in wave1d but now solved directly on a 1+1 dimensional spacetime mesh
Using minkowski metric the wave equation d_tt f = d_xx f is equal to
*d*d f = 0
initial values are given to f so that f = [u,i] vector, with values located on primal
vertices, i corresponds to known initial values at t=0 and u to unknown.
boundary values f = 0 could not be encorporated to initial values, but in the hodge h0i
mapping dual 2 forms to primal 0 forms, such that h0i=0 for boundary dual 2 forms
The circumcenters are calculated in the sense of minkowski metric by modifying circumcenter.py
h1 is modified by multiplying timelike primal edge values with -1
    timelike if (dt)**2-(dx)**2>0
Representing the matrix *d*d as [A1,A2] *d*d f=0 leads to
    [A1,A2]*[u,i].T=[A1u+A2i]=0,
where u is solved as the least squares solution to
    A1u = -A2i = b

why do boundary values need to be in h0i?
why is the solution right only it bt = n*pi?
    implicit von neumann boundary condition at the end

requires overriding circumcenter in circumcenter.py with
    if len(pts)==3:
        pts = asarray(pts)
        A = pts[1:,:]-pts[0,:]
        A[:,1]=-A[:,1]
        b = (A[:,0]*A[:,0]-A[:,1]*A[:,1])/2
        x = solve(A,b) + pts[0,:]
        dst = sqrt(abs((x[1]-pts[0][1])**2 - (x[0]-pts[0][0])**2))
        return(x,dst)
"""

import numpy as np
import matplotlib.pyplot as plt
from gmesh_recmesh import makerec,reqrecMsh
from pydec import simplicial_complex
from scipy.sparse import diags

Nx = 10
bx = 2*np.pi
bt = 2*np.pi
dx = bx/Nx

circ = 15

V,E = makerec(bx,bt,dx) #makerec
V = np.loadtxt("C:\MyTemp\cpp\Samples\MeshGeneration2\\build\\vert.txt")*2*np.pi
E = np.loadtxt("C:\MyTemp\cpp\Samples\MeshGeneration2\\build\\tria.txt",dtype='int32')
# E = np.loadtxt("C:\MyTemp\python\mod.txt",dtype='int32')
V[:,[1,0]]=V[:,[0,1]]
sc = simplicial_complex(V,E)

plt.figure(figsize=(8, 8), dpi=80)
plt.triplot(V[:,0], V[:,1], E)
plt.plot(sc[2].circumcenter[:,0],sc[2].circumcenter[:,1],'.g')
# r = np.linalg.norm(V[sc[0].simplices[sc[2].simplices[circ][1]][0]]-sc[2].circumcenter[circ,:])
# circle2 = plt.Circle((sc[2].circumcenter[circ,0],sc[2].circumcenter[circ,1]), r, color='b', fill=False)
# ax = plt.gca()
# ax.add_patch(circle2)
edgeVectors = V[sc[1].simplices[:,0]]-V[sc[1].simplices[:,1]]
timelikeEdges = 2*(edgeVectors[:,1]**2-edgeVectors[:,0]**2>0)-1
drawtimelike = np.where((edgeVectors[:,1]**2-edgeVectors[:,0]**2>0))
plt.plot(sc[1].circumcenter[drawtimelike,0],sc[1].circumcenter[drawtimelike,1],'.r')
plt.show()

xmin = np.min(V[:,0])
xmax = np.max(V[:,0])
tmin = np.min(V[:,1])
tmax = np.max(V[:,1])

eps = 1e-10
tbord = np.where((abs(V[:,0]-xmin)<eps) | (abs(V[:,0]-xmax)<eps))
# plt.plot(V[tbord,0],V[tbord,1],'.')
xinit = np.where(abs(V[:,1]-tmin)<eps)
xend = np.where(abs(V[:,1]-tmax)<eps)
init = np.where((abs(V[:,0]-xmin)<eps) | (abs(V[:,0]-xmax)<eps) | (abs(V[:,1]-tmin)<eps)) #np.unique(np.concatenate((tbord,xinit),1))
initC = np.where( ~ ((abs(V[:,0]-xmin)<eps) | (abs(V[:,0]-xmax)<eps) | (abs(V[:,1]-tmin)<eps)))

h1 = sc[1].star.diagonal()
h1 = diags(h1*timelikeEdges)
h0i = sc[0].star_inv.diagonal()
h0i[tbord]=0
h0i = diags(h0i)

lap = h0i * sc[0].d.T * h1 * sc[0].d
print(f"h1  min: {np.min(h1.diagonal())} \r\n\tmax: {np.max(h1.diagonal())}")
print(f"h0i min: {np.min(h0i.diagonal())} \r\n\tmax: {np.max(h0i.diagonal())}")
print(f"lap min: {np.min(lap)} \r\n\tmax: {np.max(lap)} ")

A2 = lap.T[init].T
A1 = lap.T[initC].T
f = np.zeros(len(V))
f[xinit] = np.sin(V[xinit,0])
# f[xend] = np.sin(V[xend,0])*np.cos(tmax);

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
plt.title("exact solution")
ax.plot_trisurf(V[:,0], V[:,1], exact,cmap='viridis', edgecolor='none')
plt.show()

print("exact solution error: "+str(np.linalg.norm(lap*exact)))

# =============================================================================
# ax = plt.axes()
# ax.tricontour(V[:,0], V[:,1], abs(exact))
# =============================================================================

# =============================================================================
# exact = np.sin(V[:,0])*np.cos(V[:,1])
# ax = plt.axes()
# plt.title("d*exact")
# ax.tricontour(sc[1].circumcenter[:,0], sc[1].circumcenter[:,1], abs(sc[0].d*exact),levels=[0.1])
# plt.show()
# 
# exact = np.sin(V[:,0])*np.cos(V[:,1])
# ax = plt.axes()
# plt.title("h1*d*exact")
# ax.tricontour(sc[1].circumcenter[:,0], sc[1].circumcenter[:,1], abs(h1*sc[0].d*exact),levels=[0.1])
# plt.show()
# 
# exact = np.sin(V[:,0])*np.cos(V[:,1])
# ax = plt.axes(projection='3d')
# ax.set_zlim( [-1,1])
# plt.title("d.T*h1*d*exact")
# ax.plot_trisurf(V[:,0], V[:,1], sc[0].d.T * h1 * sc[0].d*exact,cmap='viridis', edgecolor='none')
# plt.show()
# =============================================================================

exact = np.sin(V[:,0])*np.cos(V[:,1])

fig = plt.figure(figsize=(20,14), dpi=80)
ax = fig.add_subplot(231)
plt.title("dxx-dtt of exact solution")
ax.tricontour(V[:,0], V[:,1], abs(lap*exact),levels = [0.0006,0.002,0.006])#
ax = fig.add_subplot(232)
plt.triplot(V[:,0], V[:,1], E)
ax = fig.add_subplot(233,projection='3d')
ax.plot_trisurf(V[:,0], V[:,1], lap*exact,cmap='viridis', edgecolor='none')

ax = fig.add_subplot(234)
plt.title("error")
ax.tricontour(V[:,0], V[:,1], abs(f-exact),levels = [0.00006,0.0002,0.0006])#,levels = [0.0001,0.0003,0.001]
ax = fig.add_subplot(235)
plt.triplot(V[:,0], V[:,1], E)
ax = fig.add_subplot(236,projection='3d')
ax.plot_trisurf(V[:,0], V[:,1], f-exact,cmap='viridis', edgecolor='none')
plt.show()

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

# =============================================================================
# fig = plt.figure(figsize=(20,8), dpi=80)
# ax = fig.add_subplot(231)
# plt.title("1d")
# ax.plot(lp1, ex1)
# ax.plot(lp1, sim1)
# ax = fig.add_subplot(234)
# plt.title("error")
# ax.plot(lp1, ex1-sim1)
# 
# ax = fig.add_subplot(232)
# plt.title("2d")
# ax.plot(lp2, ex2)
# ax.plot(lp2, sim2)
# ax = fig.add_subplot(235)
# plt.title("error")
# ax.plot(lp2, ex2-sim2)
# 
# ax = fig.add_subplot(233)
# plt.title("3d")
# ax.plot(lp3, ex3)
# ax.plot(lp3, sim3)
# ax = fig.add_subplot(236)
# plt.title("error")
# ax.plot(lp3, ex3-sim3)
# plt.show()
# =============================================================================
