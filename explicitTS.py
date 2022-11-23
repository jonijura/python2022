# -*- coding: utf-8 -*-
"""
Created on Fri Nov 11 15:01:26 2022

@author: Joona Räty

Solving pde system
    de=0
    d*e=0
explicitly on a 1+1-dimensional space-time mesh. e is located at the
mixed space-time edges and represented as a discrete one form.
This system corresponds to the wave equation and has a solution
    e = sin(x)cos(t)dx + cos(x)sin(t)dt
e is initilized at t=Lx, x=[0,Ly], by integrating the smooth solution (de rham map).
The solution is chosen such that it conforms with the natural von neumann boundary *e|_bound = 0
condition, which comes from incoplete dual forms at the boundary. No additional
consideration at the boundary is needed. The timestepping uses discrete versions
of de (when the value of all but one boundary edge of a face is known) and
d*e (when the value of all but one edge linked to a vertice is known / dual face is missing one edge)
For visualisation e is interpolated at mesh vertices and the
x-component sin(x)cos(t) is plotted after each update based on a dual face

requires overriding circumcenter in circumcenter.py with (cc according to minkowski metric)
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
from gmesh_recmesh import makerec, reqrecMsh
from pydec import simplicial_complex
import matplotlib.pyplot as plt
from scipy.sparse import diags

Lt = np.pi #length of mesh in the spatial direction, multiple of pi
Lx = np.pi #time period simulation
stretch  = 1
N=9#discretization level, meshpoints per Lx
visualise = False

#find out which edges have enough information to be updated via de=0 or d*e=0
def updateComplete(edg, marks, updates):
    #faces that are only missing one edge de=0
    for face in sc[1].d.getcol(edg).indices:
        faceEdges = sc[2].index_to_simplex[face].boundary()
        edgeInds = [sc[1].simplex_to_index[i] for i in faceEdges]
        if np.count_nonzero(marks[edgeInds,0]==3) == 2:#two of the three sides are marked as ready, see marks
            edgeToUpdate = np.where(marks[edgeInds,0]==0)
            if edgeToUpdate[0].size == 1:
                updates.append(edgeInds[edgeToUpdate[0][0]])
                marks[edgeInds[edgeToUpdate[0][0]],:] = [2,face]#2 -> upgrade with face rule de=0
    #nodes that are only missing one edge d*e=0
    for vert in sc[1].simplices[edg]:
        vertEdges = sc[0].d.getcol(vert).indices
        if np.count_nonzero(marks[vertEdges,0]==3) == vertEdges.size-1:#all but one of the node boundary edges are marked as ready
            edgeToUpdate = np.where(marks[vertEdges,0]==0)
            if edgeToUpdate[0].size == 1:
                updates.append(vertEdges[edgeToUpdate[0]][0])
                marks[vertEdges[edgeToUpdate],:] = [1,vert]#1 -> upgrade with vert rule d*e=0


#regular rectangular triangle mesh
V,E = reqrecMsh(Lt,Lx,Lx/N)
V[:,[1,0]]=V[:,[0,1]]#rotate mesh 90deg, otherwise there is no explicit timestepping scheme
V[:,1]*=stretch
sc = simplicial_complex(V,E)

initialVals = np.where(sc[1].circumcenter[:,1]>Lt*stretch-0.01)
initialEdgeBoundaries = sc[1].simplices[initialVals]

e = np.zeros(sc[1].num_simplices)
e[initialVals]=np.cos(V[initialEdgeBoundaries[:,0]][:,0])-np.cos(V[initialEdgeBoundaries[:,1]][:,0])

ax = plt.axes(projection='3d')
plt.title("Initial values")
ax.plot_trisurf(sc[1].circumcenter[:,0], sc[1].circumcenter[:,1], e,cmap='viridis', edgecolor='none')
plt.show()
'''
mark what state each edge is in:
    marks[i,0] = 3: value of edge i is known
    marks[i,0] = 2: all edges of face marks[i,1] except i are known
    marks[i,0] = 1: all edges of vertice marks[i,1] except i are known
'''
marks = np.zeros((sc[1].num_simplices,2), dtype='int32')
marks[initialVals,0] = 3

updates = [] #which edges can be solved
for edg in initialVals[0]:
    updateComplete(edg, marks, updates)

d0T = sc[0].d.T
d1 = sc[1].d
edgeVectors = V[sc[1].simplices[:,0]]-V[sc[1].simplices[:,1]]
timelikeEdges = 2*(edgeVectors[:,1]**2-edgeVectors[:,0]**2>0)-1
h1 = sc[1].star.diagonal()
h1d = h1*timelikeEdges#adjust h1 as h1[i,i]=\kappa |*e|/|e|
h1 = diags(h1d)

v = np.zeros(sc[0].num_simplices)
#Keep track of nodes which have all edges solved, for visualisation
solved = []
while len(updates)>0:
    i = updates.pop(0)
    j = marks[i,1]
    if marks[i, 0] == 1:
        e[i] = -d0T[j,i]/h1d[i]*d0T[j,:]*h1*e
        solved.append(j)
        #interpolate 1-form to vector field on vertice j, see räbinä chp 6
        left=0
        right=0
        for a in sc[0].d.getcol(j).indices:
            k=sc[1].simplices[a]
            I = V[k[1]]-V[k[0]]
            left+=h1d[a]*np.outer(I,I)
            right+=h1d[a]*I*e[a]
        R=np.linalg.solve(left,right)
        
        v[j] = R[0] #visualise dx component of one form e
        if visualise:
            fig = plt.figure(figsize=(10, 5), dpi=80)
            ax = fig.add_subplot(121,projection='3d')
            # ax.plot_trisurf(sc[1].circumcenter[:,0], sc[1].circumcenter[:,1], e,cmap='viridis', edgecolor='none')
            ax.plot_trisurf(V[:,0], V[:,1], v, cmap='viridis', edgecolor='none')
            ax = fig.add_subplot(122)
            plt.triplot(V[:,0], V[:,1], E)
            plt.plot(V[solved,0], V[solved,1], 'r.')
            plt.show()
    if marks[i,0] == 2:
        e[i] = -d1[j,i]*d1[j,:]*e
    marks[i,0] = 3#mark edge as solved
    updateComplete(i, marks, updates)
    
#for visualisation, some of the last vertices are solved via the face rule and dont get interpolated
remaining = np.where(V[:,1]<0.01)

for j in remaining[0]:
    left=0
    right=0
    for a in sc[0].d.getcol(j).indices:
        k=sc[1].simplices[a]
        I = V[k[1]]-V[k[0]]
        left+=h1d[a]*np.outer(I,I)
        right+=h1d[a]*I*e[a]
    R=np.linalg.solve(left,right)
    v[j] = R[0]
    
fig = plt.figure(figsize=(10, 5), dpi=80)
ax = fig.add_subplot(121,projection='3d')
ax.plot_trisurf(V[:,0], V[:,1], v, cmap='viridis', edgecolor='none')
ax = fig.add_subplot(122)
plt.triplot(V[:,0], V[:,1], E)
plt.plot(V[:,0], V[:,1], 'r.')
plt.show()
    
ax = plt.axes()
ax.tricontour(V[:,0], V[:,1], v)









