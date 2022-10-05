# -*- coding: utf-8 -*-
"""
Created on Mon Sep 26 15:11:06 2022

@author: jonijura
"""

import gmsh
import numpy as np


def makerec(h,w,r):
    gmsh.initialize()
    gmsh.model.add("rec")
    
    gmsh.model.geo.addPoint(0, 0, 0, r, 1)
    gmsh.model.geo.addPoint(w, 0, 0, r, 2)
    gmsh.model.geo.addPoint(w, h, 0, r, 3)
    gmsh.model.geo.addPoint(0, h, 0, r, 4)

    gmsh.model.geo.addLine(1, 2, 1)
    gmsh.model.geo.addLine(2, 3, 2)
    gmsh.model.geo.addLine(3, 4, 3)
    gmsh.model.geo.addLine(4, 1, 4)

    gmsh.model.geo.addCurveLoop([1,2,3,4], 1)

    gmsh.model.geo.addPlaneSurface([1], 1)
    
    gmsh.model.geo.synchronize()

    gmsh.model.mesh.generate(2)
    #gmsh.model.mesh.refine()
    #gmsh.model.mesh.optimize("Laplace2D")

    b = gmsh.model.mesh.getNodes()[1]
    B = np.reshape(b,(int(len(b)/3),3))
    V = np.delete(B,2,1)

    c = gmsh.model.mesh.getElements(2)[2]
    E = np.reshape(c, (int(len(c[0])/3),3))-1

    gmsh.finalize()
    return (V,E.astype('int32'))

def makecube(w,h,d,r):
    gmsh.initialize()
    gmsh.model.add("cub")
    
    m = gmsh.model.geo
    
    m.addPoint(0, 0, 0, r, 1)
    m.addPoint(w, 0, 0, r, 2)
    m.addPoint(w, h, 0, r, 3)
    m.addPoint(0, h, 0, r, 4)
    m.addPoint(0, 0, d, r, 5)
    m.addPoint(w, 0, d, r, 6)
    m.addPoint(w, h, d, r, 7)
    m.addPoint(0, h, d, r, 8)

    m.addLine(1, 2, 1)
    m.addLine(2, 3, 2)
    m.addLine(3, 4, 3)
    m.addLine(4, 1, 4)
    m.addLine(1, 5, 5)
    m.addLine(2, 6, 6)
    m.addLine(3, 7, 7)
    m.addLine(4, 8, 8)
    m.addLine(6, 5, 9)
    m.addLine(7, 6, 10)
    m.addLine(8, 7, 11)
    m.addLine(5, 8, 12)

    m.addCurveLoop([1,2,3,4], 1)
    m.addPlaneSurface([1], 1)
    m.addCurveLoop([6,-10,-7,-2], 2)
    m.addPlaneSurface([2], 2)
    m.addCurveLoop([7,-11,-8,-3], 3)
    m.addPlaneSurface([3], 3)
    m.addCurveLoop([12,11,10,9], 4)
    m.addPlaneSurface([4], 4)
    m.addCurveLoop([8,-12,-5,-4], 5)
    m.addPlaneSurface([5], 5)
    m.addCurveLoop([-1,5,-9,-6], 6)
    m.addPlaneSurface([6], 6)
    
    m.addSurfaceLoop([1,2,3,4,5,6],1)
    m.addVolume([1],1)
    
    gmsh.model.geo.synchronize()

    gmsh.model.mesh.generate(3)
    #gmsh.model.mesh.refine()
    #gmsh.model.mesh.optimize("Netgen")
    #gmsh.fltk.run()

    b = gmsh.model.mesh.getNodes()[1]
    V = np.reshape(b,(int(len(b)/3),3))

    c = gmsh.model.mesh.getElements(3)[2]
    E = np.reshape(c, (int(len(c[0])/4),4))-1

    gmsh.finalize()
    return (V,E.astype('int32'))

def reqrecMsh(w,h,r):
    nx = int(np.ceil(w/r))+1
    dx = w/(nx-1)
    ny = int(np.round(h/(0.866*dx)))+1
    dy = h/(ny-1)
    evenr = np.linspace(0, w, nx)
    oddr = np.linspace(-dx/2, w+dx/2,nx+1)
    oddr[0]+=dx/2
    oddr[-1]-=dx/2
    V = np.empty((0,2))
    for i in range(ny):
        if i%2==0:
            row = np.vstack((evenr, i*dy*np.ones(nx))).T
        else:
            row = np.vstack((oddr, i*dy*np.ones(nx+1))).T
        V = np.vstack((V,row))
    T=np.array([[0,nx,nx+1],[0,nx+1,1]])
    RE=np.empty((0,3))
    for i in range(nx):
        RE=np.vstack((RE,T+i))
    E=np.empty((0,3))
    for i in range(ny-1):
        if i%2==0:
            row = RE[:-1]+int(i/2)*(2*nx+1)
        else:
            row = RE[1:]+int((i-1)/2)*(2*nx+1)+nx
        E=np.vstack((E,row))
    return (V,E.astype('int32'))

# =============================================================================
# plt.figure(figsize=(8, 8), dpi=80)    
# V,E = reqrecMsh(1,1,0.3)
# plt.triplot(V[:,0], V[:,1], E)
# =============================================================================
    

# =============================================================================
# V,E = makerec(np.pi,np.pi, 0.3)
# sc = SimplicialMesh(V,E)
# figure(figsize=(8, 8), dpi=80)
# plt.triplot(sc.vertices[:,0], sc.vertices[:,1], sc.indices)
# plt.show()
# =============================================================================


# =============================================================================
# gmsh.initialize()
# gmsh.model.add("rec")
# a = gmsh.model.occ.addRectangle(0,0,0,np.pi,np.pi)
# 
# gmsh.model.occ.synchronize()
# 
# gmsh.model.mesh.generate(2)
# #gmsh.model.mesh.refine()
# #gmsh.model.mesh.optimize("Laplace2D")
# 
# 
# 
# 
# b = gmsh.model.mesh.getNodes()[1]
# B = np.reshape(b,(int(len(b)/3),3))
# V = np.delete(B,2,1)
# 
# c = gmsh.model.mesh.getElements(2)[2]
# E = np.reshape(c, (int(len(c[0])/3),3))-1
# 
# sc = SimplicialMesh(V,E)
# figure(figsize=(8, 8), dpi=80)
# plt.triplot(sc.vertices[:,0], sc.vertices[:,1], sc.indices)
# plt.show()
# #gmsh.fltk.run()
# gmsh.finalize()
# =============================================================================


