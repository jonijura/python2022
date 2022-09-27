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


