# -*- coding: utf-8 -*-
"""
Created on Mon Sep 26 15:11:06 2022

@author: jonijura
"""

import gmsh
import numpy as np

import matplotlib.pyplot as plt
from pydec import SimplicialMesh
from matplotlib.pyplot import figure


gmsh.initialize()
gmsh.model.add("rec")
gmsh.model.occ.addRectangle(0,0,0,1,1)
gmsh.model.occ.synchronize()

gmsh.model.mesh.generate(2)
gmsh.model.mesh.refine()
gmsh.model.mesh.optimize("Laplace2D")




b = gmsh.model.mesh.getNodes()[1]
B = np.reshape(b,(int(len(b)/3),3))
V = np.delete(B,2,1)

c = gmsh.model.mesh.getElements(2)[2]
E = np.reshape(c, (int(len(c[0])/3),3))-1

sc = SimplicialMesh(V,E)
figure(figsize=(8, 8), dpi=80)
plt.triplot(sc.vertices[:,0], sc.vertices[:,1], sc.indices)
plt.show()
gmsh.fltk.run()
gmsh.finalize()