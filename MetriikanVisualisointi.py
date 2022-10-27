# -*- coding: utf-8 -*-
"""
Created on Fri Oct  7 11:51:34 2022

@author: jonijura
"""

import numpy as np
import matplotlib.pyplot as plt


# N = 20
# lx = np.linspace(-1,1,N)
# lt = np.linspace(-1, 1,N)
# lx,lt = np.meshgrid(lx,lt);
# d = lt**2-lx**2

# fig = plt.figure(figsize=(10,10), dpi=80)
# levels=np.sign(np.linspace(-1,1,15))*np.linspace(-1,1,15)**2
# plt.contour(lx,lt, d, levels=levels)


from scipy.spatial import Voronoi, voronoi_plot_2d
from scipy.spatial import Delaunay
from gmesh_recmesh import reqrecMsh

points = np.random.rand(20,2)
V,E = reqrecMsh(1,1,0.4)
# points = V

vor = Voronoi(points)
plt.figure(figsize=(10,10), dpi=80)
voronoi_plot_2d(vor, plt.axes())

# tri = Delaunay(points)
# plt.triplot(points[:,0], points[:,1], tri.simplices)
# plt.plot(points[:,0], points[:,1], 'o')
plt.show()

vert = [0<a[0] and a[0]<1 and 0<a[1] and a[1]<1 for a in vor.vertices]
vor = Voronoi(vor.vertices[vert])
plt.figure(figsize=(10,10), dpi=80)
voronoi_plot_2d(vor, plt.axes())
