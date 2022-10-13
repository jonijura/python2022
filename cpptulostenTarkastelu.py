# -*- coding: utf-8 -*-
"""
Created on Thu Oct 13 12:06:48 2022

@author: jonijura
"""


import numpy as np
import matplotlib.pyplot as plt

V = np.loadtxt("C:\MyTemp\cpp\Samples\MeshGeneration2\\build\\vert.txt")
E = np.loadtxt("C:\MyTemp\cpp\Samples\MeshGeneration2\\build\\tria.txt",dtype='int32')
f = np.loadtxt("C:\MyTemp\cpp\Samples\MeshGeneration2\\build\\solutions.txt")

plt.plot(f)

# fig = plt.figure()
# ax = fig.add_subplot(121,projection='3d')
# ax.set_zlim( [-1,1])
# plt.title("Simulation result")
# ax.plot_trisurf(V[:,0], V[:,1], f,cmap='viridis', edgecolor='none')