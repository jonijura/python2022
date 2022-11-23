# -*- coding: utf-8 -*-
"""
Created on Thu Oct 13 12:06:48 2022

@author: jonijura
"""


import numpy as np
import matplotlib.pyplot as plt

# V = np.loadtxt("C:\MyTemp\cpp\Samples\MeshGeneration2\\build\\vert.txt")
# E = np.loadtxt("C:\MyTemp\cpp\Samples\MeshGeneration2\\build\\tria.txt",dtype='int32')
# f = np.loadtxt("C:\MyTemp\cpp\Samples\MeshGeneration2\\build\\solutions.txt")

# plt.plot(f)
# plt.show()

# fig = plt.figure()
# ax = fig.add_subplot(121,projection='3d')
# ax.set_zlim( [-1,1])
# plt.title("Simulation result")
# ax.plot_trisurf(V[:,0], V[:,1], f,cmap='viridis', edgecolor='none')

# -*- coding: utf-8 -*-



solGlob = np.loadtxt("C:\MyTemp\cpp\Samples\MeshGeneration2\\build\\globalTimestepEnergy.txt")
solLoc1 = np.loadtxt("C:\MyTemp\cpp\Samples\MeshGeneration2\\build\\localTimestepEnergy.txt")

plt.plot(np.linspace(0,1,solGlob.size),solGlob)
plt.plot(np.linspace(0,1,solLoc1.size),solLoc1)
plt.legend(["global timestep", "local timestep1"])
plt.show()

# from scipy.sparse import csr_matrix
# from scipy.sparse.linalg import eigs
# mat = np.loadtxt("C:\MyTemp\cpp\Samples\MeshGeneration2\\build\\SystemMatrix.txt")
# x = mat[2].astype(np.int32)
# y = mat[1].astype(np.int32)
# size = max(x)+1
# matcsr = csr_matrix((mat[0], (x,y)), shape=(size,size))
# val, vec = eigs(matcsr,1)
# print(f"eigenvalue: {abs(val)[0]}")

# dt = 0.000384441
# i = np.arange(20000)
# fm = 5.0;
# lev = 20/(np.pi*fm);
# pl = np.exp(-(i*dt-3*lev)*(i*dt-3*lev)/(lev*lev))*np.sin(2*np.pi*fm*(i*dt-3*lev));
# plt.plot(pl)

#Fourier analyysi
refphist = np.loadtxt("C:\MyTemp\cpp\Samples\MeshGeneration2\\build\\GTSrefphistory.txt")
fmin=0.01
fmax=1.5
dt = 0.000224792
a = 1.0
N=refphist.shape[0]

f = (1/dt)*np.linspace(0,1,N) #Diskreettiä muunnosta vastaavat taajuudet
st = int(np.round(fmin*dt*N))
end = int(np.round(fmax*dt*N))

egmds = np.sqrt(np.array((3,6,9,11,12,14,19)))/(2*a)

p = abs(np.fft.fft(refphist[:,0]))**2
plt.plot(f[st:end],p[st:end])
# plt.show()
p = abs(np.fft.fft(refphist[:,1]))**2
plt.plot(f[st:end],p[st:end])
# plt.show()
p = abs(np.fft.fft(refphist[:,2]))**2
plt.plot(f[st:end],p[st:end])
# plt.plot(egmds, np.ones(egmds.size),'r.')
plt.show()

aoi = p[st:end]
foi = f[np.argwhere(aoi>1e7)]


tsmsh = np.loadtxt("C:\MyTemp\cpp\Samples\MeshGeneration2\\build\\tsmsh2.txt")
fig = plt.figure(figsize=(10, 8), dpi=80)
ax = fig.add_subplot(231,projection='3d')
ax.plot_trisurf(tsmsh[:,0], tsmsh[:,1], tsmsh[:,2],cmap='viridis', edgecolor='none')
plt.title("Simulation results")

ax = fig.add_subplot(232,projection='3d')
ax.plot_trisurf(tsmsh[:,0], tsmsh[:,1], tsmsh[:,3],cmap='viridis', edgecolor='none')

ax = fig.add_subplot(233,projection='3d')
ax.plot_trisurf(tsmsh[:,0], tsmsh[:,1], tsmsh[:,4],cmap='viridis', edgecolor='none')

ax = fig.add_subplot(234,projection='3d')
e1 = 1/np.sqrt(2) * np.cos(tsmsh[:,0]) * np.sin(tsmsh[:,1]) * np.sin(np.sqrt(2)*np.pi)
ax.plot_trisurf(tsmsh[:,0], tsmsh[:,1], tsmsh[:,2]-e1 ,cmap='viridis', edgecolor='none')
plt.title("Error")

e2 = 1/np.sqrt(2) * np.sin(tsmsh[:,0]) * np.cos(tsmsh[:,1]) * np.sin(np.sqrt(2)*np.pi)
ax = fig.add_subplot(235,projection='3d')
ax.plot_trisurf(tsmsh[:,0], tsmsh[:,1], tsmsh[:,3]-e2 ,cmap='viridis', edgecolor='none')

e3 = np.sin(tsmsh[:,0]) * np.sin(tsmsh[:,1]) * np.cos(np.sqrt(2)*np.pi)
ax = fig.add_subplot(236,projection='3d')
ax.plot_trisurf(tsmsh[:,0],tsmsh[:,1],  tsmsh[:,4]-e3,cmap='viridis', edgecolor='none')
plt.show()






















