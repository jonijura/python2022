# -*- coding: utf-8 -*-
"""
Created on Thu Oct 13 12:06:48 2022

@author: jonijura

jukan kirjastolla tuotettujen tuloksien analysointi 2+1 ja 1+1 aika-paikka verkko simuloinneille
"""


import numpy as np
import matplotlib.pyplot as plt

resultFolder = "C:\MyTemp\cpp\Samples\MeshGeneration2\\build\\"

#1+1 interpolated 1-form
tsmsh = np.loadtxt(resultFolder+"conv\\16_interpolated.txt")
fig = plt.figure(figsize=(8, 8), dpi=120)
ax = fig.add_subplot(221,projection='3d')
ax.plot_trisurf(tsmsh[:,0], tsmsh[:,1], tsmsh[:,2],cmap='viridis', edgecolor='none')
ax.set_xlabel('$x$')
ax.set_ylabel('$t$')
ax.set_zlabel('$F_x$')

ax = fig.add_subplot(222,projection='3d')
ax.plot_trisurf(tsmsh[:,0], tsmsh[:,1], tsmsh[:,3],cmap='viridis', edgecolor='none')
ax.set_xlabel('$x$')
ax.set_ylabel('$t$')
ax.set_zlabel('$F_t$')

e1 = np.sin(tsmsh[:,0]) * np.cos(tsmsh[:,1])
ax = fig.add_subplot(223,projection='3d')
ax.plot_trisurf(tsmsh[:,0], tsmsh[:,1], e1- tsmsh[:,2] ,cmap='viridis', edgecolor='none')
ax.set_xlabel('$x$')
ax.set_ylabel('$t$')
ax.set_zlabel('$\delta F_x$')

e2 = np.cos(tsmsh[:,0]) * np.sin(tsmsh[:,1])
ax = fig.add_subplot(224,projection='3d')
ax.plot_trisurf(tsmsh[:,0], tsmsh[:,1], e2- tsmsh[:,3] ,cmap='viridis', edgecolor='none')
ax.set_xlabel('$x$')
ax.set_ylabel('$t$')
ax.set_zlabel('$\delta F_t$')
plt.show()



#1+1 without interpolation as differential form (visualised at edge centers)
tsmsh2 = np.loadtxt(resultFolder+"conv\\16_1+1form.txt")
sol = np.cos(tsmsh2[:,0]) * np.cos(tsmsh2[:,1]) - np.cos(tsmsh2[:,2]) * np.cos(tsmsh2[:,3])
h1 = tsmsh2[:,0:2]-tsmsh2[:,2:4]
h1 = np.sum(np.abs(h1)**2,axis=-1)**(1./2)
diff = sol-tsmsh2[:,4]
print((diff*diff).mean())

ccs = tsmsh2[:,2:4] + 0.5*(tsmsh2[:,0:2]-tsmsh2[:,2:4])
fig = plt.figure(figsize=(12, 6), dpi=120)
ax = fig.add_subplot(121,projection='3d')
ax.plot_trisurf(ccs[:,0], ccs[:,1], diff*h1*diff,cmap='viridis', edgecolor='none')
ax.set_xlabel('$x$')
ax.set_ylabel('$t$')
ax.title.set_text("sqr diff hodge normed")

ax = fig.add_subplot(122,projection='3d')
ax.plot_trisurf(ccs[:,0], ccs[:,1], diff*diff,cmap='viridis', edgecolor='none')
ax.set_xlabel('$x$')
ax.set_ylabel('$t$')
ax.title.set_text("sqr diff")
plt.show()


#2+1 form with interpolation, dont forget to update reference time
tsmsh = np.loadtxt(resultFolder+"2+1interpolated.txt")
#reference solution, to compare to exact vector field solution, replace commented lines
exact_interpolated = np.loadtxt(resultFolder+"2+1exact_solution_interpolated.txt") 
ref_time = np.pi

fig = plt.figure(figsize=(11, 8), dpi=120)
ax = fig.add_subplot(231,projection='3d')
ax.plot_trisurf(tsmsh[:,0], tsmsh[:,1], tsmsh[:,2],cmap='viridis', edgecolor='none')
ax.set_xlabel('$x$')
ax.set_ylabel('$y$')
ax.set_zlabel('$F_x$')

ax = fig.add_subplot(232,projection='3d')
ax.plot_trisurf(tsmsh[:,0], tsmsh[:,1], tsmsh[:,3],cmap='viridis', edgecolor='none')
ax.set_xlabel('$x$')
ax.set_ylabel('$y$')
ax.set_zlabel('$F_y$')

ax = fig.add_subplot(233,projection='3d')
ax.plot_trisurf(tsmsh[:,0], tsmsh[:,1], tsmsh[:,4],cmap='viridis', edgecolor='none')
ax.set_xlabel('$x$')
ax.set_ylabel('$y$')
ax.set_zlabel('$F_t$')

e1 = exact_interpolated[:,2]#1/np.sqrt(2) * np.cos(tsmsh[:,0]) * np.sin(tsmsh[:,1]) * np.sin(np.sqrt(2)*ref_time)
ax = fig.add_subplot(234,projection='3d')
ax.plot_trisurf(tsmsh[:,0], tsmsh[:,1], e1 - tsmsh[:,2] ,cmap='viridis', edgecolor='none')
ax.set_xlabel('$x$')
ax.set_ylabel('$y$')
ax.set_zlabel('$\delta F_x$')


e2 = exact_interpolated[:,3]#1/np.sqrt(2) * np.sin(tsmsh[:,0]) * np.cos(tsmsh[:,1]) * np.sin(np.sqrt(2)*ref_time)
ax = fig.add_subplot(235,projection='3d')
ax.plot_trisurf(tsmsh[:,0], tsmsh[:,1], e2 - tsmsh[:,3] ,cmap='viridis', edgecolor='none')
ax.set_xlabel('$x$')
ax.set_ylabel('$y$')
ax.set_zlabel('$\delta F_y$')


e3 = exact_interpolated[:,4]#np.sin(tsmsh[:,0]) * np.sin(tsmsh[:,1]) * np.cos(np.sqrt(2)*ref_time)
ax = fig.add_subplot(236,projection='3d')
ax.plot_trisurf(tsmsh[:,0],tsmsh[:,1],  e3 - tsmsh[:,4],cmap='viridis', edgecolor='none')
ax.set_xlabel('$x$')
ax.set_ylabel('$y$')
ax.set_zlabel('$\delta F_t$')

plt.show()



#2+1 simulation values without interpolation
from scipy.stats import binned_statistic
noip = np.loadtxt(resultFolder+"2+1nointerpolation.txt")
sq2 = np.sqrt(2)
exact = -(np.sin(noip[:,0]) * np.sin(noip[:,1]) * np.sin(sq2*noip[:,2]) 
        - np.sin(noip[:,3]) * np.sin(noip[:,4]) * np.sin(sq2*noip[:,5]))/sq2
diff = exact - noip[:,6]
h1 = noip[:,0:3]-noip[:,3:6]
h1 = np.sum(np.abs(h1)**2,axis=-1)**(1./2)
aproxt = (noip[:,2]+noip[:,5])/2

binnum = 15
mean_stat = binned_statistic(aproxt, diff*diff, 
                             statistic='mean', 
                             bins=binnum , 
                             range=(0, np.pi))
plt.plot(aproxt, diff*diff, '.')
plt.title("diff squared")
plt.plot(np.linspace(0, np.pi, binnum ), 10*mean_stat.statistic)
plt.show()
mean_stat = binned_statistic(aproxt, diff*h1*diff, 
                             statistic='mean', 
                             bins=binnum , 
                             range=(0, np.pi))
plt.plot(aproxt, diff*h1*diff, '.')
plt.title("diff squared hodge normed")
plt.plot(np.linspace(0, np.pi, binnum ), 10*mean_stat.statistic)
plt.show()



#1+1 convergence test
arr = [(i+3)**2 for i in range(9)]
solt = np.zeros((len(arr),3))
for i,val in enumerate(arr):
    tsmsh2 = np.loadtxt(resultFolder+"conv\\"+str(val)+"_1+1form.txt")
    sol = np.cos(tsmsh2[:,0]) * np.cos(tsmsh2[:,1]) - np.cos(tsmsh2[:,2]) * np.cos(tsmsh2[:,3])
    h1 = tsmsh2[:,0:2]-tsmsh2[:,2:4]
    h1 = np.sum(np.abs(h1)**2,axis=-1)**(1./2)
    diff = sol-tsmsh2[:,4]
    solt[i][0] = max(abs(h1*diff))
    solt[i][1] = (diff*h1*diff).mean()
    solt[i][2] = np.dot(diff, h1*diff)
dx = np.array([0.541509,
                0.312709,
                0.202599,
                0.141636,
                0.104479,
                0.080202,
                0.063483,
                0.051488,
                0.042592])
fig = plt.figure(figsize=(13, 4), dpi=120)

ax = fig.add_subplot(131)
plt.yscale("log")
plt.xscale("log")
ax.plot(dx, solt[:,0])
ax.plot(dx,1e-2*dx**4)
ax.set_ylabel("error")
ax.legend(["maximum error","$O(n^4)$"])
ax = fig.add_subplot(132)
plt.yscale("log")
plt.xscale("log")
ax.plot(dx, solt[:,1])
ax.plot(dx,1e-2*dx**7)
ax.set_xlabel("Longest mesh edge")
ax.legend(["mean square error","$O(n^7)$"])
ax = fig.add_subplot(133)
plt.yscale("log")
plt.xscale("log")
ax.plot(dx, solt[:,2])
ax.plot(dx,1e-2*dx**5)
ax.legend(["error norm","$O(n^5)$"])

plt.show()















