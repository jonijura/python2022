# -*- coding: utf-8 -*-
"""
Created on Thu Sep 29 15:53:20 2022

@author: jonijura
"""

import numpy as np
import matplotlib.pyplot as plt
from gmesh_recmesh import makerec
from pydec import simplicial_complex

Nxy = 5
bxy = np.pi
dxy = bxy/Nxy

circ = 15

V,E = makerec(bxy,bxy,dxy)
sc = simplicial_complex(V,E)

plt.figure(figsize=(8, 8), dpi=80)

plt.triplot(V[:,0], V[:,1], E)
plt.plot(sc[2].circumcenter[:,0],sc[2].circumcenter[:,1],'.')

r = np.linalg.norm(V[sc[0].simplices[sc[2].simplices[circ][1]][0]]-sc[2].circumcenter[circ,:])
circle2 = plt.Circle((sc[2].circumcenter[circ,0],sc[2].circumcenter[circ,1]), r, color='b', fill=False)
ax = plt.gca()
ax.add_patch(circle2)

plt.show()