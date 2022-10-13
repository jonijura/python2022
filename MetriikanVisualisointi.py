# -*- coding: utf-8 -*-
"""
Created on Fri Oct  7 11:51:34 2022

@author: jonijura
"""

import numpy as np
import matplotlib.pyplot as plt


N = 20
lx = np.linspace(-1,1,N)
lt = np.linspace(-1, 1,N)
lx,lt = np.meshgrid(lx,lt);
d = lt**2-lx**2

fig = plt.figure(figsize=(10,10), dpi=80)
levels=np.sign(np.linspace(-1,1,15))*np.linspace(-1,1,15)**2
plt.contour(lx,lt, d, levels=levels)
