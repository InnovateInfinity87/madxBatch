# -*- coding: utf-8 -*-
"""
Created on Thu Aug 11 09:53:03 2016

@author: wvandepo
"""
import numpy as np
from matplotlib import rcParams

from scipy.interpolate import griddata
import matplotlib.pyplot as plt

rcParams.update({'font.size': 12})

def calculate_aspect(shape, extent):
    dx = (extent[1] - extent[0]) / float(shape[1])
    dy = (extent[3] - extent[2]) / float(shape[0])
    return dx / dy
    
data = np.loadtxt("./scanplus2",skiprows=1,unpack=True)

grid_x, grid_y = np.mgrid[0:700:120j, 0:1e-4:60j]

data[0]=data[0]**(-1)/0.0000233
#data[2]/=3.35e-8
#data one is in 10e-8 dk values, 1A = 8.4e-6 dk
data[1]/=(840.*1777.)

points=(data[0],data[1])
values=data[4]/data[5]
grid = griddata(points, values, (grid_x, grid_y), method='linear')


shape = (60,120)
extent = [0,700,0,1e-4]
    
plt.figure()

#plt.imshow(grid.T,extent=extent,aspect=calculate_aspect(shape, extent),origin='lower')
#
##plt.xscale("log")
##plt.yscale("log")
#cb1  = plt.colorbar(orientation='horizontal')
#cb1.set_label('Fraction of losses at first septum')
#plt.grid(True)
##plt.matshow(grid)
#plt.scatter(data[0],data[1],alpha=0.5,color="k")
##plt.gcf().set_size_inches(6, 6)
#plt.title(r"Fraction of losses at first septum")
#plt.xlabel("Hz")
#plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
#plt.ylabel(r"$\frac{\Delta I}{I}$",rotation=0)
#plt.tight_layout()

plt.scatter(data[1],data[4])
plt.xlabel("relative emittance")
plt.ylabel("Percent losses at ZS")

