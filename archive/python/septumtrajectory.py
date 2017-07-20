# -*- coding: utf-8 -*-
"""
Created on Thu Aug 11 17:09:32 2016

@author: wvandepo
"""
import numpy as np
import itertools
import matplotlib.pyplot as plt

steps=1877
dz=0.01

def hitSeptum(xlist,pxlist):
    global steps
    global dz
    a=2.7e-5
    xlist=np.asarray(xlist)
    pxlist=np.asarray(pxlist)
    steps=1877
    sepcorstep=0.0041/3.13
    dz=0.01

    
    for x0,px0 in itertools.izip(xlist,pxlist):
        x=x0-.0679-0.0001
        xtrack=[]
        for z in range(1,steps+1):
            b=(sepcorstep+a*z*dz+px0)
            x+=b*dz
            xtrack.append(x)
        return xtrack

plt.figure()
plt.plot(np.linspace(0,steps*dz,steps),hitSeptum([0.07],[-0.0016]))
plt.plot(np.linspace(0,steps*dz,steps),hitSeptum([0.069],[-0.0016]))
plt.plot(np.linspace(0,steps*dz,steps),hitSeptum([0.069],[-0.0013]))
plt.plot(np.linspace(0,steps*dz,steps),hitSeptum([0.075],[-0.00165]))
plt.axhline(0,color="k",linestyle ="--",lw=3)
plt.ylabel("[m]")
plt.xlabel("[m]")
plt.title("Trajectory in septum")
plt.ylim(-0.005,0.01)
plt.xlim(0,19.)
plt.show()       
    