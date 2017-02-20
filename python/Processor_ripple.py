# -*- coding: utf-8 -*-
"""
Created on Mon Jul 11 13:46:53 2016

@author: wvandepo
"""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rcParams
from MADXreader import SingleLossFile
import Constants as c
import os  
import itertools

rcParams.update({'font.size': 7})

def hitSeptum(xlist,pxlist):
    a=2.7e-5
    xlist=np.asarray(xlist)
    pxlist=np.asarray(pxlist)
    steps=313
    sepcorstep=0.0041/3.13
    dz=0.01
    
    lossesx=[]
    lossespx=[]
    
    for x0,px0 in itertools.izip(xlist,pxlist):
        x=x0-.0679-0.0002
        for z in range(1,steps+1):
            b=(sepcorstep+a*z*dz+px0)
            x+=b*dz
            if(x<0):
                lossesx.append(x0)
                lossespx.append(px0)
                #print "lost a particle"
                break
            elif(b>0):
                break
    return lossesx,lossespx

def losses():
    f, axarr = plt.subplots(3, 2,sharex=False, sharey=False)
        
    turns=[]    
    x=np.zeros(0)
    px=np.zeros(0)

    lossdir="../Data/ripple_t69.0_a126.0_p10k/losses/"
    print "number of batches: " +str(len(os.listdir(lossdir)))
    
    for fn in os.listdir(lossdir): 
        lossbatch=SingleLossFile(lossdir+fn)
        if(lossbatch ==-1):
            print lossdir+fn
        else:
            
            turns+=lossbatch["TURN"].tolist()
            x=np.append(x,lossbatch['X'])
            px=np.append(px,lossbatch['PX'])
            
            boollist = [True if i == 1668.9776 else False for i in lossbatch["S"].tolist()]
            if(not all(boollist)):
                print fn
    print len(turns)
    lossesx,lossespx = hitSeptum(x,px)
    print "percentage of extr hits septum: "+str(float(len(lossesx))/len(turns)*100)
    print sum(i > 1000 for i in turns)
    turns=np.asarray(turns)
    
    xm=np.mean(x)
    pxm=np.mean(px)
    
    dx=(x-xm)
    dpx=(px-pxm)
    eps2=np.mean(dx**2)*np.mean(dpx**2)-np.mean((dx)*(dpx))**2
    print eps2
    
    y,binEdges = np.histogram(turns,bins=int(c.Nturns/c.turnmultiplicity/10.))
    bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
    menStd     = np.sqrt(y)
    width      = c.turnmultiplicity*10
    axarr[0,0].bar(bincenters, y, width=width,alpha=0.5, color='g', yerr=menStd,linewidth=0)
    axarr[0,0].set_title('Without glitch 10k particles')
    axarr[0,0].set_xlabel("Turn number of particle extraction")
    axarr[0,0].set_ylabel("Numer of particles lost")
    
    axarr[1,0].plot(x,px,'.',label=r'$\epsilon$ ='+str(np.sqrt(eps2)))
    axarr[1,0].plot(lossesx,lossespx,'y.',label=r'$losses:$ ='+str(float(len(lossesx))/len(turns)*100)+'%')
    axarr[1,0].set_xlabel("X [m]")
    axarr[1,0].set_ylabel("PX")
    axarr[1,0].axvline(67.9/1000, color='black', label="ZS.UP")
    
    axarr[2,0].hist(turns, int(c.Nturns/c.turnmultiplicity/10.), normed=1,cumulative=-1, alpha = 0.5, color='g',linewidth=0, label= str(2.*len(turns)/len(os.listdir(lossdir)))+'% extracted')

    
    axarr[2,0].set_xlabel("Turn number of particle extraction")
    axarr[2,0].set_ylabel("1- cumulative particles lost")
    axarr[1,0].legend()
    axarr[2,0].legend()
    
    ############################################################################
    
    turns=[]    
    x=np.zeros(0)
    px=np.zeros(0)
    
    lossdir="../Data/ripple_t69.0_a126.0_p10k/losses/"
    print "number of batches: " +str(len(os.listdir(lossdir)))
    for fn in os.listdir(lossdir): 
        lossbatch=SingleLossFile(lossdir+fn)
        if(lossbatch ==-1):
            print lossdir+fn
        else:
            
            turns+=lossbatch["TURN"].tolist()
            x=np.append(x,lossbatch['X'])
            px=np.append(px,lossbatch['PX'])
            
            boollist = [True if i == 1668.9776 else False for i in lossbatch["S"].tolist()]
            if(not all(boollist)):
                print fn
    print len(turns)
    lossesx,lossespx = hitSeptum(x,px)
    print "percentage of extr hits septum: "+str(float(len(lossesx))/len(turns)*100)
    print sum(i > 1000 for i in turns)
    turns=np.asarray(turns)
    
    xm=np.mean(x)
    pxm=np.mean(px)
    
    dx=(x-xm)
    dpx=(px-pxm)
    eps2=np.mean(dx**2)*np.mean(dpx**2)-np.mean((dx)*(dpx))**2
    print eps2

    y,binEdges = np.histogram(turns,bins=int(c.Nturns/c.turnmultiplicity/10.))
    bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
    menStd     = np.sqrt(y)
    axarr[0,1].bar(bincenters, y, width=width,alpha=0.5, color='g', yerr=menStd,linewidth=0)
    axarr[0,1].set_title('realglitch_10k_0Comp')
    axarr[0,1].set_xlabel("Turn number of particle extraction")
    axarr[0,1].set_ylabel("Numer of particles lost")
    
    axarr[1,1].plot(x,px,'.',label=r'$\epsilon$ ='+str(np.sqrt(eps2)))
    axarr[1,1].plot(lossesx,lossespx,'y.',label=r'$losses:$ ='+str(float(len(lossesx))/len(turns)*100)+'%')
    axarr[1,1].set_xlabel("X [m]")
    axarr[1,1].set_ylabel("PX")
    axarr[1,1].axvline(67.9/1000, color='black', label="ZS.UP")
    
    axarr[2,1].hist(turns, int(c.Nturns/c.turnmultiplicity/10.), normed=1,cumulative=-1, alpha = 0.5, color='g',linewidth=0, label= str(2.*len(turns)/len(os.listdir(lossdir)))+'% extracted')
    axarr[1,1].legend()
    axarr[2,1].legend()
    axarr[2,1].set_xlabel("Turn number of particle extraction")
    axarr[2,1].set_ylabel("1- cumulative particles lost")
    
    ###########################################################################

losses()