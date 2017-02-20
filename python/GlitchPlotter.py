# -*- coding: utf-8 -*-
"""
Created on Fri Aug 19 09:22:04 2016

@author: wvandepo
"""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rcParams
from MADXreader import SingleLossFile
import Constants as c
import os  
import itertools

rcParams.update({'font.size': 9})

def hitSeptum(xlist,pxlist):
    a=2.7e-5
    xlist=np.asarray(xlist)
    pxlist=np.asarray(pxlist)
    steps=1877#313£1877
    sepcorstep=(67.7-63.34)/3130
    dz=0.01
    
    lossesx=[]
    lossespx=[]
    lostfirst=0
    
    for x0,px0 in itertools.izip(xlist,pxlist):
        x=x0-.0679-0.0001
        for z in range(1,steps+1):
            b=(sepcorstep+a*z*dz+px0)
            x+=b*dz
            if(x<0):
                lossesx.append(x0)
                lossespx.append(px0)
                #print "lost a particle"
                if(z*dz<3.13):
                    lostfirst+=1
                break
            elif(b>0):
                break
    return lossesx,lossespx,lostfirst

    
def losses():
    f, axarr = plt.subplots(3, 4,sharex=False, sharey=False)
        
    turns=[]    
    x=np.zeros(0)
    px=np.zeros(0)

    lossdir="../Data/LinSweep_p100k_turns34095/losses/"
    #lossdir="../Data/ripple_t84.0_a107.0_p50k/losses/"
    print "number of batches: " +str(len(os.listdir(lossdir)))
    
    for fn in os.listdir(lossdir): 
        lossbatch=SingleLossFile(lossdir+fn)
        
        if(lossbatch==-1):
            print lossdir+fn
            break
        turns+=lossbatch["TURN"].tolist()
        x=np.append(x,lossbatch['X'])
        px=np.append(px,lossbatch['PX'])
        
        boollist = [True if i == 1668.9776 else False for i in lossbatch["S"].tolist()]
        if(not all(boollist)):
            print fn

    print "number extracted particles: " +str(len(turns))
    print 'extracted percentage'+str(2.*len(turns)/len(os.listdir(lossdir)))+'%'
    lossesx,lossespx,lostfirst = hitSeptum(x,px)
    print "percentage of extr hits first septum: "+str(float(lostfirst)/len(turns)*100)
    print "percentage of extr hits a septum: "+str(float(len(lossesx))/len(turns)*100)
    print 'number ox particles extracted after first 1000 turns'+ str(sum(i > 1000 for i in turns))
    turns=np.asarray(turns)
    
    xm=np.mean(x)
    pxm=np.mean(px)
    
    dx=(x-xm)
    dpx=(px-pxm)
    eps2=np.mean(dx**2)*np.mean(dpx**2)-np.mean((dx)*(dpx))**2
    print eps2
    
    y,binEdges = np.histogram(turns,bins=int(c.Nturns/c.turnmultiplicity/8.))
    bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
    menStd     = np.sqrt(y)
    widthh      = c.turnmultiplicity*8.
    axarr[0,0].bar(bincenters, y, width=widthh,alpha=0.5, color='g', yerr=menStd,linewidth=0,label='Slow extraction rate')
    axarr[0,0].set_title('Linear tune sweep with 100k particles')
    axarr[0,0].set_xlabel("Beam intensity")
    axarr[0,0].set_ylabel("Numer of particles lost")
    axarr[0,0].legend()
    
    axarr[1,0].plot(x,px,'.',label=r'$\epsilon$ = %.3g' % float(np.sqrt(eps2)))
    axarr[1,0].plot(lossesx,lossespx,'y.',label=r'losses: = %.2f %%' % float(float(len(lossesx))/len(turns)*100))
    axarr[1,0].set_xlabel("X [m]")
    axarr[1,0].set_ylabel("PX")
    axarr[1,0].axvline(67.9/1000, color='black', label="Extraction aperture")
    axarr[1,0].set_ylim(-0.0018,-0.0013)
    axarr[1,0].set_xlim(0.067,0.085)
    axarr[1,0].legend()
    
    axarr[2,0].hist(turns, int(c.Nturns/c.turnmultiplicity/10.), normed=1,cumulative=-1, alpha = 0.5, color='g',linewidth=0, label= 'Extracted: %.2f %%' % float(2.*len(turns)/len(os.listdir(lossdir))))   
    axarr[2,0].set_xlabel("Turn number of particle extraction")
    axarr[2,0].set_ylabel("1- cumulative particles losses")
    axarr[2,0].legend()
    
    ###########################################################################
    
    turns=[]    
    x=np.zeros(0)
    px=np.zeros(0)

    lossdir="../Data/glitch_t500e500_a2000_p50k/losses/"
    #lossdir="../Data/ripple_t84.0_a107.0_p50k/losses/"
    print "number of batches: " +str(len(os.listdir(lossdir)))
    
    for fn in os.listdir(lossdir): 
        lossbatch=SingleLossFile(lossdir+fn)
        
        if(lossbatch==-1):
            print lossdir+fn
            break
        turns+=lossbatch["TURN"].tolist()
        x=np.append(x,lossbatch['X'])
        px=np.append(px,lossbatch['PX'])
        
        boollist = [True if i == 1668.9776 else False for i in lossbatch["S"].tolist()]
        if(not all(boollist)):
            print fn

    print "number extracted particles: " +str(len(turns))
    print 'extracted percentage'+str(2.*len(turns)/len(os.listdir(lossdir)))+'%'
    lossesx,lossespx,lostfirst = hitSeptum(x,px)
    print "percentage of extr hits first septum: "+str(float(lostfirst)/len(turns)*100)
    print "percentage of extr hits a septum: "+str(float(len(lossesx))/len(turns)*100)
    print 'number ox particles extracted after first 1000 turns'+ str(sum(i > 1000 for i in turns))
    turns=np.asarray(turns)
    
    xm=np.mean(x)
    pxm=np.mean(px)
    
    dx=(x-xm)
    dpx=(px-pxm)
    eps2=np.mean(dx**2)*np.mean(dpx**2)-np.mean((dx)*(dpx))**2
    print eps2
    
    y,binEdges = np.histogram(turns,bins=int(c.Nturns/c.turnmultiplicity/8.))
    bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
    menStd     = np.sqrt(y)
    widthh      = c.turnmultiplicity*8.
    axarr[0,1].bar(bincenters, y, width=widthh,alpha=0.5, color='g', yerr=menStd,linewidth=0,label='Slow extraction rate')
    axarr[0,1].set_title('glitch_t500e500_a2000_p50k')
    axarr[0,1].set_xlabel("Beam intensity")
    axarr[0,1].set_ylabel("Numer of particles lost")
    axarr[0,1].legend()
    
    axarr[1,1].plot(x,px,'.',label=r'$\epsilon$ = %.3g' % float(np.sqrt(eps2)))
    axarr[1,1].plot(lossesx,lossespx,'y.',label=r'losses: = %.2f %%' % float(float(len(lossesx))/len(turns)*100))
    axarr[1,1].set_xlabel("X [m]")
    axarr[1,1].set_ylabel("PX")
    axarr[1,1].axvline(67.9/1000, color='black', label="Extraction aperture")
    axarr[1,1].set_ylim(-0.0018,-0.0013)
    axarr[1,1].set_xlim(0.067,0.085)
    axarr[1,1].legend()
    
    axarr[2,1].hist(turns, int(c.Nturns/c.turnmultiplicity/10.), normed=1,cumulative=-1, alpha = 0.5, color='g',linewidth=0, label= 'Extracted: %.2f %%' % float(2.*len(turns)/len(os.listdir(lossdir))))   
    axarr[2,1].set_xlabel("Turn number of particle extraction")
    axarr[2,1].set_ylabel("1- cumulative particles losses")
    axarr[2,1].legend()
    
    
    ###########################################################################
    
    turns=[]    
    x=np.zeros(0)
    px=np.zeros(0)

    lossdir="../Data/glitch_t500e1000_a1000_p50k/losses/"
    #lossdir="../Data/ripple_t84.0_a107.0_p50k/losses/"
    print "number of batches: " +str(len(os.listdir(lossdir)))
    
    for fn in os.listdir(lossdir): 
        lossbatch=SingleLossFile(lossdir+fn)
        
        if(lossbatch==-1):
            print lossdir+fn
            break
        turns+=lossbatch["TURN"].tolist()
        x=np.append(x,lossbatch['X'])
        px=np.append(px,lossbatch['PX'])
        
        boollist = [True if i == 1668.9776 else False for i in lossbatch["S"].tolist()]
        if(not all(boollist)):
            print fn

    print "number extracted particles: " +str(len(turns))
    print 'extracted percentage'+str(2.*len(turns)/len(os.listdir(lossdir)))+'%'
    lossesx,lossespx,lostfirst = hitSeptum(x,px)
    print "percentage of extr hits first septum: "+str(float(lostfirst)/len(turns)*100)
    print "percentage of extr hits a septum: "+str(float(len(lossesx))/len(turns)*100)
    print 'number ox particles extracted after first 1000 turns'+ str(sum(i > 1000 for i in turns))
    turns=np.asarray(turns)
    
    xm=np.mean(x)
    pxm=np.mean(px)
    
    dx=(x-xm)
    dpx=(px-pxm)
    eps2=np.mean(dx**2)*np.mean(dpx**2)-np.mean((dx)*(dpx))**2
    print eps2
    
    y,binEdges = np.histogram(turns,bins=int(c.Nturns/c.turnmultiplicity/8.))
    bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
    menStd     = np.sqrt(y)
    widthh      = c.turnmultiplicity*8.
    axarr[0,3].bar(bincenters, y, width=widthh,alpha=0.5, color='g', yerr=menStd,linewidth=0,label='Slow extraction rate')
    axarr[0,3].set_title('glitch_t500e1000_a1000_p50k')
    axarr[0,3].set_xlabel("Beam intensity")
    axarr[0,3].set_ylabel("Numer of particles lost")
    axarr[0,3].legend()
    
    axarr[1,3].plot(x,px,'.',label=r'$\epsilon$ = %.3g' % float(np.sqrt(eps2)))
    axarr[1,3].plot(lossesx,lossespx,'y.',label=r'losses: = %.2f %%' % float(float(len(lossesx))/len(turns)*100))
    axarr[1,3].set_xlabel("X [m]")
    axarr[1,3].set_ylabel("PX")
    axarr[1,3].axvline(67.9/1000, color='black', label="Extraction aperture")
    axarr[1,3].set_ylim(-0.0018,-0.0013)
    axarr[1,3].set_xlim(0.067,0.085)
    axarr[1,3].legend()
    
    axarr[2,3].hist(turns, int(c.Nturns/c.turnmultiplicity/10.), normed=1,cumulative=-1, alpha = 0.5, color='g',linewidth=0, label= 'Extracted: %.2f %%' % float(2.*len(turns)/len(os.listdir(lossdir))))   
    axarr[2,3].set_xlabel("Turn number of particle extraction")
    axarr[2,3].set_ylabel("1- cumulative particles losses")
    axarr[2,3].legend()
    
        ###########################################################################
    
    turns=[]    
    x=np.zeros(0)
    px=np.zeros(0)

    lossdir="../Data/glitch_t500e500_a1000_p50k/losses/"
    #lossdir="../Data/ripple_t84.0_a107.0_p50k/losses/"
    print "number of batches: " +str(len(os.listdir(lossdir)))
    
    for fn in os.listdir(lossdir): 
        lossbatch=SingleLossFile(lossdir+fn)
        
        if(lossbatch==-1):
            print lossdir+fn
            break
        turns+=lossbatch["TURN"].tolist()
        x=np.append(x,lossbatch['X'])
        px=np.append(px,lossbatch['PX'])
        
        boollist = [True if i == 1668.9776 else False for i in lossbatch["S"].tolist()]
        if(not all(boollist)):
            print fn

    print "number extracted particles: " +str(len(turns))
    print 'extracted percentage'+str(2.*len(turns)/len(os.listdir(lossdir)))+'%'
    lossesx,lossespx,lostfirst = hitSeptum(x,px)
    print "percentage of extr hits first septum: "+str(float(lostfirst)/len(turns)*100)
    print "percentage of extr hits a septum: "+str(float(len(lossesx))/len(turns)*100)
    print 'number ox particles extracted after first 1000 turns'+ str(sum(i > 1000 for i in turns))
    turns=np.asarray(turns)
    
    xm=np.mean(x)
    pxm=np.mean(px)
    
    dx=(x-xm)
    dpx=(px-pxm)
    eps2=np.mean(dx**2)*np.mean(dpx**2)-np.mean((dx)*(dpx))**2
    print eps2
    
    y,binEdges = np.histogram(turns,bins=int(c.Nturns/c.turnmultiplicity/8.))
    bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
    menStd     = np.sqrt(y)
    widthh      = c.turnmultiplicity*8.
    axarr[0,2].bar(bincenters, y, width=widthh,alpha=0.5, color='g', yerr=menStd,linewidth=0,label='Slow extraction rate')
    axarr[0,2].set_title('glitch_t500e500_a1000_p50k')
    axarr[0,2].set_xlabel("Beam intensity")
    axarr[0,2].set_ylabel("Numer of particles lost")
    axarr[0,2].legend()
    
    axarr[1,2].plot(x,px,'.',label=r'$\epsilon$ = %.3g' % float(np.sqrt(eps2)))
    axarr[1,2].plot(lossesx,lossespx,'y.',label=r'losses: = %.2f %%' % float(float(len(lossesx))/len(turns)*100))
    axarr[1,2].set_xlabel("X [m]")
    axarr[1,2].set_ylabel("PX")
    axarr[1,2].axvline(67.9/1000, color='black', label="Extraction aperture")
    axarr[1,2].set_ylim(-0.0018,-0.0013)
    axarr[1,2].set_xlim(0.067,0.085)
    axarr[1,2].legend()
    
    axarr[2,2].hist(turns, int(c.Nturns/c.turnmultiplicity/10.), normed=1,cumulative=-1, alpha = 0.5, color='g',linewidth=0, label= 'Extracted: %.2f %%' % float(2.*len(turns)/len(os.listdir(lossdir))))   
    axarr[2,2].set_xlabel("Turn number of particle extraction")
    axarr[2,2].set_ylabel("1- cumulative particles losses")
    axarr[2,2].legend()
    
    plt.tight_layout()
    plt.show()
    ############################################################################




losses()