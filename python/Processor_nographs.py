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
    sepcorstep=0.0041/3.13
    dz=0.01
    
    lossesx=[]
    lossespx=[]
    lostfirst=0
    
    for x0,px0 in itertools.izip(xlist,pxlist):
        x=x0-.0679-0.0002
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

def losses(lossdir):
    
    turns=[]    
    x=np.zeros(0)
    px=np.zeros(0)
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
            
    lossesx,lossespx,lostfirst = hitSeptum(x,px)
    xm=np.mean(x)
    pxm=np.mean(px)
    
    dx=(x-xm)
    dpx=(px-pxm)
    em= str(np.sqrt(np.mean(dx**2)*np.mean(dpx**2)-np.mean((dx)*(dpx))**2))
    extr= str(2.*len(turns)/len(os.listdir(lossdir)))
    sep1=str(float(lostfirst)/len(turns)*100)
    sepall=str(float(len(lossesx))/len(turns)*100)

    return em,extr,sep1,sepall
#    return 'a','b','c','d'
    
def scanlosses():
    rootdir='../Data/'
    f = open('scanplus2','w')
    f.write('turn\tamp\tem\textr\tsep1\tsepall\n')
    
    dirs=next(os.walk(rootdir))[1]
    dirs=[d for d in dirs if d.startswith('ripple')]
    for d in dirs:
        print d
        t=d.split('.')[0][8:]
        a=d.split('.')[1][3:]
        lossdir=os.path.join(rootdir,d)+'/losses/'
        if(len(os.listdir(lossdir))==0):
            print "No losses files in "+str(lossdir)
        else:
            em,extr,sep1,sepall = losses(lossdir)
            line=t+'\t'+a+'\t'+em+'\t'+extr+'\t'+sep1+'\t'+sepall+'\n'
            print line
            f.write(line)
        
    f.close()
        


scanlosses()