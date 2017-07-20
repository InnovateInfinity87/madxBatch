# -*- coding: utf-8 -*-
"""
Created on Mon Jul 11 13:46:53 2016

@author: wvandepo
"""


import matplotlib
#matplotlib.use("Agg")
import matplotlib.pyplot as plt
#import matplotlib.animation as animation
import numpy as np
from matplotlib import rcParams
#from mpl_toolkits.mplot3d import Axes3D
from MADXreader import TrackRead,TwissRead,SingleTrack,SingleLossFile
import Constants as c
import os  
import itertools

rcParams.update({'font.size': 9})

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
        x=x0-.0679-0.0001
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
#            else:
#                print "approaching septum"+str(b)
    return lossesx,lossespx
            
    
def Reader():
    elements=['"AP.UP.ZS21633"','"ZS.21633"','"AP.DO.ZS21633"']
    plotplace=0 #0=start
    
    initial = np.loadtxt(c.data+"initial_distribution_gauss.txt",skiprows=7,unpack=True)
    inidpp = initial[7]
    
    variable_data, parameters = TwissRead('twiss_after_thinning.prt')
    
    elementindex=variable_data['NAME'].index(elements[plotplace])
    sqrtbetax = np.sqrt(variable_data['BETX'][elementindex])
    alfax = variable_data['ALFX'][elementindex]
    
    amp=[]
    dpp=[]
    turns=[]   
    
    for fn in os.listdir(c.tracksdir): 
        place =  int(fn.split(".")[2][3:])
        if(place==plotplace+1):  #at start place=1
            #print fn
            batch =  int(fn.split(".")[1][5:])
            particle =  int(fn.split(".")[3][1:])
            track=SingleTrack(fn)
            turns.append(track["TURN"][-1])
    
            xallm=np.mean(np.nan_to_num(track["X"]))
            pxallm=np.mean(np.nan_to_num(track["PX"]))
            
            normx=(track['X'][0]-xallm)/sqrtbetax
            normpx=normx*alfax + sqrtbetax*(track['PX'][0]-pxallm)
                
            amp.append(np.sqrt(normx**2+normpx**2))
            dpp.append(inidpp[particle-1+batch*c.Nparperbatch])
    
    return dpp, amp, turns 
    
def Plotter(dpp, amp, turns):
    fig = plt.figure()
    ax3d = fig.add_subplot(111, projection='3d')
    
     
    
    for x,y,z in zip(dpp,amp,turns):
        ax3d.scatter(x, y, z)
    ax3d.set_xlabel(r'$\Delta p$')
    ax3d.set_ylabel('Normalised amplitude')
    ax3d.set_zlabel('Turns before extraction')
    
    
    plt.title("Slow extraction: 34095 turns, dpp=0.5e-3, 7500 particles")
    plt.show()
    
    f, ax = plt.subplots(1, 3,sharex=False, sharey=False)
    
    ax[0].plot(dpp,amp,".")
    ax[0].set_xlabel(r'$\Delta p$')
    ax[0].set_ylabel('Amplitude')
    
    ax[1].plot(dpp,turns,".")
    ax[1].set_xlabel(r'$\Delta p$')
    ax[1].set_ylabel('Turns')
    
    ax[2].plot(amp,turns,".")
    ax[2].set_xlabel('Amplitude')
    ax[2].set_ylabel('Turns')
    
    plt.title("Slow extraction: 34095 turns, dpp=0.5e-3, 7500 particles")
    plt.show()
    
    
    
def plotPhaseSpace():
    
    elements=['"AP.UP.ZS21633"','"ZS.21633"','"AP.DO.ZS21633"']
    places=len(elements)
    Nbatches=1
    Nparperbatch=10
    
    
    ##### Load the intialdistributions #####
    initial = np.loadtxt("../input/initial_distribution_gauss.txt",skiprows=7,unpack=True)
    inidpp = initial[7]
    
    ##### Load all tracks #####
    if(c.whichbatch==-1):
        Tracks = TrackRead()
        batchloop=range(Nbatches)
        print 'Tracks loaded in RAM'
    
    ##### Load some tracks ####
    else:
        Tracks = {}
        batchloop=[c.whichbatch]
        for fn in os.listdir(c.tracksdir):
            batch =  int(fn.split(".")[1][5:])
            if(batch==c.whichbatch):
                place =  int(fn.split(".")[2][3:])
                particle =  int(fn.split(".")[3][1:])
                Tracks[(batch,place,particle)]=SingleTrack(fn)
        print 'Tracks loaded in RAM'
    
    
    #Tracks in phase space
    f, axarr = plt.subplots(places, 2,sharex=True, sharey=True)
    f.suptitle("Tracks in phase space", fontsize=12)
    print "Start plotting of tracks"
    
    
    for place in range(places):
        for b in batchloop:
            for p in range(Nparperbatch):
                this=Tracks[(b,place+1,p+1)]
                axarr[place, 0].plot(this['X'],this['PX'],".")
                axarr[place, 0].set_title('Tracking at '+elements[place])
                axarr[place, 1].plot(this['Y'],this['PY'],".",label=str("dpp= '%.5g , t: %d" % (inidpp[b*c.Nparperbatch+p],len(this['X'])*c.turnmultiplicity) ))
                axarr[place,0].set_ylabel(r'$P_x$')
                axarr[place,1].set_ylabel(r'$P_y$')
                axarr[place, 0].axhline(0, color='black')
                axarr[place, 0].axvline(0, color='black')
                axarr[place,1].axhline(0, color='black')
                axarr[place,1].axvline(0, color='black')
    
    axarr[place,0].set_xlabel('x')
    axarr[place,1].set_xlabel('y')
    plt.legend()
    plt.tight_layout()
    
    #Normalised coordinates
    variable_data, parameters = TwissRead('twiss_after_thinning.prt')
    #print variable_data['NAME']
    
    #Tracks in normalised phase space
    
    f, axarr = plt.subplots(places, 2,sharex=False, sharey=False)
    f.suptitle("Tracks in normalised phase space", fontsize=14)
    print "Start plotting of normalised tracks"
    
    for place in range(places):
        
        elementindex=variable_data['NAME'].index(elements[place])
        sqrtbetax = np.sqrt(variable_data['BETX'][elementindex])
        sqrtbetay = np.sqrt(variable_data['BETY'][elementindex])
        alfax = variable_data['ALFX'][elementindex]
        alfay = variable_data['ALFY'][elementindex]
        
        xall =  np.empty(0)
        yall =  np.empty(0)
        pxall = np.empty(0)
        pyall = np.empty(0)
        
        for b in batchloop:
            for p in range(Nparperbatch):
                this=Tracks[(b,place+1,p+1)]
                xall = np.append(xall,this['X'])
                yall = np.append(yall,this['Y'])
                pxall = np.append(pxall,this['PX'])
                pyall = np.append(pyall,this['PY'])
        xallm=np.mean(np.nan_to_num(xall))
        yallm=np.mean(np.nan_to_num(yall))
        pxallm=np.mean(np.nan_to_num(pxall))
        pyallm=np.mean(np.nan_to_num(pyall))
        
        for b in batchloop:
            for p in range(Nparperbatch):
                this=Tracks[(b,place+1,p+1)]
                
                normx=(this['X']-xallm)/sqrtbetax
                normy=(this['Y']-yallm)/sqrtbetay
                normpx=normx*alfax + sqrtbetax*(this['PX']-pxallm)
                normpy=normy*alfay + sqrtbetay*(this['PY']-pyallm)
                
                axarr[place, 0].plot(normx+xallm,normpx+pxallm,".")
                axarr[place, 0].set_title('Tracking at '+elements[place])
                axarr[place, 1].plot(normy+yallm,normpy+pyallm,".")
                axarr[place,0].set_ylabel(r'norm $P_X$')
                axarr[place,1].set_ylabel(r'norm $P_Y$')
                axarr[place, 0].axhline(0, color='black')
                axarr[place, 0].axvline(0, color='black')
                axarr[place,1].axhline(0, color='black')
                axarr[place,1].axvline(0, color='black')
    
    axarr[place,0].set_xlabel('norm X')
    axarr[place,1].set_xlabel('norm Y')
    
    plt.tight_layout()


def losses():
    
    turns=[]    
    normx=np.zeros(0)
    normpx=np.zeros(0)
    x=np.zeros(0)
    px=np.zeros(0)

    
    lossdir="../Data/ripple_t69.0_a146.0_p50k/losses/"
    print "number of batches: " +str(len(os.listdir(lossdir)))
    f, axarr = plt.subplots(3, 2,sharex=False, sharey=False)
    
    for fn in os.listdir(lossdir): 
        lossbatch=SingleLossFile(lossdir+fn)
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
    #xglitch = x[np.intersect1d(np.where(turns>5000),np.where(turns<7000))]	
    #pxglitch = px[np.intersect1d(np.where(turns>5000),np.where(turns<7000))]
    #print len(xglitch)
    
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
    #axarr[1,0].plot(xglitch,pxglitch,'r.')
    axarr[1,0].plot(lossesx,lossespx,'y.',label=r'$losses:$ ='+str(float(len(lossesx))/len(turns)*100)+'%')
    axarr[1,0].set_xlabel("X [m]")
    axarr[1,0].set_ylabel("PX")
    axarr[1,0].axvline(67.9/1000, color='black', label="ZS.UP")
    
    axarr[2,0].hist(turns, int(c.Nturns/c.turnmultiplicity/10.), normed=1,cumulative=-1, alpha = 0.5, color='g',linewidth=0, label= str(2.*len(turns)/len(os.listdir(lossdir)))+'% extracted')

    
    axarr[2,0].set_xlabel("Turn number of particle extraction")
    axarr[2,0].set_ylabel("1- cumulative particles lost")
    axarr[1,0].legend()
    axarr[2,0].legend()
    
    
    turns=[]    
    x=np.zeros(0)
    px=np.zeros(0)
    
    lossdir="../Data/ripple_t69.0_a126.0_p10k/losses/"
    print "number of batches: " +str(len(os.listdir(lossdir)))
    for fn in os.listdir(lossdir): 
        lossbatch=SingleLossFile(lossdir+fn)
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
    #xglitch = x[np.intersect1d(np.where(turns>5000),np.where(turns<7000))]	
    #pxglitch = px[np.intersect1d(np.where(turns>5000),np.where(turns<7000))]
    #print len(xglitch)
    
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
    axarr[0,1].set_title('glitch_t500e500_a1000')
    axarr[0,1].set_xlabel("Turn number of particle extraction")
    axarr[0,1].set_ylabel("Numer of particles lost")
    
    axarr[1,1].plot(x,px,'.',label=r'$\epsilon$ ='+str(np.sqrt(eps2)))
    #axarr[1,1].plot(xglitch,pxglitch,'r.',alpha=0.5)
    axarr[1,1].plot(lossesx,lossespx,'y.',label=r'$losses:$ ='+str(float(len(lossesx))/len(turns)*100)+'%')
    axarr[1,1].set_xlabel("X [m]")
    axarr[1,1].set_ylabel("PX")
    axarr[1,1].axvline(67.9/1000, color='black', label="ZS.UP")
    
    axarr[2,1].hist(turns, int(c.Nturns/c.turnmultiplicity/10.), normed=1,cumulative=-1, alpha = 0.5, color='g',linewidth=0, label= str(2.*len(turns)/len(os.listdir(lossdir)))+'% extracted')
    axarr[1,1].legend()
    axarr[2,1].legend()
    axarr[2,1].set_xlabel("Turn number of particle extraction")
    axarr[2,1].set_ylabel("1- cumulative particles lost")

#    turns=[]    
#    x=np.zeros(0)
#    px=np.zeros(0)
#    
#    lossdir="../Data/glitch_t500e500_a2000_take3/losses/"
#    print "number of batches: " +str(len(os.listdir(lossdir)))
#    for fn in os.listdir(lossdir): 
#        lossbatch=SingleLossFile(lossdir+fn)
#        turns+=lossbatch["TURN"].tolist()
#        x=np.append(x,lossbatch['X'])
#        px=np.append(px,lossbatch['PX'])
#        
#        boollist = [True if i == 1668.9776 else False for i in lossbatch["S"].tolist()]
#        if(not all(boollist)):
#            print fn
#            
#    #lossesx,lossespx = hitSeptum(x,px)
#    print len(turns)
#    print sum(i > 1000 for i in turns)
#    turns=np.asarray(turns)
#    xglitch = x[np.intersect1d(np.where(turns>5000),np.where(turns<7000))]	
#    pxglitch = px[np.intersect1d(np.where(turns>5000),np.where(turns<7000))]
#    print len(xglitch)
#    #print "percentage of extr hits septum: "+str(float(len(lossesx))/len(turns)*100)
#    xm=np.mean(x)
#    pxm=np.mean(px)
#    
#    dx=(x-xm)
#    dpx=(px-pxm)
#    eps2=np.mean(dx**2)*np.mean(dpx**2)-np.mean((dx)*(dpx))**2
#    print np.sqrt(eps2)
#
#    y,binEdges = np.histogram(turns,bins=int(c.Nturns/c.turnmultiplicity/10.))
#    bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
#    menStd     = np.sqrt(y)
#    axarr[0,2].bar(bincenters, y, width=width,alpha=0.5, color='g', yerr=menStd,linewidth=0)
#    axarr[0,2].set_title('glitch_t500e500_a2000')
#    axarr[0,2].set_xlabel("Turn number of particle extraction")
#    axarr[0,2].set_ylabel("Numer of particles lost")
#    
#    axarr[1,2].plot(x,px,'.',label=r'$\epsilon$ ='+str(np.sqrt(eps2)) )
#    #axarr[1,2].plot(lossesx,lossespx,'y.')
#    axarr[1,2].plot(xglitch,pxglitch,'r.',alpha=0.5)
#    axarr[1,2].set_xlabel("X [m]")
#    axarr[1,2].set_ylabel("PX")
#    axarr[1,2].axvline(67.9/1000, color='black', label="ZS.UP")
#    
#    axarr[2,2].hist(turns, int(c.Nturns/c.turnmultiplicity/10.), normed=1,cumulative=-1, alpha = 0.5, color='g',linewidth=0, label= str(2.*len(turns)/len(os.listdir(lossdir)))+'% extracted')
#    axarr[1,2].legend()
#    axarr[2,2].legend()
#    axarr[2,2].set_xlabel("Turn number of particle extraction")
#    axarr[2,2].set_ylabel("1- cumulative particles lost")
#    
#    
#    
#    turns=[]    
#    x=np.zeros(0)
#    px=np.zeros(0)
#    
#    lossdir="../Data/glitch_t500e1000_a1000_take2/losses/"
#    print "number of batches: " +str(len(os.listdir(lossdir)))
#    for fn in os.listdir(lossdir): 
#        lossbatch=SingleLossFile(lossdir+fn)
#        turns+=lossbatch["TURN"].tolist()
#        x=np.append(x,lossbatch['X'])
#        px=np.append(px,lossbatch['PX'])
#        
#        boollist = [True if i == 1668.9776 else False for i in lossbatch["S"].tolist()]
#        if(not all(boollist)):
#            print fn
#            
#    #lossesx,lossespx = hitSeptum(x,px)
#    print len(turns)
#    print sum(i > 1000 for i in turns)
#    turns=np.asarray(turns)
#    xglitch = x[np.intersect1d(np.where(turns>5000),np.where(turns<7000))]	
#    pxglitch = px[np.intersect1d(np.where(turns>5000),np.where(turns<7000))]
#    print len(xglitch)
#    #print "percentage of extr hits septum: "+str(float(len(lossesx))/len(turns)*100)
#    xm=np.mean(x)
#    pxm=np.mean(px)
#    
#    dx=(x-xm)
#    dpx=(px-pxm)
#    eps2=np.mean(dx**2)*np.mean(dpx**2)-np.mean((dx)*(dpx))**2
#    print np.sqrt(eps2)
#
#    y,binEdges = np.histogram(turns,bins=int(c.Nturns/c.turnmultiplicity/10.))
#    bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
#    menStd     = np.sqrt(y)
#    axarr[0,3].bar(bincenters, y, width=width,alpha=0.5, color='g', yerr=menStd,linewidth=0)
#    axarr[0,3].set_title('glitch_t500e1000_a1000')
#    axarr[0,3].set_xlabel("Turn number of particle extraction")
#    axarr[0,3].set_ylabel("Numer of particles lost")
#    
#    axarr[1,3].plot(x,px,'.',label=r'$\epsilon$ ='+str(np.sqrt(eps2)) )
#    #axarr[1,2].plot(lossesx,lossespx,'y.')
#    axarr[1,3].plot(xglitch,pxglitch,'r.',alpha=0.5)
#    axarr[1,3].set_xlabel("X [m]")
#    axarr[1,3].set_ylabel("PX")
#    axarr[1,3].axvline(67.9/1000, color='black', label="ZS.UP")
#    
#    axarr[2,3].hist(turns, int(c.Nturns/c.turnmultiplicity/10.), normed=1,cumulative=-1, alpha = 0.5, color='g',linewidth=0, label= str(2.*len(turns)/len(os.listdir(lossdir)))+'% extracted')
#    axarr[1,3].legend()
#    axarr[2,3].legend()
#    axarr[2,3].set_xlabel("Turn number of particle extraction")
#    axarr[2,3].set_ylabel("1- cumulative particles lost")
#    plt.tight_layout()
#    plt.show()

losses()
         








































##### ANIMATION

#global dpp, amp, turns 
#dpp, amp, turns = Reader()
#dpp=np.asarray(dpp)
#amp=np.asarray(amp)
#turns=np.asarray(turns)
#tuples = zip(dpp,amp)
#
#    
#    
#def update_plot(i, fig, scat):
#    global dpp, amp, turns 
#    amp=np.where(turns<50*i,0.0009,amp)
#    #amp[turns>10*i ]=0.09
#    tuples = zip(dpp,amp)
#    scat.set_offsets(tuples)
#    print('Frames: %d' %i)
#
#    return scat,
#
##def ScatterAni(dpp, amp, turns):
#    #### ANIMATION
#
#Plotter(dpp, amp, turns)
# 
#fig = plt.figure()
#
#ax = fig.add_subplot(111, autoscale_on=False, xlim=(-0.0003, 0.0003), ylim=(-0.0001, 0.001))
#ax.grid()
#ax.set_title("Slow extraction: 34095 turns, dpp=0.5e-3, 7500 particles")
#ax.set_ylabel(r'Normalised amplitude')
#ax.set_xlabel(r'$\Delta p$')
#
#scat = plt.scatter(dpp, amp)
#anim=animation.FuncAnimation(fig, update_plot, fargs = (fig, scat), frames = 680, interval = 25,repeat_delay=30000, blit=True)
#
#Writer = animation.writers['ffmpeg']
#writer = Writer(fps=15, metadata=dict(artist='Wouter Van De Pontseele'), bitrate=1800)
#
#anim.save('im.mp4', writer=writer)




#### BACKUP


#print  variable_data['BETX']
#print  variable_data['S']
#
#plt.plot(variable_data['S'],variable_data['BETX'])
#plt.plot(variable_data['S'],variable_data['BETY'])
#plt.title('Simulation output')
#plt.xlabel('S')
#plt.ylabel(r'$\beta_x$ and $\beta_y$')