# -*- coding: utf-8 -*-
"""
Created on Fri Aug 05 14:25:11 2016

@author: wvandepo
"""

import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
from matplotlib import rcParams
from MADXreader import TwissRead,SingleTrack,SingleLossFile
import Constants as c
import os 
import time
import math
import cPickle as pickle
rcParams.update({'font.size': 9})

###############################################################################

plotplace=2 #UP AP ZS
batchrange=3

def TrackRead():
    global batchrange
    global plotplace
    Tracks = {}
    current=-1
    filelist=[os.path.join(dp, f) for dp, dn, fn in os.walk(os.path.expanduser(c.tracksdir)) for f in fn]
    #filelist=os.listdir(c.tracksdir+'0') #only batch zero
    for fn in filelist:
        fnn=fn.split("/")[-1]
        batch =  int(fnn.split(".")[1][5:])
        if(current!=batch and batch <batchrange):
            print 'Processing batch '+str(batch)
            current=batch
        place =  int(fnn.split(".")[2][3:])
        if(place==plotplace and batch <batchrange):
            #print fn
            particle =  int(fnn.split(".")[3][1:])
            Tracks[(batch,particle)]=SingleTrack(fn)
    with open('Tracks','wb') as fp:
        pickle.dump(Tracks,fp)
    return Tracks

def LossRead():
    global batchrange
    bigdictloss={}
    
    for fn in os.listdir(c.lossdir):
          
        batch =  int(fn[5:])
        if(batch<batchrange):
            print fn
            lossbatch=SingleLossFile(c.lossdir+fn)    
            nrinbatch=lossbatch["NUMBER"]-1
            absnr=np.asarray(nrinbatch)+(batch*c.Nparperbatch)
            xexttra=lossbatch["X"]
            lossx=dict(zip(absnr, xexttra))
            bigdictloss.update(lossx)
    print bigdictloss
    with open('Losses','wb') as fp:
        pickle.dump(bigdictloss,fp)
    return lossbatch

        
#initial = np.loadtxt(c.data+"initial_distribution_gauss.txt",skiprows=7,unpack=True)
#with open('initial','wb') as fp:
#        pickle.dump(initial[7],fp)
#
#Losses=LossRead()
#print "losses done"
#Tracks=TrackRead()
#raise SystemExit


###############################################################################

start = time.time()

with open('Tracks','rb') as fp:
    Tracks=pickle.load(fp)
print 'Tracks loaded in RAM'

with open('Losses','rb') as fp:
    Losses=pickle.load(fp)
print 'Loss files loaded in RAM'


with open('initial','rb') as fp:
    inidpp=pickle.load(fp)
print 'Initial particle momenta loaded in RAM'
##############################################################################

print 'Converting tracks to normalised amplitudes'

elements=['"AP.UP.ZS21633"','"ZS.21633"','"AP.DO.ZS21633"']

variable_data, parameters = TwissRead('twiss_after_thinning.prt')
elementindex=variable_data['NAME'].index(elements[plotplace-2])
sqrtbetax = np.sqrt(variable_data['BETX'][elementindex])
alfax = variable_data['ALFX'][elementindex]

#xall =  np.empty(0)
#yall =  np.empty(0)
#pxall = np.empty(0)
#pyall = np.empty(0)
#
#xamp=[]
#
#print 'Calculating average deviation at AP UP ZS'
#for b in range(batchrange):#c.Nbatches):
#    for p in range(c.Nparperbatch):
#        this=Tracks[(b,p+1)]
#        xall = np.append(xall,this['X'])
#        xamp.append(this['X'])
#        pxall = np.append(pxall,this['PX'])
#xall=np.mean(np.nan_to_num(xall))
#pxall=np.mean(np.nan_to_num(pxall))
#print xall,pxall

xall  = 0.0441 
pxall = -0.000969

xlist=[]
amp=[]
dpp=[]
turns=[]
    
for b in range(batchrange):
    for p in range(c.Nparperbatch):
        #print 'batch'+str(b)+" , particle "+str(p)
        this=Tracks[(b,p+1)] 
        if(len(this["TURN"])>0):
            normx=(this['X']-xall)/sqrtbetax
            normpx=normx*alfax + sqrtbetax*(this['PX']-pxall)
            turns.append(this["TURN"][-1])
            xlist.append(this['X'])
            amp.append(np.sqrt(normx**2+normpx**2))
            dpp.append(inidpp[p+b*c.Nparperbatch])
           
print 'Conversion done'

################################################################################
#amp=np.asarray(amp)

maxl = math.ceil(c.Nturns/c.turnmultiplicity)+2
amp=xlist
for p in range(len(dpp)):
    if p in Losses:
        amp[p]=np.append(amp[p],Losses[p])
    print 'particle '+str(p) + " amp len: "+str(len(amp[p]))
    l=maxl-len(amp[p])
    amp[p]=np.lib.pad(amp[p],(0,int(l)),'maximum')
    amp[p]=amp[p].tolist()

#print amp[0]
amp=map(list, zip(*amp))
stepamp=amp[0]
#print len(amp)
#print len(amp[0])

tuples = zip(dpp,stepamp)
  
print 'Plotting started'
    
    
def update_plot(i, fig, scat):
    global dpp, amp, turns 
    stepamp=amp[i]
    #amp=np.where(turns<50*i,0.0009,amp)
    tuples = zip(dpp,stepamp)
    scat.set_offsets(tuples)
    print('Frames: %d' %i)

    return scat,
 
fig = plt.figure()

ax = fig.add_subplot(111, autoscale_on=False, xlim=(-0.0003, 0.0003), ylim=(0.03, 0.085))
ax.axhline(y=0.0679)
ax.grid()
ax.set_title("Slow extraction: 34095 turns, dpp=0.5e-3, 150 particles")
ax.set_ylabel(r'X amplitude')
ax.set_xlabel(r'$\Delta p$')

scat = plt.scatter(dpp, stepamp)
anim=animation.FuncAnimation(fig, update_plot, fargs = (fig, scat), frames = int(math.floor(c.Nturns/c.turnmultiplicity)), interval = 25,repeat_delay=30000, blit=True)

Writer = animation.writers['ffmpeg']
writer = Writer(fps=15, metadata=dict(artist='Wouter Van De Pontseele'), bitrate=1800)

anim.save('im.mp4', writer=writer)


end = time.time()
elapsed = end - start
print "CPU time: "+str(elapsed)+"s"