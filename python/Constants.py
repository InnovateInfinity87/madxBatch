# -*- coding: utf-8 -*-
"""
Created on Tue Jul 26 10:26:33 2016

@author: wvandepo
"""

import os
import sys

### For simulation ###

elements=['AP.UP.ZS21633','ZS.21633','AP.DO.ZS21633']


Nturns=34095
Nbatches=1000
Nparperbatch=50
turnmultiplicity=50 #ffile MADX

createTwiss=False
trackingBool=True
LXplus=True
writetrack=True

ripple=False
startturn=1000
amplitude=1000 #e-8 units k
period=500 #turns
expperiod=0 #turns 0-> pure sine

dataripple=True
ripplefile="ripple"


user = os.environ["USER"]
home = sys.path[0]
twissdir = home+'/twiss/'
inputdir = home+'/input/'
madxdir  = home+'/madx/'

name="name"
#name="t80_a35"

data=home+'Data/'+str(name)+"/"

tracksdir= data+'tracks/'
lossdir  = data+'losses/'
jobsdir  = data+'jobs/'


def SetGeneral(f_Nturns=34095,f_Nbatches=200,f_Nparperbatch=50,f_turnmultiplicity=50):
    global Nturns
    global Nbatches
    global Nparperbatch
    global turnmultiplicity

    Nturns=              f_Nturns
    Nbatches=            f_Nbatches
    Nparperbatch=        f_Nparperbatch
    turnmultiplicity=    f_turnmultiplicity

def SetBools(f_ripple=True,f_dataripple=False,f_writetrack=True):
    global ripple
    global dataripple
    global writetrack

    ripple=f_ripple
    dataripple=f_dataripple
    writetrack = f_writetrack

def SetRipple(f_startturn=1000,f_amplitude=1000,f_period=500,f_expperiod=0):
    global startturn
    global amplitude #e-8 units k
    global period #turns
    global expperiod#turns 0-> pure sine

    startturn=f_startturn
    amplitude=f_amplitude
    period=f_period
    expperiod=f_expperiod

def SetDirs(f_name=None):
    global name
    global data
    global tracksdir
    global lossdir
    global jobsdir

    if f_name is None:
        if( not ripple and not dataripple):
            name='LinSweep_p'+str(int(Nbatches*Nparperbatch/1000))+'k_turns'+str(Nturns)
        elif(expperiod==0):
            name='ripple_t'+str(period)+'_a'+str(amplitude)+'_p'+str(int(Nbatches*Nparperbatch/1000))+'k'
        elif(dataripple):
            print "you need to submit a name for the dataripple dir"
        else:
            name='glitch_t'+str(period)+'e'+str(expperiod)+'_a'+str(amplitude)+'_p'+str(int(Nbatches*Nparperbatch/1000))+'k'
    else:
        name=f_name
    data=home+'Data/'+str(name)+"/"
    tracksdir= data+'tracks/'
    lossdir  = data+'losses/'
    jobsdir  = data+'jobs/'

def setripplefile(f_ripple):
    global ripplefile
    ripplefile=f_ripple
