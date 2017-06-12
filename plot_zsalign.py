# -*- coding: utf-8 -*-
"""
Example of a script to make some standard plots with output from madxBatch

@author: Linda Stoel
"""
import os
import sys
import subprocess
import matplotlib as mpl
mpl.use('Agg') #To use plotting routines on machines without display
import python.dataprocessing as datproc

if __name__ == "__main__":
    folder = sys.argv[1]
    lossfolder = folder+'/losses'
    trackfolder = folder+'/tracks'
    errfolder = folder+'/error'
    name = folder.split("/")[-1]
    if name=="":
        name = folder.split("/")[-2]

    if "_pc" in name:
        pycoll=True
    else:
        pycoll=False

    if pycoll:
        lossloc = "AP.UP.TPST21760"
    else:
        lossloc = "AP.UP.ZS21633"

    zsdo = float(name[3:8])/1E6 + 0.0002

    datproc.efficiency(lossfolder, aperturex=[0.06815,0.08815], aperturey=[-0.023,0.023], zs_len=18.77, zs_an=4.1635E-4, aperturex2=[zsdo,zsdo+0.02])
    print ""
    datproc.lossstats(lossfolder)
