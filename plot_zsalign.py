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

thickness = 100E-6

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

    zsup = 0.06795 + thickness
    zsdo = float(name.split("_")[1])/1E6 + thickness

    if pycoll:
        datproc.efficiency(lossfolder, pycoll=True, aperturex=[0.04219,0.08219], aperturey=[-0.0100,0.0100])
    else:
        datproc.efficiency(lossfolder, aperturex=[zsup, zsup+0.02], aperturey=[-0.023,0.023], zs_len=18.77, zs_an=4.1635E-4, aperturex2=[zsdo,zsdo+0.02])
    print ""
    datproc.lossstats(lossfolder)
