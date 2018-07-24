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

    plotfolder = folder+"/../plots/"+name

    if not os.path.exists(plotfolder):
        os.makedirs(plotfolder)

    lossloc = "AP.UP.ZS21633"


    for width in [200,500,1000,2000]:
        print width
        datproc.wireangle(lossfolder, wiremax=(0.06795+width*1E-6))
        datproc.losshistscatter(lossfolder, lossloc=lossloc, cax="PT",
                                xlim=[0.06795,0.06995], ylim=[-0.0016, -0.00135],
                                xbin=0.00005, ybin=0.000005,
                                datalim=[[None,0.06795+width*1E-6],[None,None],[None,None]],
                                log=False, save=plotfolder+"/losshist_"+str(width)+".png")

    datproc.efficiency(lossfolder, aperturex=[0.06815,0.08815], aperturey=[-0.023,0.023], zs_len=18.77, zs_an=4.1635E-4, aperturex2=[0.04265,0.06265])
    datproc.lossstats(lossfolder)
