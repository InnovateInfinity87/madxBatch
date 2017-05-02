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

pycoll=False

if __name__ == "__main__":
    folder = sys.argv[1]
    lossfolder = folder+'/losses'
    trackfolder = folder+'/tracks'
    errfolder = folder+'/error'
    plotfolder = folder+"/../plots"#"/../sweepspeed_plots"#
    name = folder.split("/")[-1]
    if name=="":
        name = folder.split("/")[-2]

    lossloc = "AP.UP.ZS21633"

    if not os.path.exists(trackfolder+'/0'):
        subprocess.check_call('./unpacker.sh', shell=True, cwd=trackfolder)

    failed, messages = datproc.errorcheck(errfolder)
    print (str(len(failed))+" processes failed")
    #failed = [int(x) for x in failed]
    #print failed
    print messages
  #  print ""
  #  datproc.beamstats(lossfolder, lossloc=lossloc)
    print ""
    datproc.efficiency(lossfolder, aperturex=[0.06815,0.08815], aperturey=[-0.023,0.023])
    print ""
    datproc.lossstats(lossfolder)
    print ""
    datproc.wireangle(lossfolder)
    datproc.lossplot(lossfolder, save=plotfolder+"/beam_"+name+".png")
    datproc.lossplot(lossfolder, xax='TURN', yax='PT', cax='PX', save=plotfolder+"/sweep_"+name+".png")
    datproc.trackplot(trackfolder, obsloc="obs0002", cax="TURN", tpt=9, save=plotfolder+"/zs_tracks_"+name+".png")

    datproc.losshistscatter(lossfolder, lossloc=lossloc, xlim=[0.06815,0.08815], ylim=[-1.2E-3, -1.8E-3], xbin=0.001, ybin=0.00005, log=False, save=plotfolder+"/losshist_"+name+".png")
   
