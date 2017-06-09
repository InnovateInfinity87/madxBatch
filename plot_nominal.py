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

    if "pc" in name and "nopc" not in name:
        pycoll=True
    else:
        pycoll=False

    if pycoll:
        lossloc = "AP.UP.TPST21760"
    else:
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
    if not pycoll:
        datproc.efficiency(lossfolder, aperturex=[0.06815,0.08815], aperturey=[-0.023,0.023], zs_len=18.77, zs_an=4.1635E-4, aperturex2=[0.04265,0.06265])
        print ""
        datproc.wireangle(lossfolder)
        datproc.lossplot(lossfolder, ylim=[-0.0018, -0.0013], save=plotfolder+"/beam.png")
        datproc.lossplot(lossfolder, ylim=[-0.0016, 0.0016], xax='TURN', yax='PT', cax='PX', save=plotfolder+"/sweep.png")
    else:
        datproc.efficiency(lossfolder, pycoll=True, aperturex=[0.04219,0.08219], aperturey=[-0.0100,0.0100])
        print ""
        datproc.lossstats(lossfolder)

        #datproc.zsbacktrack(trackfolder, cax="TURN", save_d=plotfolder+"/zs_down_backtrack_"+name+".png", save_u=plotfolder+"/zs_up_backtrack_"+name+".png")

        datproc.trackplot(trackfolder, obsloc="obs0002", cax="TURN", tpt=1, save=plotfolder+"/zs_up_"+name+".png")

        datproc.trackplot(trackfolder, obsloc="obs0003", cax="TURN", tpt=1, save=plotfolder+"/zs_down_"+name+".png")

        datproc.losshistscatter(lossfolder, lossloc=lossloc, xlim=[0.03759, 0.13], ylim=[0,0.004], xbin=0.001, ybin=0.00005, log=True, save=plotfolder+"/tpst_losshist_"+name+".png")
        datproc.losshistscatter(lossfolder, lossloc=lossloc, xax='Y', yax='PY', xlim=[-0.04, 0.04], ylim=[-0.002,0.002], xbin=0.001, ybin=0.00005, log=True, save=plotfolder+"/tpst_losshist_v_.png")
        datproc.losshistscatter(lossfolder, lossloc=lossloc, xax='X', yax='Y', xlim=[0.03759, 0.13], ylim=[-0.04, 0.04], xbin=0.001, ybin=0.001, log=True, save=plotfolder+"/tpst_losshist_s.png")

        datproc.lossplot(lossfolder, lossloc="AP.DO.ZS21676_M", xax='X', yax='Y', cax='PT', save=plotfolder+"/zs_do_loss.png")

        datproc.lossplot(lossfolder, lossloc=lossloc, xax='TURN', yax='PT', cax='X', save=plotfolder+"/sweep.png")
        datproc.lossplot(lossfolder, lossloc=None, xlim=[1668,1688], xax='S', yax='X', cax='PT', clim=(-0.0015,0.0015), save=plotfolder+"/zs_loss.png")
