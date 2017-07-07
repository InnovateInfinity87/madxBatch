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
        datproc.efficiency(lossfolder, aperturex=[0.06815,0.08815], aperturey=[-0.023,0.023], zs_len=18.77, zs_an=4.1635E-4)
        print ""
        datproc.wireangle(lossfolder)
        datproc.lossplot(lossfolder, xlim=[0.0675, 0.084], ylim=[-0.0018, -0.0013], cax="PT", save=plotfolder+"/beam.png")
        datproc.lossplot(lossfolder, xlim=[0.0675, 0.070], ylim=[-0.00150, -0.00141], cax="PT", save=plotfolder+"/beam_start.png")
        datproc.lossplot(lossfolder, xlim=[0.0679, 0.0685], ylim=[-0.001465, -0.00145], cax="PT", save=plotfolder+"/beam_start_zoom.png")
        datproc.lossplot(lossfolder, xlim=[0.079, 0.083], ylim=[-0.00175, -0.0017], cax="PT", save=plotfolder+"/beam_end.png")
        datproc.lossplot(lossfolder, xlim=[0.0795, 0.0815], ylim=[-0.00176, -0.0017], cax="PT", save=plotfolder+"/beam_end_zoom.png")
    else:
        datproc.efficiency(lossfolder, pycoll=True, aperturex=[0.04219,0.08219], aperturey=[-0.0100,0.0100])
        print ""
        datproc.lossstats(lossfolder)

        #datproc.zsbacktrack(trackfolder, cax="TURN", save_d=plotfolder+"/zs_down_backtrack_"+name+".png", save_u=plotfolder+"/zs_up_backtrack_"+name+".png")

        datproc.losshistscatter(lossfolder, lossloc=lossloc, xlim=[0.03759, 0.13], ylim=[0,0.004], xbin=0.001, ybin=0.00005, log=True, save=plotfolder+"/tpst_losshist_"+name+".png")
        datproc.losshistscatter(lossfolder, lossloc=lossloc, xax='Y', yax='PY', xlim=[-0.04, 0.04], ylim=[-0.002,0.002], xbin=0.001, ybin=0.00005, log=True, save=plotfolder+"/tpst_losshist_v_.png")
        datproc.losshistscatter(lossfolder, lossloc=lossloc, xax='X', yax='Y', xlim=[0.03759, 0.13], ylim=[-0.04, 0.04], xbin=0.001, ybin=0.001, log=True, save=plotfolder+"/tpst_losshist_s.png")

        datproc.lossplot(lossfolder, lossloc=lossloc, xax='TURN', yax='PT', cax='X', save=plotfolder+"/sweep.png")
        datproc.lossplot(lossfolder, lossloc=None, xlim=[1668,1688], xax='S', yax='X', cax='PT', clim=(-0.0015,0.0015), save=plotfolder+"/zs_loss.png")


    datproc.losshistscatter(lossfolder, lossloc=lossloc, cax="PT", xlim=[0.06815,0.08315], ylim=[-0.0018, -0.00135], xbin=0.0002, ybin=0.00001, log=False, save=plotfolder+"/losshist.png")
    datproc.trackplot(trackfolder, obsloc="obs0002", cax="PT", tpt=300, save=plotfolder+"/zs_up.png")
    datproc.trackplot(trackfolder, obsloc="obs0003", cax="PT", tpt=300, save=plotfolder+"/zs_down.png")
