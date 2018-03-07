# -*- coding: utf-8 -*-
"""
Script to make nominal plots

@author: Linda Stoel
"""
import os
import sys
import subprocess
import matplotlib as mpl
mpl.use('Agg') #To use plotting routines on machines without display
import python.dataprocessing as datproc

tpstcap = 0.03729401
tpsteap = 0.04
tpstblade = 0.0046

def makeplots(folder):
    lossfolder = folder+'/losses'
    trackfolder = folder+'/tracks'
    errfolder = folder+'/error'
    name = folder.split("/")[-1]
    if name=="":
        name = folder.split("/")[-2]

    settings = datproc.getsettings(folder)

    sliced = False if settings['slices']=='None' else True
    pycoll = eval(settings['pycollimate'])
    lossloc = "AP.UP.TPST21760" if pycoll else "AP.UP.ZS21633"

    failed, messages = datproc.errorcheck(errfolder)
    print (str(len(failed))+" processes failed")
    print "errors"
    for message in messages:
        print message
    print ""

    plotfolder = folder+"/../plots/"+name

    if not os.path.exists(plotfolder):
        os.makedirs(plotfolder)

    if not os.path.exists(trackfolder+'/0'):
        subprocess.check_call('./unpacker.sh', shell=True, cwd=trackfolder)

    # Print efficiency and loss stats
    stdout = sys.stdout
    with open(plotfolder+"/stats.txt", 'w') as sys.stdout:
        if pycoll:
            datproc.efficiency(lossfolder, pycoll=True, aperturex=[tpstcap+tpstblade, tpstcap+tpstblade+tpsteap], aperturey=[-0.0100,0.0100])
            print ""
            datproc.lossstats(lossfolder)
            print ""
            datproc.emittance(lossfolder, lossloc=lossloc, ap=[tpstcap+tpstblade, tpstcap+tpstblade+tpsteap], betagamma=426.3167)
        else:
            datproc.efficiency(lossfolder, aperturex=[0.06815,0.08815], aperturey=[-0.023,0.023], zs_len=18.77, zs_an=4.1635E-4, aperturex2=[0.04160,0.06160])
            print ""
            datproc.wireangle(lossfolder)
            print ""
            datproc.emittance(lossfolder, lossloc=lossloc, ap=[0.06815,0.08815], betagamma=426.3167)
    sys.stdout = stdout

    # Make plots
    datproc.lossplot(lossfolder, lossloc=lossloc, xax='TURN', yax='PT', cax='PT', ylim=[-0.0025, 0.0020], clim=[-0.0025, 0.0020], save=plotfolder+"/sweep_"+name+".png")
    datproc.trackplot(trackfolder, obsloc="obs0002", cax="PT", clim=[-0.0025, 0.0020], tpt=1, save=plotfolder+"/zs_up_"+name+".png")
    datproc.trackplot(trackfolder, obsloc="obs0003", cax="PT", clim=[-0.0025, 0.0020], tpt=1, save=plotfolder+"/zs_down_"+name+".png")
    datproc.trackplot(trackfolder, obsloc="obs0004", cax="PT", clim=[-0.0025, 0.0020], tpt=1, save=plotfolder+"/tpst_circ_"+name+".png")

    if sliced:
        datproc.trackplot(trackfolder, obsloc="obs0002", cax="PT", clim=[-0.0025, 0.0020], tpt=300, save=plotfolder+"/zs_up_full_"+name+".png")
        datproc.trackplot(trackfolder, obsloc="obs0003", cax="PT", clim=[-0.0025, 0.0020], tpt=300, save=plotfolder+"/zs_down_full_"+name+".png")
        datproc.trackplot(trackfolder, obsloc="obs0004", cax="PT", clim=[-0.0025, 0.0020], tpt=300, save=plotfolder+"/tpst_circ_full_"+name+".png")

    if pycoll:
        datproc.losshistscatter(lossfolder, lossloc=lossloc, xlim=[tpstcap, 0.095], ylim=[0.0003,0.0027], cax="PT", clim=[-0.0025, 0.0020], xbin=0.001, ybin=0.00005, log=True, save=plotfolder+"/tpst_losshist_h_"+name+".png")
        datproc.losshistscatter(lossfolder, lossloc=lossloc, xax='Y', yax='PY', xlim=[-0.02, 0.02], ylim=[-0.0006,0.0006], cax="PT", clim=[-0.0025, 0.0020], xbin=0.0006, ybin=0.00003, log=True, save=plotfolder+"/tpst_losshist_v_"+name+".png")
        datproc.losshistscatter(lossfolder, lossloc=lossloc, xax='X', yax='Y', xlim=[0.03759, 0.095], ylim=[-0.02, 0.02], cax="PT", clim=[-0.0025, 0.0020], xbin=0.001, ybin=0.0006, log=True, save=plotfolder+"/tpst_losshist_s_"+name+".png")
        datproc.losshistcombo(lossfolder, lossloc=lossloc, xlim=[tpstcap, 0.095], ylim=[0.0003,0.0027], xbin=0.0005, ybin=0.000025, cm='viridis', log=True, save=plotfolder+"/losshist2_"+name+".png")
        datproc.lossplot(lossfolder, lossloc=None, xlim=[1668,1688], xax='S', yax='X', cax="PT", clim=[-0.0025, 0.0020], save=plotfolder+"/zs_loss_"+name+".png")
    else:
        datproc.lossplot(lossfolder, xlim=[0.0675, 0.084], ylim=[-0.0018, -0.0013], cax="PT", clim=[-0.0025, 0.0020], save=plotfolder+"/beam_"+name+".png")
        datproc.losshistscatter(lossfolder, lossloc=lossloc, xlim=[0.06815,0.08515], ylim=[-0.00185, -0.00135],cax="PT", clim=[-0.0025, 0.0020], xbin=0.0002, ybin=0.00001, log=False, save=plotfolder+"/losshist_"+name+".png")
        datproc.losshistcombo(lossfolder, lossloc=lossloc, xlim=[0.06815,0.08515], ylim=[-0.00185, -0.00135], xbin=0.0002, ybin=0.000005, cm='viridis', log=False, save=plotfolder+"/losshist2_"+name+".png")
        datproc.losshistcombo(lossfolder, lossloc=lossloc, xlim=[0.06815,0.07215], ylim=[-0.00155, -0.0014], xbin=0.00008, ybin=0.000004, cm='viridis', log=False, save=plotfolder+"/losshist2_zoom_"+name+".png")


if __name__ == "__main__":
    try:
        study = sys.argv[1]
        makeplots(study)
    except IndexError:
        studies = []
        root = "/afs/cern.ch/project/sloex/nominal/"
        for name in os.listdir(root):
            if os.path.isdir(root+study) and not study=="plots":
                studies += [name]
        for i, study in enumerate(studies):
            print study+' (study '+str(i+1)+' of '+str(len(studies))+')'
            makeplots(root+study)
        

