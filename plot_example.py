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

pycoll=True

if __name__ == "__main__":
    folder = sys.argv[1]
    lossfolder = folder+'/losses'
    trackfolder = folder+'/tracks'
    errfolder = folder+'/error'
    plotfolder = folder+"/../plots"#"/../sweepspeed_plots"#
    name = folder.split("/")[-1]
    if name=="":
        name = folder.split("/")[-2]

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
        datproc.efficiency(lossfolder, aperturex=[0.06815,0.08815], aperturey=[-0.023,0.023])
        print ""
        datproc.wireangle(lossfolder)
        datproc.lossplot(lossfolder, ylim=[-0.0018, -0.0013], save=plotfolder+"/beam_"+name+".png")
        #datproc.lossplot(lossfolder, cax="PT", ylim=[-0.0018, -0.0013], save=plotfolder+"/beam_pt_"+name+".png")
        datproc.lossplot(lossfolder, ylim=[-0.0016, 0.0016], xax='TURN', yax='PT', cax='PX', save=plotfolder+"/sweep_"+name+".png")
    else:
        datproc.efficiency(lossfolder, pycoll=True, aperturex=[0.04219,0.08219], aperturey=[-0.0100,0.0100])
        print ""
        datproc.lossstats(lossfolder)

        #datproc.fullplot(folder, lossloc="TCE.21695", obsloc="obs0005", tpt=1, save=plotfolder+"/tce_beam_"+name+".png")
        #datproc.fullplot(folder, lossloc="TCE.21695", obsloc="obs0005", tpt=1, xax='Y', yax='PY', save=plotfolder+"/tce_beam_v_"+name+".png")
        #datproc.fullplot(folder, lossloc="TCE.21695", obsloc="obs0005", tpt=1, xax='X', yax='Y', save=plotfolder+"/tce_beam_s_"+name+".png")
        #datproc.zsbacktrack(trackfolder, cax="TURN", save_d=plotfolder+"/zs_down_backtrack_"+name+".png", save_u=plotfolder+"/zs_up_backtrack_"+name+".png")

       # datproc.trackplot(trackfolder, obsloc="obs0002", cax="TURN", tpt=1, save=plotfolder+"/zs_up_"+name+".png")
       # datproc.trackplot(trackfolder, obsloc="obs0002", cax="TURN", tpt=9, save=plotfolder+"/zs_up_more_"+name+".png")
        #datproc.trackplot(trackfolder, obsloc="obs0002", xax='Y', yax='PY', cax="TURN", tpt=1, save=plotfolder+"/zs_up_v_"+name+".png")
        #datproc.trackplot(trackfolder, obsloc="obs0002", xax='X', yax='Y', cax="TURN", tpt=1, save=plotfolder+"/zs_up_s_"+name+".png")
        #datproc.trackplot(trackfolder, obsloc="obs0002", cax="TURN", tpt=None, batches=[0,1,2,3,4], save=plotfolder+"/full_zs_up_"+name+".png")

        #datproc.trackplot(trackfolder, obsloc="obs0002", cax="TURN", tpt=50, save=plotfolder+"/zs_up_big_"+name+".png")

        #datproc.trackplot(trackfolder, obsloc="obs0003", cax="TURN", tpt=1, save=plotfolder+"/zs_down_"+name+".png")
        #datproc.trackplot(trackfolder, obsloc="obs0003", xax='Y', yax='PY', cax="TURN", tpt=1, save=plotfolder+"/zs_down_v_"+name+".png")
        #datproc.trackplot(trackfolder, obsloc="obs0003", xax='X', yax='Y', cax="TURN", tpt=1, save=plotfolder+"/zs_down_s_"+name+".png")
        #datproc.trackplot(trackfolder, obsloc="obs0003", cax="TURN", tpt=None, batches=[0,1,2,3,4], save=plotfolder+"/full_zs_down_"+name+".png")

        #datproc.lossplot(lossfolder, lossloc="TCE.21695", cax="TURN", save=plotfolder+"/tce_loss_"+name+".png")
        #datproc.lossplot(lossfolder, lossloc="TCE.21695", xax='Y', yax='PY', cax="TURN", save=plotfolder+"/tce_loss_v_"+name+".png")
        #datproc.lossplot(lossfolder, lossloc="TCE.21695", xax='X', yax='Y', cax="TURN", save=plotfolder+"/tce_loss_s_"+name+".png")

       # datproc.fullplot(folder, lossloc=lossloc, obsloc="obs0004", tpt=1, save=plotfolder+"/tpst_beam_"+name+".png")
       # datproc.fullplot(folder, lossloc=lossloc, obsloc="obs0004", xlim=[0.01, 0.08], ylim=[-0.0015,0.0015], tpt=1, save=plotfolder+"/tpst_cut_"+name+".png")
        #datproc.fullplot(folder, lossloc=lossloc, obsloc="obs0004", tpt=1, xax='Y', yax='PY', save=plotfolder+"/tpst_beam_v_"+name+".png")
        #datproc.fullplot(folder, lossloc=lossloc, obsloc="obs0004", tpt=1, xax='X', yax='Y', save=plotfolder+"/tpst_beam_s_"+name+".png")
       # datproc.lossplot(lossfolder, lossloc=lossloc, cax="TURN", save=plotfolder+"/tpst_loss_"+name+".png")
       # datproc.lossplot(lossfolder, lossloc=lossloc, xax='Y', yax='PY', cax="TURN", save=plotfolder+"/tpst_loss_v_"+name+".png")
       # datproc.lossplot(lossfolder, lossloc=lossloc, xax='X', yax='Y', cax="TURN", save=plotfolder+"/tpst_loss_s_"+name+".png")
       # datproc.lossplot(lossfolder, lossloc=lossloc, xlim=[0.03759, 0.13], ylim=[0,0.004], cax="TURN", save=plotfolder+"/tpst_exbeam_"+name+".png")
       # datproc.lossplot(lossfolder, lossloc=lossloc, xlim=[0.03759, 0.07759+0.0046], ylim=[0.0004,0.0015], cax="TURN", save=plotfolder+"/tpst_exbeam2_"+name+".png")
        #datproc.lossplot(lossfolder, lossloc=lossloc, xlim=[0.03759, 0.13], ylim=[0,0.004], cax="PT", save=plotfolder+"/tpst_exbeam_pt_"+name+".png")
       # datproc.losshistscatter(lossfolder, lossloc=lossloc, xlim=[0.03759, 0.13], ylim=[0,0.004], xbin=0.001, ybin=0.00005, log=True, save=plotfolder+"/tpst_losshist_"+name+".png")
       # datproc.losshistscatter(lossfolder, lossloc=lossloc, xax='Y', yax='PY', xlim=[-0.04, 0.04], ylim=[-0.002,0.002], xbin=0.001, ybin=0.00005, log=True, save=plotfolder+"/tpst_losshist_v_"+name+".png")
       # datproc.losshistscatter(lossfolder, lossloc=lossloc, xax='X', yax='Y', xlim=[0.03759, 0.13], ylim=[-0.04, 0.04], xbin=0.001, ybin=0.001, log=True, save=plotfolder+"/tpst_losshist_s_"+name+".png")

        #datproc.lossplot(lossfolder, lossloc='AP.UP.MSE21832', cax="TURN", save=plotfolder+"/mseloss_"+name+".png")

       # datproc.lossplot(lossfolder, lossloc=lossloc, xax='TURN', yax='PT', cax='X', save=plotfolder+"/sweep_"+name+".png")
       # datproc.lossplot(lossfolder, lossloc=None, xlim=[1668,1688], xax='S', yax='X', cax='PT', clim=(-0.0015,0.0015), save=plotfolder+"/zs_loss_"+name+".png")
       # datproc.losschart(lossfolder, xax="S", binwidth=0.1, save=plotfolder+"/lossmap_"+name+".png")

    #testhead, testtab = datproc.readsingletrack(trackfolder+'/0/track.batch0.obs0001.p0001')

    #print "testhead:\n"
    #print testhead
    #print "\ntesttab:\n"
    #print testtab
