# -*- coding: utf-8 -*-
"""
Example of a script to make some standard plots with output from madxBatch

@author: Linda Stoel
"""
import os
import sys
import subprocess
import python.dataprocessing as datproc

pycol=False

if __name__ == "__main__":
    folder = sys.argv[1]
    lossfolder = folder+'/losses'
    trackfolder = folder+'/tracks'
    errfolder = folder+'/error'
    plotfolder = folder+"/../sweepspeed_plots"#"/../plots"#
    name = folder.split("/")[-1]
    if name=="":
        name = folder.split("/")[-2]

    if pycol:
        lossloc = "AP.UP.TPST21760"
    else:
        lossloc = "AP.UP.ZS21633"

    #if not os.path.exists(trackfolder+'/0'):
    #    subprocess.check_call('./unpacker.sh', shell=True, cwd=trackfolder)

    failed, messages = datproc.errorcheck(errfolder)
    print (str(len(failed))+" processes failed")
    #failed = [int(x) for x in failed]
    #print failed
    #print messages
    #datproc.beamstats(lossfolder, lossloc=lossloc)
    if not pycol:
        #datproc.wireangle(lossfolder)
        #datproc.lossplot(lossfolder, ylim=[-0.0018, -0.0013], save=plotfolder+"/beamt_"+name+".png")
        datproc.lossplot(lossfolder, cax="PT", ylim=[-0.0018, -0.0013], save=plotfolder+"/beam_"+name+".png")
        datproc.lossplot(lossfolder, ylim=[-0.0016, 0.0016], xax='TURN', yax='PT', cax='PX', save=plotfolder+"/sweep_"+name+".png")
    else:
        datproc.lossplot(lossfolder, lossloc=None, xlim=[1668,1688], xax='S', yax='X', cax='PT', save=plotfolder+"/check_"+name+".png")
        datproc.lossplot(lossfolder, lossloc=lossloc, cax="PT", save=plotfolder+"/fullbeam_"+name+".png")
        datproc.lossplot(lossfolder, lossloc=lossloc, xlim=[0.03759, 0.12], cax="PT", save=plotfolder+"/exbeam_"+name+".png")
        datproc.lossplot(lossfolder, lossloc=lossloc, xax='TURN', yax='PT', cax='PX', save=plotfolder+"/sweep_"+name+".png")
        datproc.losschart(lossfolder, xax="S", binwidth=0.1, save=plotfolder+"/lossmap_"+name+".png")

    #testhead, testtab = datproc.readsingletrack(trackfolder+'/0/track.batch0.obs0001.p0001')

    #print "testhead:\n"
    #print testhead
    #print "\ntesttab:\n"
    #print testtab
