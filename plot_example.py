# -*- coding: utf-8 -*-
"""
Example of a script to make some standard plots with output from madxBatch

@author: Linda Stoel
"""
import os
import sys
import subprocess
import python.dataprocessing as datproc

if __name__ == "__main__":
    folder = sys.argv[1]
    lossfolder = folder+'/losses'
    trackfolder = folder+'/tracks'
    errfolder = folder+'/error'
    plotfolder = folder+"/../plots"
    name = folder.split("/")[-1]
    if name=="":
        name = folder.split("/")[-2]

    #if not os.path.exists(trackfolder+'/0'):
    #    subprocess.check_call('./unpacker.sh', shell=True, cwd=trackfolder)

    failed, messages, missout = datproc.errorcheck(errfolder)
    print (str(len(failed))+" processes failed, among successfull jobs "+str(len(missout))+" output files missing")
    #failed = [int(x) for x in failed]
    #print failed
    #print messages
    datproc.beamstats(lossfolder)
    datproc.wireangle(lossfolder)
    datproc.lossplot(lossfolder, ylim=[-0.0018, -0.0013], save=plotfolder+"/turnnum_"+name+".png")
    datproc.lossplot(lossfolder, cax="PT", ylim=[-0.0018, -0.0013], save=plotfolder+"/mom_"+name+".png")
    datproc.lossplot(lossfolder, ylim=[-0.0016, 0.0016], xax='TURN', yax='PT', cax='PX', save=plotfolder+"/turnptpx_"+name+".png")
    datproc.losschart(lossfolder, xax="S", binwidth=0.1, save=plotfolder+"/lossmap_"+name+".png")

    #testhead, testtab = datproc.readsingletrack(trackfolder+'/0/track.batch0.obs0001.p0001')

    #print "testhead:\n"
    #print testhead
    #print "\ntesttab:\n"
    #print testtab
