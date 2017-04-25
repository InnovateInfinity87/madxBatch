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

    #if not os.path.exists(trackfolder+'/0'):
    #    subprocess.check_call('./unpacker.sh', shell=True, cwd=trackfolder)

    failed, messages = datproc.errorcheck(errfolder)
    print (str(len(failed))+" processes failed")
    #failed = [int(x) for x in failed]
    #print failed
    #print messages
    datproc.beamstats(lossfolder)
    datproc.wireangle(lossfolder)
    datproc.lossplot(lossfolder, ylim=[-0.0018, -0.0013])
    datproc.lossplot(lossfolder, cax="PT", ylim=[-0.0018, -0.0013])
    datproc.lossplot(lossfolder, ylim=[-0.0016, 0.0016], xax='TURN', yax='PT', cax='PX')
    #datproc.losschart(lossfolder, xax="S", binwidth=0.1)

    #testhead, testtab = datproc.readsingletrack(trackfolder+'/0/track.batch0.obs0001.p0001')

    #print "testhead:\n"
    #print testhead
    #print "\ntesttab:\n"
    #print testtab
