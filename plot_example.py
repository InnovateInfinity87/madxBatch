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

    if not os.path.exists(trackfolder+'/0'):
        subprocess.check_call('./unpacker.sh', shell=True, cwd=trackfolder)

    datproc.lossplot(lossfolder)
    datproc.losschart(lossfolder, xax="S", binwidth=1)

    testhead, testtab = datproc.readsingletrack(trackfolder+'/0/track.batch0.obs0001.p0001')

    print "testhead:\n"
    print testhead
    print "\ntesttab:\n"
    print testtab
