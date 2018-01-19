# -*- coding: utf-8 -*-
"""
Batched MAD-X slow extraction

Re-submit jobs that failed (obviously only use this if they failed due
to the cluster, not if there's a mistake in your code...) to be used
from your terminal as:
python resub.py foldername
where you have checked that _all_ jobs in that folder have compleded,
and some had errors. (Will fail at >11 resubs.)

(This script was made to solve black holing issues in madxBatch and may
not apply in other cases)

@author: Linda Stoel
"""
import os
import sys
import fileinput
import subprocess
from shutil import copyfile
from dataprocessing import errorcheck

if __name__ == "__main__":
    folder = sys.argv[1]
    subfiles = []

    for file_ in os.listdir(folder):
        if file_.endswith(".sub"):
            subfiles += [file_]
    subcount = len(subfiles)

    if subcount==0:
        print "no submit file found"
        sys.exit()
    elif subcount==1:
        toclone = subfiles[0]
        resub = folder+"/"+toclone[:-4]+"_r"+str(subcount)+".sub"
    else:
        for sub in subfiles:
            if sub.endswith(str(subcount-1)+".sub"):
                toclone = sub
                resub = folder+"/"+toclone[:-5]+str(subcount)+".sub"

    failed, _ = errorcheck(folder+"/error", checkerr=False)

    if len(failed)==0:
        print "no failed jobs found"
        sys.exit()

    with open(resub, 'w') as out:
        for line in fileinput.input(folder+"/"+toclone, inplace=0):
            line = line.replace("ProcId", "ReId")
            out.write(line)

    for line in fileinput.input(resub, inplace=1):
        if line.startswith("log"):
            line = "log = log_r"+str(subcount)+".txt\n"
        elif line.startswith("queue"):
            line = "queue ReId in "+", ".join(failed)+"\n"
        sys.stdout.write(line)

    subprocess.check_call("condor_submit "+resub, shell=True, cwd=folder)
