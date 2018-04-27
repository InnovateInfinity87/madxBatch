# -*- coding: utf-8 -*-
"""
Track madxBatch output to another point in the machine (<1 turn)

WARNING: Will not apply correct bumps for sliced simulations
WARNING: Backtrack does not actually work.

@author: Linda Stoel
"""
from version import __version__

import os
import sys
import subprocess

locations = {"handover": ["QDA->L/2+0.004724617306961", "QDA.21910"],
             "zsup": ["0", "AP.UP.ZS21633"], "tpstup": ["0", "AP.UP.TPST21760"],
             "Q216start": ["-QFA->L/2", "QFA.21610"],
             "Q216": ["0", "QFA.21610"], "Q219": ["0", "QDA.21910"],
             "grid218": ["0", "AP.UP.MSE21857"]}

def findsloexcodedir(madcode):
    for line in madcode:
        if "madxBatch" in line:
            if line.startswith("CALL, FILE = '"):
                codedir = line[14:-2].split("/")
                i = codedir.index("madxBatch")
                codedir = "/".join(codedir[:i])
                return codedir
    raise IOError('Could not locate slow extraction code directory...')
    
def installmarkers(startlocation, endlocation, twisstrack):
    startspec = locations[startlocation]
    endspec = locations[endlocation]
    if twisstrack:
        ap = ""
    else:
       ap = " , APERTYPE=CIRCLE, APERTURE=1E-9"
    madcode = ("trackfrom_m: MARKER;\n"+
               "trackto_m: MARKER"+ap+";\n\n"+

               "SEQEDIT, SEQUENCE=sps;\n"+
               " FLATTEN;\n"+
               " INSTALL, ELEMENT=trackfrom_m, AT="+startspec[0]+", FROM="+startspec[1]+";\n"+
               " INSTALL, ELEMENT=trackto_m, AT="+endspec[0]+", FROM="+endspec[1]+";\n"+
               " FLATTEN;\n"+
               "ENDEDIT;\n\n"+

               "USE, SEQUENCE = sps;\n\n")
    
    return madcode

def cutsequence():
    madcode = ("SEQEDIT, SEQUENCE=sps;\n"+
               " FLATTEN;\n"+
               " CYCLE, START=trackfrom_m;\n"+
               " FLATTEN;\n"+
               " EXTRACT, SEQUENCE=sps, FROM=trackfrom_m, TO=trackto_m, NEWNAME=trackseq;\n"+
               " FLATTEN;\n"+
               "ENDEDIT;\n\n"+
               
               "USE, SEQUENCE = trackseq;\n\n")
    
    return madcode


def tracklossto(location, lossfolder, startpoint, backtrack=False, sloexcodedir=None, madxexe=None,
                twisstrack=False, static=False, batch=0):
    '''Track particles in lossfolder to location, save in lossfolder/../twisstrackloss_startpoint_to_location'''
    if not location in locations:
        locations['custom'] = ["0", location]
        location = 'custom'

    if twisstrack:
        savefolder = lossfolder+"/../twisstrackloss_"+startpoint+"_to_"+location+"/"+str(batch)
    else:
        savefolder = lossfolder+"/../trackloss_"+startpoint+"_to_"+location
    jobfile = lossfolder+"/../jobs/"+str(batch)+".madx"
    madfile = savefolder+"/track.madx"
    hastrmacro = False

    if not os.path.isfile(jobfile):
        print "Could not find expected jobfile at "+jobfile+"\n"
        return

    if not os.path.exists(savefolder):
        os.makedirs(savefolder)
    os.chdir(savefolder)

    # Identify and re-use code to build the initial sequence
    # (Expects default format of madxBatch jobs)
    with open(jobfile, 'r') as infile:
        madstart = infile.readlines()
    for i, line in enumerate(madstart):
        if static:
            if line.startswith("! Tracking code"):
                break
        else:
            if line.startswith("TRACK,"):
               break
    madstart = madstart[:i-1]

    try:
        madstart.remove("SYSTEM, 'mkdir "+str(batch)+"';\n")
    except ValueError:
        pass
    
    for i, line in enumerate(madstart):
        if "tr$macro" in line:
            hastrmacro = True
            madstart[i] = line.replace("tr$macro", "turnmacro")

    if sloexcodedir is None:    
        sloexcodedir = findsloexcodedir(madstart)

    # Construct the new MAD-X code
    with open(madfile, 'w') as outf:
        outf.writelines(madstart)
        
        outf.write("SELECT, FLAG=ERROR, FULL;\n"+
                      "ESAVE, FILE=errfile;\n\n")
                      
        outf.write(installmarkers(startpoint, location, twisstrack))

        outf.write("READMYTABLE, file=errfile, table=errtab;\n"+
                        "SETERR, TABLE=errtab;\n\n")

        if twisstrack:
            outf.write("SAVEBETA, LABEL=starttrack, PLACE=trackfrom_m, SEQUENCE=sps;\n")
        
        outf.write("TWISS;\n"+
                      "trackto_start_s = TABLE(TWISS, trackfrom_m, S);\n\n")

        if twisstrack:
            outf.write("lss2_noapp = 1;\n")
        outf.write("CALL, FILE='"+sloexcodedir
                      +"/madxBatch/madx/lss2extraction.cmdx';\n\n")

        if twisstrack:
            outf.write("SELECT, FLAG=ERROR, CLEAR;\n")
            outf.write("SELECT, FLAG=ERROR, RANGE=trackto_m;\n")
            outf.write("EALIGN, AREX=-5, AREY=-5;\n\n")
        
        outf.write("SELECT, FLAG=ERROR, FULL;\n"+
                      "ESAVE, FILE=errfile2;\n\n")

        outf.write(cutsequence())

        outf.write("READMYTABLE, file=errfile2, table=errtab;\n"+
                    "SETERR, TABLE=errtab;\n\n")

        outf.write("hastrmacro = "+("1" if hastrmacro else "0")+";\n")
        outf.write("backtrack = "+("1" if backtrack else "0")+";\n\n")

        outf.write("CALL, FILE='"+sloexcodedir
                   +"/madxBatch/madx/setseptamacro.cmdx';\n\n")

        if twisstrack:
            outf.write("CALL, FILE='"+sloexcodedir
                          +"/madxBatch/madx/twisstracktomacro.cmdx';\n\n")

            outf.write("PTOT_GEV := (BEAM->PC)/(1-PT) + TABLE(PTC_TWISS, X)*0.0;\n"+
                       "X_CM := TABLE(PTC_TWISS, X)*100;\n"+
                       "Y_CM := TABLE(PTC_TWISS, Y)*100;\n"+
                       "Z_CM := TABLE(PTC_TWISS, S)*100;\n"+
                       "COSX := TABLE(PTC_TWISS, PX)*(1-PT);\n"+
                       "COSY := TABLE(PTC_TWISS, PY)*(1-PT);\n\n")

            outf.write("SELECT, FLAG=ptc_twiss, CLEAR;\n"+
                       "SELECT, FLAG=ptc_twiss, COLUMN=NAME, PTOT_GEV, X_CM, Y_CM, Z_CM, COSX, COSY, S, X, PX, Y, PY, T, PT, TURN;\n\n")

            outf.write(" EXEC, trackto('"+lossfolder+"/"+str(batch)+".tfs', "+str(batch)+");\n")
        else:
            outf.write("CALL, FILE='"+sloexcodedir
                          +"/madxBatch/madx/tracktomacro.cmdx';\n\n")
            outf.write("EXEC, trackto('"+lossfolder+"/"+str(batch)+".tfs', '"+str(batch)+".tfs');\n\n")
        
        outf.flush()
        
    if madxexe is None:
        madxexe = "/afs/cern.ch/user/m/mad/bin/madx"
    with open(os.devnull, 'w') as quiet:
        subprocess.check_call(madxexe+" track.madx", shell=True, stdout=quiet)
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
