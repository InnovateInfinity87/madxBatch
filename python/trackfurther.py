# -*- coding: utf-8 -*-
"""
Track madxBatch output to another point in the machine (<1 turn)

WARNING: Will not apply correct bumps for sliced simulations
WARNING: Backtrack does not actually work.

@author: Linda Stoel
"""
import os
import sys
import subprocess

locations = {"handover": ["QDA->L/2+0.00456325", "QDA.21910"],
             "zsup": ["0", "AP.UP.ZS21633"], "tpstup": ["0", "AP.UP.TPST21760"],
             "Q216start": ["-QFA->L/2", "QFA.21610"],
             "Q216": ["0", "QFA.21610"], "Q219": ["0", "QDA.21910"]}

def findsloexcodedir(madcode):
    for line in madcode:
        if "madxBatch" in line:
            if line.startswith("CALL, FILE = '"):
                codedir = line[14:-2].split("/")
                i = codedir.index("madxBatch")
                codedir = "/".join(codedir[:i-1])
                return codedir
    raise IOError('Could not locate slow extraction code directory...')
    
def installmarkers(startlocation, endlocation):
    startspec = locations[startlocation]
    endspec = locations[endlocation]
    madcode = ("trackfrom_m: MARKER;\n"+
               "trackto_m: MARKER, APERTYPE=circle, APERTURE={7.0E-5};\n\n"+
                    
               "SEQEDIT, SEQUENCE=sps;\n"+
               " FLATTEN;\n"+
               " INSTALL, ELEMENT=trackfrom_m, AT="+startspec[0]+", FROM="+startspec[1]+";\n"+
               " INSTALL, ELEMENT=trackto_m, AT="+endspec[0]+", FROM="+endspec[1]+";\n"+
               " FLATTEN;\n"+
               "ENDEDIT;\n\n"+
               
               "USE, SEQUENCE = sps;\n\n")
    
    return madcode


def tracklossto(location, lossfolder, startpoint, backtrack=False, sloexcodedir=None, madxexe=None,
                twisstrack=False):
    "Track particles in lossfolder to location, save in ../location"
    if not location in locations:
        locations['custom'] = ["0", location]
        location = 'custom'

    if twisstrack:
        savefolder = lossfolder+"/../twisstrackloss_"+startpoint+"_to_"+location
    else:
        savefolder = lossfolder+"/../trackloss_"+startpoint+"_to_"+location
    jobfile = lossfolder+"/../jobs/0.madx"
    madfile = savefolder+"/track.madx"
    hastrmacro = False

    if not os.path.isfile(jobfile):
        print "Could not find expected jobfile at "+jobfile+"\n"
        return

    if not os.path.exists(savefolder):
        os.mkdir(savefolder)
    os.chdir(savefolder)

    # Identify and re-use code to build the initial sequence
    # (Expects default format of madxBatch jobs)
    with open(jobfile, 'r') as infile:
        madstart = infile.readlines()
    for i, line in enumerate(madstart):
        if line.startswith("TRACK,"):
           break
    madstart = madstart[:i-1]
    madstart.remove("SYSTEM, 'mkdir 0';\n")
    
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
                      
        outf.write(installmarkers(startpoint, location))

        if twisstrack:
            outf.write("READMYTABLE, file=errfile, table=errtab;\n"+
                        "SETERR, TABLE=errtab;\n"+
                        "SAVEBETA, LABEL=starttrack, PLACE=trackfrom_m, SEQUENCE=sps;\n\n")
        
        outf.write("TWISS;\n"+
                      "trackto_start_s = TABLE(TWISS, trackfrom_m, S);\n\n")
                      
        outf.write("hastrmacro = "+("1" if hastrmacro else "0")+";\n")
        outf.write("backtrack = "+("1" if backtrack else "0")+";\n\n")
        
        outf.write("CALL, FILE='"+sloexcodedir
                      +"/madxBatch/madx/lss2extraction.cmdx';\n")
        if twisstrack:
            outf.write("CALL, FILE='"+sloexcodedir
                          +"/madxBatch/madx/twisstracktomacro.cmdx';\n\n")
            outf.write("SELECT, FLAG=twiss, CLEAR;\n"+
                       "SELECT, FLAG=twiss, COLUMN=NAME, S, X, PX, Y, PY, T, PT, TURN;\n\n")
            for lossfile in os.listdir(lossfolder):
                outf.write("EXEC, trackto('"+lossfolder+"/"+lossfile+"', "+lossfile.split(".")[0]+");\n\n")
        else:
            outf.write("CALL, FILE='"+sloexcodedir
                          +"/madxBatch/madx/tracktomacro.cmdx';\n\n")
            for lossfile in os.listdir(lossfolder):
                outf.write("EXEC, trackto('"+lossfolder+"/"+lossfile+"', '"+lossfile+"');\n\n")
        
        outf.flush()
        
    if madxexe is None:
        madxexe = "madx_dev"
    subprocess.check_call(madxexe+" track.madx", shell=True, stdout=sys.stdout)
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
