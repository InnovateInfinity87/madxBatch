# -*- coding: utf-8 -*-
"""
Batched MAD-X slow extraction

Modifies files for running MAD-X simulations of slow extraction with a
tune ripple and then submits them to the CERN LSF batching service.

@author: Wouter van de Pontseele, Linda Stoel
"""
import fileinput
import subprocess
import python.make_ps_distribution as dis
import os
import sys
import numpy as np
import python.Constants as c

def whichQueue():
    """Determines appropriate LSF queue based on particles*turns."""
    n = c.Nturns*c.Nparperbatch
    queue = "1nh"
    if n > 10e6:
        queue = "8nh"
    if n > 10e6:
        queue = "1nd"
    elif n > 2e7:
        queue = "2nd"
    elif n > 7e7:
        queue = "1nw"
    elif n > 1e8:
        queue = "2nw"
    return queue


def replacer(filename,searchExp,replaceExp):
    """In-place replacement in text files.

    Finds all occurences of searchExp in the file and replaces each one
    with replaceExp, overwriting the original file.
    """
    for line in fileinput.input(filename, inplace=1):
        if searchExp in line:
            line = line.replace(searchExp,replaceExp)
        sys.stdout.write(line)

def adapter(filename,searchExp,replaceExp):
    """Creates a new file adapted from the tracker.madx template."""
    with open(filename, 'w') as out:
        for line in fileinput.input(c.madxdir+'tracker.madx', inplace=0):
            if searchExp in line:
                line = line.replace(searchExp,replaceExp)
            out.write(line)

def writeTracker(k,data):
    """Creates the text to replace ADAPTTHISPART in tracker.madx.

    To simulate a linear sweep plus a perturbation in the form of a
    damped sinusoid.

    Some comments:
    - Changing Nturns without changing start and end tune actually
      changes the tune sweep speed, is that intended?
    - Commented out lines of code can be removed, right?
    - Most replacement seems rather unnecessary with the exception
      of specifying tracksdir and initial coordinates. The other
      constants can just be written to a file that's read by MADX. If we
      would split the input coordinates into different files there would
      be no more need to generate the full script for every batch, we
      could just generate (shorter) input files. This would also make
      the MAD-X templates easier to maintain.
    """
    line= "m_f = (kqf1_start - kqf1_end) / (1 - "+str(c.Nturns)+");\n"
    line+="m_d = (kqd_start -  kqd_end ) / (1 - "+str(c.Nturns)+");\n\n"

    line+="n_f = kqf1_start - m_f;\n"
    line+="n_d = kqd_start - m_d;\n\n"

    if(c.dataripple):
        line+= "READTABLE, FILE='"+c.inputdir+c.ripplefile+".prt';\n"

    line+="tr$macro(turn): macro = {"
    line+='\n'
    if(c.ripple):
        line+='\n'
        line+="IF(turn>"+str(c.startturn)+")"
        line+='\n'
        #line+="{ pert=3e-7*SIN(2*PI*turn/500)*EXP((5000-turn)/1000);}"
        line+="{\n"
        if(c.expperiod==0):
            line+="pert=-"+str(c.amplitude)+"e-8*SIN(2*PI*turn/"+str(c.period)+");\n"
        else:
            line+="pert=-"+str(c.amplitude)+"e-8*SIN(2*PI*turn/"+str(c.period)+")*EXP(("+str(c.startturn)+"-turn)/"+str(c.expperiod)+");\n"
        line+="};\n"
        #line+='\n'
        line+="ELSE{pert=0;};"
        line+='\n'
        line+="kqf1= m_f * turn + n_f + pert;\n"
    elif(c.dataripple):
        line+='\n'
        line+="SETVARS, TABLE=RIPPLE, ROW=turn;\n\n"
        line+="kqf1= m_f * turn + n_f + RIPPLEVAL;\n"
        line+='\n'
    else:
        line+='\n'
        line+="kqf1= m_f * turn + n_f;"
        line+='\n'
    line+="kqd = m_d * turn + n_d;"
    line+='\n'
    line+="};"
    line+='\n\n'
    line+="track,onepass, aperture,  recloss, file='"+c.tracksdir+str(k)+"/track.batch"+str(k)+"', update;"
    line+='\n\n'
    for i in range(k*c.Nparperbatch,(k+1)*c.Nparperbatch):
        line+="start, x="+str(data[0][i])+", px= "+str(data[1][i])+", y= "+str(data[2][i])+", py= "+str(data[3][i])+", t= "+str(0)+", pt= "+str(data[4][i])+";"
        line+="\n"
    for place in c.elements:
        line+= "observe, place="+place+';';
        line+='\n'
    #line+= "run,turns="+str(Nturns)+", ffile =" +str(int(Nturns/100))+";"
    line+= "run,turns="+str(c.Nturns)+",ffile="+str(c.turnmultiplicity)+";" #every turn logging
    line+='\n\n'
    line+='endtrack;'
    line+='\n\n'
    if(c.writetrack):
        line+="write, table = trackloss, file='"+c.lossdir+"batch"+str(k)+"';"

    return line

def writeTabler():
    """Seems unused. To be deleted?"""
    return 0


def onerun():
    """Submits jobs for simulation defined in main.

    Needs some fixing:
    - The replacers inside the createTwiss part will have worked only
      the first time he ran it, since replace() is defined as in-place
      replacement. After that they were useless.
      Moreover, it looks like the difference between initial and final
      tune did not take number of turns into account. Probably best to
      set initial and final tune as standard and then change
      "run,turns=6;" -> "run,turns=Nturns;", and then save a replaced
      file somewhere else.
      Also, it seems like a bit of a waste of time to run a tracking to
      create the tune_var table. If you agree the commented tables
      aren't useful then we can find a quicker way. (A simple loop
      should do the trick.)
    - Replacer in the tracking has the same problem, but with the
      relative path 'homedir' replacement has become redundant anyhow.
      We can delete that line.
    - The MAD-X code used should somehow be linked to the
      slowExtractionMADX code. (Can be done with an out-of-place
      replacer and the c.madxdir variable.)
    - We may want to pick non-random seeds for comparing loss-reduction
      methods.

    And then the adapter needs to be made much more flexible, or
    different template options added, so that we can run new simulations
    instead of just tune ripple.
    """
    #Redirect the output to DEVNULL
    oldstdout = sys.stdout
    FNULL = open(os.devnull, 'w')
    LOG = open("log.txt","w")

    #Generate twiss files before and after thinning, used to make initial distributions
    if(c.createTwiss):
        print c.madxdir, c.home
        replacer(c.madxdir+'TableTemplate.madx','homedir',c.home)
        replacer(c.madxdir+'TableTemplate.madx','Nturns',str(c.Nturns))
        print 'running: ', "madxColl<"+c.madxdir+"TableTemplate.madx"
        subprocess.Popen("madxColl<"+c.madxdir+"TableTemplate.madx", stdout=LOG, shell=True).wait()
        print 'Table creation is finished!'

    if(c.trackingBool):
        if not os.path.exists(c.data):
                os.makedirs(c.data)
        #Generate initial particle distributions
        data=dis.get_gauss_distribution(output=c.data+'initial_distribution_gauss', input=c.twissdir+'twiss_after_thinning.prt',n_part=c.Nbatches*c.Nparperbatch, sigmas=6, beam_t='FT', seed=np.random.randint(9999))
        print 'Initial particles distributions created!'

        if not os.path.exists(c.tracksdir):
                os.makedirs(c.tracksdir)
        if not os.path.exists(c.lossdir):
                os.makedirs(c.lossdir)
        if not os.path.exists(c.jobsdir):
                os.makedirs(c.jobsdir)

        for i in range(c.Nbatches):
            if not os.path.exists(c.tracksdir+str(i)):
                os.makedirs(c.tracksdir+str(i))

        #Prepare madx files for batching
        replacer(c.madxdir+'tracker.madx','homedir',c.home)
        for nr in range(c.Nbatches):
            file=c.jobsdir+'batch'+str(nr)+'.madx'
            searchExp='ADAPTTHISPART'
            replaceExp=writeTracker(nr,data)
            adapter(file,searchExp,replaceExp)
        print 'Input madx files created!'

        #Actual particle tracking
        for nrb in range(c.Nbatches):
            print 'Job '+str(nrb)+' of '+str(c.Nbatches)+', containing '+ str(c.Nparperbatch)+' particles, started.'

            if(c.LXplus):
                subprocess.Popen("bsub -q "+whichQueue()+" madx "+c.jobsdir+"batch"+str(nrb)+".madx", stdout=LOG, shell=True).wait()
                print "Job submitted!"
            else:
                print 'running: ', "madxColl<"+c.jobsdir+"batch"+str(nrb)+".madx"
                subprocess.Popen("madxColl<"+c.jobsdir+"batch"+str(nrb)+".madx", stdout=LOG, shell=True).wait()
    print "Done!"

def main():
    """Set parameters for the simulation and let it be batched.

    At the moment you have to change and comment settings for different
    studies. At some point we can change this to take command line input.
    """

#    print '##########################################################'
#    print '#####   Glitches with 50k particles and 34095 turns  #####'
#    print '##########################################################'
#    #Glitches
#    c.SetGeneral(f_Nturns=34095,f_Nbatches=1000,f_Nparperbatch=50,f_turnmultiplicity=50)
#    c.SetBools(f_ripple=True,f_dataripple=False, f_writetrack=True)
#
#    c.SetRipple(f_startturn=5000,f_amplitude=1000,f_period=500,f_expperiod=500)
#    c.SetDirs(f_name="misalligned")
#    onerun()
#
#    c.SetRipple(f_startturn=5000,f_amplitude=2000,f_period=500,f_expperiod=500)
#    c.SetDirs(f_name=None)
#    onerun()
#
#    c.SetRipple(f_startturn=5000,f_amplitude=1000,f_period=500,f_expperiod=1000)
#    c.SetDirs(f_name=None)
#    onerun()
#
#    print '##########################################################'
#    print '## Ripple from data with 50k particles and 34095 turns  ##'
#    print '#############  Without harmonics suppression  ############'
#    print '##  or (50,100,150)Hz or (50,100,150,300)Hz supression  ##'
#    print '##########################################################'
#
#    c.SetGeneral(f_Nturns=34095,f_Nbatches=1000,f_Nparperbatch=50,f_turnmultiplicity=10000)
#    c.SetBools(f_ripple=False,f_dataripple=True, f_writetrack=True)
#    c.setripplefile('ripple0')
#    c.SetDirs(f_name='realripple_50k_pure')
#    onerun()
#
#    c.setripplefile('ripple50100150')
#    c.SetDirs(f_name='realripple_50k_comp50100150')
#    onerun()
#
#    c.setripplefile('ripple300')
#    c.SetDirs(f_name='realripple_50k_comp50100150300')
#    onerun()

#    print '##########################################################'
#    print '#######  Emittance scan with 50 (freq,amp)-points  #######'
#    print '############  50k particles and 34095 turns  #############'
#    print '##########################################################'
#    #Emittance scan
#    freq=np.random.uniform(10,700, size=100)
#    freq=[10,10,10,350,350,350,700,700,700]
#    freq=np.asarray(freq)
#    amp=np.random.uniform(0.67,100, size=100)*1e-6
#    turnsperiod=np.rint(freq**(-1.)/0.0000233)
#    quadamp=np.rint(amp*1777.*840) #units e-8
#    quadamp=[1,50,150,1,50,150,1,50,150]
#    quadamp=np.asarray(quadamp)
#    for i,j in zip(turnsperiod,quadamp):
#        c.SetGeneral(f_Nturns=34095,f_Nbatches=200,f_Nparperbatch=50,f_turnmultiplicity=10000)
#        c.SetBools(f_ripple=True,f_dataripple=False, f_writetrack=True)
#        c.SetRipple(f_startturn=1000,f_amplitude=j,f_period=i,f_expperiod=0)
#        c.SetDirs(f_name=None)
#        onerun()
#
    print '##########################################################'
    print '###  Linear Sweep with 50k particles and 34095*6 turns ###'
    print '##########################################################'
    #Long Lin Sweep
    c.SetGeneral(f_Nturns=34095*6,f_Nbatches=2000,f_Nparperbatch=50,f_turnmultiplicity=10000)
    c.SetBools(f_ripple=False,f_dataripple=False, f_writetrack=True)
    c.SetDirs(f_name='runwithoutvirtual200k')
    onerun()



def tester():
    """Test of the main program functionality.

    Runs a simple test to see if everything is working as it should. We
    should think about replacing this with our (still to be defined)
    test-case(s).
    """
    print 'test'

    c.createTwiss=True
    c.trackingBool=True
    c.LXplus=False

    c.SetGeneral(f_Nturns=3,f_Nbatches=1,f_Nparperbatch=50,f_turnmultiplicity=10000)
    c.SetBools(f_ripple=False,f_dataripple=False, f_writetrack=True)
    # c.setripplefile('ripple0')
    c.SetDirs(f_name='ATEST_fra')
    onerun()

tester()
# main()
