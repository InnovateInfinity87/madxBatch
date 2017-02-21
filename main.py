# -*- coding: utf-8 -*-
"""
Created on Thu Jul 14 17:31:19 2016

@author: wvandepo
"""
import fileinput
import subprocess
import python.make_ps_distribution as dis
import os
import sys
import numpy as np
import python.Constants as c

def whichQueue(turns):
    queue = "1nh"
    if turns > 10e6:
        queue = "8nh"
    if turns > 10e6:
        queue = "1nd"
    elif turns > 2e7:
        queue = "2nd"
    elif turns > 7e7:
        queue = "1nw"
    elif turns > 1e8:
        queue = "2nw"
    return queue
 

def replacer(filename,searchExp,replaceExp):
    for line in fileinput.input(filename, inplace=1):
        if searchExp in line:
            line = line.replace(searchExp,replaceExp)
        sys.stdout.write(line)
        
def adapter(filename,searchExp,replaceExp):
    out=open(filename, 'w')
    #print filename
    for line in fileinput.input(c.madxdir+'tracker.madx', inplace=0):
        if searchExp in line:
            line = line.replace(searchExp,replaceExp)
        out.write(line)
        
def writeTracker(k,data):
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
        line+="}\n"
        #line+='\n'
        line+="ELSE{pert=0;}"
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
    return 0
    
    
def onerun():
    #Redirect the output to DEVNULL
    oldstdout = sys.stdout
    FNULL = open(os.devnull, 'w')
    LOG = open("log.txt","w")
    
    #Generate twiss files before and after thinning, used to make initial distributions
    if(c.createTwiss):
        replacer(c.madxdir+'TableTemplate.madx','homedir',c.home)
        replacer(c.madxdir+'TableTemplate.madx','Nturns',str(c.Nturns)) 
        subprocess.Popen("madx<"+c.madxdir+"TableTemplate.madx", stdout=LOG, shell=True).wait()
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
                subprocess.Popen("bsub -q "+whichQueue(c.Nturns*c.Nparperbatch)+" madx "+c.jobsdir+"batch"+str(nrb)+".madx", stdout=LOG, shell=True).wait()
                print "Job submitted!"
            else:
                subprocess.Popen("madx<"+c.jobsdir+"batch"+str(nrb)+".madx", stdout=LOG, shell=True).wait()
    print "Done!"
    
def main():
    
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
    print 'test'
    c.SetGeneral(f_Nturns=34095,f_Nbatches=1,f_Nparperbatch=50,f_turnmultiplicity=10000)
    c.SetBools(f_ripple=False,f_dataripple=True, f_writetrack=True)
    c.setripplefile('ripple0')
    c.SetDirs(f_name='ATEST_realripple_50k_pure')
    onerun()
    
#tester()   
main()
