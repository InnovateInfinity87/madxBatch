# -*- coding: utf-8 -*-
"""
Batched MAD-X slow extraction

Modifies files for running MAD-X simulations of slow extraction with a
tune ripple and then submits them to the CERN HTCondor batching service.

@author: Linda Stoel, Wouter van de Pontseele
"""
import os
import stat
import sys
import fileinput
import subprocess
from shutil import copyfile

import make_ps_distribution as dis
from localsub import localsub

class Settings:
    """Settings for the batched simulations"""
    def __init__(self, name, outputdir=None, studygroup='', disk='eos'):
        if not disk in ['eos', 'afsprivate', 'afspublic']:
            print("Disk setting not valid. valid options are 'eos'/'afsprivate'"+
                  "/'afspublic'. Defaulting to eos, unless desired outputdir was given.")
            disk='eos'
        if outputdir==None:
            user = os.environ["USER"]
            if disk=='eos':
                outputdir = '/eos/user/'+user[0]+'/'+user+'/madxBatchData/'+studygroup
            elif disk=='afsprivate':
                outputdir = '/afs/cern.ch/work/'+user[0]+'/'+user+'/private/madxBatchData/'+studygroup
            else:
                outputdir = '/afs/cern.ch/work/'+user[0]+'/'+user+'/public/madxBatchData/'+studygroup

        self.name = name
        self.datadir = outputdir+"/"+name+"/"
        self.home = '/'.join(os.path.dirname(os.path.abspath(__file__)).split('/')[:-1])
        self.pycolldir = self.home+'/../pycollimate/'
        self.slowexdir = self.home+'/../slowExtractionMADX'
        self.slowexfiles = (self.home+'/input/sps.ele, '+
                            self.home+'/input/aperturedb_1.dbx, '+
                            self.home+'/input/aperturedb_3.dbx, '+
                            self.slowexdir+'/cmd/matchtune.cmdx, '+
                            self.slowexdir+'/cmd/matchchroma.cmdx')

        self.trackertemplate = self.home+"/madx/tracker_nominal_template.madx"
        self.myreplace = {}
        self.twissfile = None
        self.local = False

        self.seed = None

        self.nturns=10000
        self.nbatches=1000
        self.nparperbatch=100
        self.ffile=100 #ffile MADX
        self.flavour=None

        self.elements=[]

        self.trackingbool = True
        self.dynamicbump=False
        self.pycollimate = True
        self.septadb = self.home+"/input/septa_DB_nominal.tfs"
        self.septadbreplace = None

        self.saveloss = True
        self.saveout = False

        self.monitor = False
        
        self._custom_datadir = False
        
    def set_name(self, name):
        self.name = name
        if not self._custom_datadir:
            outputdir = '/'.join(self.datadir.split('/')[:-2])
            self.datadir = outputdir+'/'+name+"/"
        
    def set_studygroup(self, studygroup):
        if not self._custom_datadir:
            temp = self.datadir.split('/')
            temp[-3] = studygroup
            self.datadir = '/'.join(temp)+'/'

    def set_datadir(self, dir):
        self.datadir = dir
        self._custom_datadir = True
        
    def set_params(nturns=None, nbatches=None, nparperbatch=None, ffile=None):
        if nturns is not None:
            self.nturns = nturns
        if nbatches is not None:
            self.nbatches = nbatches
        if nparperbatch is not None:
            self.nparperbatch = nparperbatch
        if ffile is not None:
            self.ffile = ffile


def flavour(nturns, parperbatch,pycollimate):
    """Determines appropriate job flavour based on particles*turns."""
    #TODO: benchmark for better estimate
    n = nturns*parperbatch
    if pycollimate:
        if n > 144e4:
            f = "nextweek"
        elif n > 48e4:
            f = "testmatch"
        elif n > 16e4:
            f = "tomorrow"
        elif n > 4e4:
            f = "workday"
        elif n > 2e4:
            f = "longlunch"
        elif n > 0.66e4:
            f = "microcentury"
        else:
            f = "espresso"
    else:
        if n > 48e6:
            f = "nextweek"
        elif n > 16e6:
            f = "testmatch"
        elif n > 5.33e6:
            f = "tomorrow"
        elif n > 1.33e6:
            f = "workday"
        elif n > 0.66e6:
            f = "longlunch"
        elif n > 0.22e6:
            f = "microcentury"
        else:
            f = "espresso"

    return f


def replacer(filename,searchExp,replaceExp):
    """In-place replacement in text files.

    Finds all occurences of searchExp in the file and replaces each one
    with replaceExp, overwriting the original file.
    """
    for line in fileinput.input(filename, inplace=1):
        line = line.replace(searchExp,replaceExp)
        sys.stdout.write(line)

def customize_tracker(filename,searchExp,replaceExp):
    """Creates a new file adapted from the tracker.madx template."""
    with open(filename, 'w') as out:
        for line in fileinput.input('tracker.madx', inplace=0):
            line = line.replace(searchExp,replaceExp)
            out.write(line)
        out.flush()
        os.fsync(out.fileno())


def track_lin(k,data,settings):
    """Creates the text to replace pyTRACKER in tracker.madx. (nominal case)"""
    line = ("m_f = (kqf1_start - kqf1_end) / (1 - "+str(settings.nturns)+");\n"+
            "m_d = (kqd_start -  kqd_end ) / (1 - "+str(settings.nturns)+");\n\n"+

            "n_f = kqf1_start - m_f;\n"+
            "n_d = kqd_start - m_d;\n\n"+
        
            "SYSTEM, 'mkdir "+str(k)+"';\n\n")

    if settings.dynamicbump:
        line += (orthogonal_bumps()+"\n\n"+

                "x_knob_start = 406*dpp_start + (-0.0765*dpp_start*1e3);\n"+
                "px_knob_start = (-20.667*dpp_start) + (-0.0022*start*1e3);\n"+
                "EXEC, obump(x_knob_start, px_knob_start);\n"+
                "dpp_inc = (dpp_end-dpp_start)/"+str(settings.nturns)+";\n"
                "x_knob_inc = 406*dpp_inc + (-0.0765*dpp_inc*1e3);\n"+
                "px_knob_inc = (-20.667*dpp_inc) + (-0.0022*dpp_inc*1e3);\n\n")
        
    line += ("tr$macro(turn): MACRO = {\n"+
             " kqf1 = m_f * turn + n_f;\n"+
             " kqd = m_d * turn + n_d;\n")
    if settings.dynamicbump:
        line += (" EXEC, obump(x_knob_inc, px_knob_inc);\n")
    line += "};\n\n"

    line += ("OPTION, -WARN;\n"+
             "TRACK, ONEPASS, APERTURE, RECLOSS, "+
                  "FILE='"+str(k)+"/track.batch"+str(k)+"', UPDATE;\n\n")

    for i in range(k*settings.nparperbatch,(k+1)*settings.nparperbatch):
        line += (" START, X ="+str(data[0][i])+", "+
                       "PX = "+str(data[1][i])+", "+
                       "Y = "+str(data[2][i])+", "+
                       "PY = "+str(data[3][i])+", "+
                       "T = "+str(0)+", "+
                       "PT = "+str(data[4][i])+";\n")
    for place in settings.elements:
        line += " OBSERVE, PLACE = "+place+';\n';
    line += (" RUN, TURNS="+str(settings.nturns)+", "+
                 "MAXAPER={0.1,0.01,0.1,0.01,"+str(0.03*settings.nturns)+",0.1}, "+
                 "FFILE="+str(settings.ffile)+";\n\n"+

             "ENDTRACK;\n"+
             "OPTION, WARN;\n\n")

    if(settings.saveloss):
        line+="WRITE, TABLE = trackloss, FILE = 'losses.tfs';\n\n"

    line += "SYSTEM, 'tar -czf tracks.tar.gz "+str(k)+"';"

    return line


def orthogonal_bumps():
    "MAD-X code for the orthogonal bump macro"
    return ("EOPTION, ADD = TRUE;\n\n"+

            "obump(x_knob_value, px_knob_value) : MACRO = {\n"+
            " SELECT, FLAG=ERROR, CLEAR;\n"+
            ' SELECT, FLAG=ERROR, RANGE=MPSH_rb.21202..1/MPSH_rb.21202..4;\n'+
            " EFCOMP, ORDER:=0, DKN= 2.331301861e-05 * x_knob_value * (-1/4) + 9.929676756e-06 * px_knob_value * (-1e2/4);\n"+
            " SELECT, FLAG=ERROR, CLEAR;\n"+
            ' SELECT, FLAG=ERROR, RANGE=MPLH_rb.21431..1/MPLH_rb.21431..4;\n'+
            " EFCOMP, ORDER:=0, DKN= -9.91281169e-06   * x_knob_value * (-1/4) + 2.600240449e-07   * px_knob_value * (-1e2/4);\n"+
            " SELECT, FLAG=ERROR, CLEAR;\n"+
            ' SELECT, FLAG=ERROR, RANGE=MPNH_rb.21732..1/MPNH_rb.21732..4;\n'+
            " EFCOMP, ORDER:=0, DKN=  2.98e-05 * x_knob_value * (-1/4) + 2.174863504e-05 * px_knob_value * (-1e2/4);\n"+
            " SELECT, FLAG=ERROR, CLEAR;\n"+
            ' SELECT, FLAG=ERROR, RANGE=MPLH_rb.21995..1/MPLH_rb.21995..4;\n'+
            " EFCOMP, ORDER:=0, DKN= -1.923181311e-05 * x_knob_value * (-1/4) + (-5.128495572e-06) * px_knob_value * (-1e2/4);\n"+
            " SELECT, FLAG=ERROR, CLEAR;\n"+
            ' SELECT, FLAG=ERROR, RANGE=MPLH_rb.22195..1/MPLH_rb.22195..4;\n'+
            " EFCOMP, ORDER:=0, DKN= 2.105219176e-05 * x_knob_value * (-1/4) + 8.373821994e-06 * px_knob_value * (-1e2/4);\n"+
            "};")


def submit_job(settings):
    """Creates and submits job cluster for simulation defined by settings.

    Needs some fixing:
    - Ideally, the MAD-X template used should somehow be linked to the
      slowExtractionMADX code.
    - We may want to pick non-random seeds for comparing loss-reduction
      methods.
    - The adapter needs to be made much more flexible, and/or different
      template options added.
    - The generated twiss table should be generated in output, not input.
    """
    if os.path.exists(settings.datadir):
        settings.datadir = settings.datadir[:-1]+'_/'
        i = 1
        while True:
            i += 1
            if not os.path.exists(settings.datadir[:-1]+str(i)+"/"):
                break
        settings.datadir = settings.datadir[:-1]+str(i)+"/"
    os.makedirs(settings.datadir)
    os.chdir(settings.datadir)

    #Generate twiss files before and after thinning, used to make initial distributions
    copyfile(settings.trackertemplate, "tracker.madx")
    replacer("tracker.madx",'pyDATADIR', settings.datadir)
    replacer("tracker.madx", 'pyHOMEDIR', settings.home)
    replacer("tracker.madx", 'pySLOWEXDIR', settings.slowexdir)
    replacer("tracker.madx", 'pyPYCOLL', str(int(settings.pycollimate)))
    for key, replacement in settings.myreplace.iteritems():
        replacer("tracker.madx", key, replacement)
    if settings.twissfile is None:
        log = open("twisslog.txt", 'w')
        copyfile("tracker.madx", "table.madx")
        replacer("table.madx", '/*pyTWISS', '')
        replacer("table.madx", 'pyTWISS*/', '')
        print "Creating Twiss table"
        if settings.pycollimate:
            subprocess.Popen(settings.pycolldir+"madxColl<table.madx", stdout=log, shell=True).wait()
        else:
            subprocess.Popen("/afs/cern.ch/user/m/mad/bin/madx_dev<table.madx", stdout=log, shell=True).wait()
        log.close()
        os.remove("table.madx")
        print 'Table creation is finished!'
    else:
        copyfile(settings.twissfile, settings.datadir+"/thin_twiss.tfs")

    if(settings.trackingbool):
        #Generate initial particle distribution
        #TODO: Non-random input
        data=dis.get_gauss_distribution(output=settings.datadir+'initial_distribution',
                                        input=settings.datadir+"/thin_twiss.tfs",
                                        n_part=settings.nbatches*settings.nparperbatch,
                                        sigmas=6, beam_t='FT',
                                        seed=settings.seed,
                                        file_head=settings.home+"/input/distributionheader.txt")
        print 'Initial particle distribution created!'

        #Create datastructure
        os.mkdir("tracks")
        os.mkdir("losses")
        os.mkdir("jobs")
        os.mkdir("output")
        os.mkdir("error")

        copyfile(settings.home+"/other/unpacker.sh", "tracks/unpacker.sh")
        st = os.stat("tracks/unpacker.sh")
        os.chmod("tracks/unpacker.sh", st.st_mode | stat.S_IEXEC)

        #Customize septa db if needed
        if settings.septadbreplace is not None:
            copyfile(settings.septadb, "septa_DB_custom.tfs")
            for key, replacement in settings.septadbreplace.iteritems():
                replacer("septa_DB_custom.tfs", key, replacement)

        #Prepare job files for each batch
        copyfile(settings.home+"/other/job.sh", "jobs/start.sh")
        if settings.pycollimate:
            replacer("jobs/start.sh", 'MADXEXE', "madxColl")
        else:
            replacer("jobs/start.sh", 'MADXEXE', "madx_dev")
        for i in range(settings.nbatches):
            outfile="jobs/"+str(i)+".madx"
            searchExp='pyTRACKER'
            replaceExp=track_lin(i,data,settings)
            customize_tracker(outfile,searchExp,replaceExp)
        os.remove("tracker.madx")
        print 'Job files created!'

        #Write submit file
        with open(settings.name+".sub",'w') as subfile:
            subfile.write("executable = jobs/start.sh\n")
            if settings.pycollimate:
                subfile.write('requirements = (OpSysAndVer =?= "CentOS7")\n')
                subfile.write("transfer_input_files = "+
                              "jobs/$(ProcId).madx, "+
                              settings.pycolldir+"madxColl, "+
                              settings.pycolldir+"track_inside_coll.py, "+
                              settings.pycolldir+"pycollimate.py, "+
                              (settings.septadb if (settings.septadbreplace is None) else "septa_DB_custom.tfs")+", "+
                              settings.home+"/other/matplotlibrc, "+
                              settings.slowexfiles+"\n")
            else:
                subfile.write("transfer_input_files = "+
                              "jobs/$(ProcId).madx, "+
                              "/afs/cern.ch/user/m/mad/bin/madx_dev, "+
                              settings.slowexfiles+"\n")
            subfile.write('arguments = "$(ProcId)"\n')
            subfile.write("initialdir = "+settings.datadir+"\n")
            if settings.saveout:
                subfile.write("output = output/$(ProcId).out\n")
            subfile.write("error = error/$(ProcId).err\n")
            subfile.write("log = log.txt\n")
            subfile.write('transfer_output_remaps = "'+
                          'tracks.tar.gz=tracks/$(ProcId).tar.gz; '+
                          'losses.tfs=losses/$(ProcId).tfs"\n')
            if settings.flavour is None:
                subfile.write("+JobFlavour = '"+flavour(settings.nturns, settings.nparperbatch, settings.pycollimate)+"'\n")
            else:
                subfile.write("+JobFlavour = '"+settings.flavour+"'\n")
            subfile.write("\nqueue "+str(settings.nbatches))
        print 'Submit file created!'

        #Submit job
        if settings.local:
            localsub(settings.name+".sub")
            print 'Local batching completed!'
        else:
            subprocess.check_call("condor_submit "+settings.name+".sub", shell=True)
            if settings.monitor:
                subprocess.check_call("condor_wait -status "+settings.datadir+"log.txt", shell=True)


def tester():
    """Test of the main program functionality.

    Runs a simple test to see if everything is working as it should. We
    should think about replacing this with our (still to be defined)
    test-case(s).
    """
    print "Running tester()"

    settings=Settings('6900', studygroup='difftest', disk='afspublic')

    settings.trackingbool=True
    settings.trackertemplate=settings.home+"/madx/tracker_diffuser_template.madx"
    settings.local=False
    settings.monitor=False
    settings.seed = 0

    settings.elements=["QFA_START",'AP.UP.ZS21633']#,'AP.DO.ZS21676','AP.UP.TPST21760']
    #settings.elements=['AP.UP.ZS21633_M']#,'AP.DO.ZS21676_M','AP.UP.TPST21760']#,'TCE.21695']

    settings.nturns=10000
    settings.nbatches=1000
    settings.nparperbatch=20
    settings.ffile=1

    settings.dynamicbump=True
    settings.pycollimate=True

    settings.flavour=None

    submit_job(settings)

    return

if __name__ == "__main__":
    tester()
