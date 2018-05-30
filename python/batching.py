# -*- coding: utf-8 -*-
"""
Batched MAD-X slow extraction

Modifies files for running MAD-X simulations of slow extraction with a
tune ripple and then submits them to the CERN HTCondor batching service.

@author: Linda Stoel, Wouter van de Pontseele
"""
from version import __version__

import os
import stat
import sys
import fileinput
import subprocess
from shutil import copyfile
import scipy.stats as stats

import make_ps_distribution as dis
from localsub import localsub

class Settings:
    """Settings for the batched simulations"""
    def __init__(self, name, outputdir=None, studygroup='', disk=None):
        if outputdir==None:
            if not disk in ['afsprivate', 'afspublic', 'afsproject']:
                print("Disk setting not valid. Valid options are afsprivate'"+
                      "/'afspublic'/'afsproject'. Defaulting to afsprivate."+
                      "(EOS is not yet supported.)")
            user = os.environ["USER"]
            if disk=='afsproject':
                outputdir = '/afs/cern.ch/project/sloex/'+studygroup
                if not os.path.exists(outputdir):
                    print("ERROR: To submit via afsproject a valid studygroup must be given. "+
                          "Take a look at the folder structure in /afs/cern.ch/project/sloex "+
                          "to find an appropriate one or contact the admins to make one. Exiting.")
                    exit()
            elif disk=='afspublic':
                outputdir = '/afs/cern.ch/work/'+user[0]+'/'+user+'/public/madxBatchData/'+studygroup
            else:
                outputdir = '/afs/cern.ch/work/'+user[0]+'/'+user+'/private/madxBatchData/'+studygroup
        else:
            if disk is not None:
                print "Disk setting is not used because outputdir was specified."

        self.name = name
        self.datadir = outputdir+"/"+name+"/"
        self.home = '/'.join(os.path.dirname(os.path.abspath(__file__)).split('/')[:-1])
        self.pycolldir = self.home+'/../pycollimate/'
        self.madxversion = '/afs/cern.ch/user/m/mad/bin/rel/5.03.06/madx-linux64-gnu'

        self.trackertemplate = self.home+"/madx/tracker_nominal_template.madx"
        self.trackerrep = track_lin
        self.thickchanges = None
        self.thinchanges = None
        self.finalchanges = None
        self.myreplace = {}
        self.twissfile = None
        self.local = False

        self.seed = None

        self.nturns=10000
        self.nbatches=1000
        self.nparperbatch=100
        self.ffile=100 #ffile MADX
        self.dppmax=None
        self.slices=None
        self.slicewidth=2.5E-4

        self.beam_t='FT'
        self.beamkwargs={'n_sigma': 6}

        self.flavour=None

        self.elements=[]

        self.trackingbool = True
        self.pycollimate = True
        self.septadb = self.home+"/input/septa_DB_nominal.tfs"
        self.septadbreplace = None
        self.pcblack = False

        # knob_x_bump = offx + cx*dpp (millimeter)
        # knob_px_bump = offpx*dpp  (microrad)
        self.dynamicbump=False
        self.dynamicbump_offx = 0.0
        self.dynamicbump_cx = -729.5
        self.dynamicbump_offpx = 0.0
        self.dynamicbump_cpx = 42600.0

        self.cose = False

        self.saveloss = True
        self.saveout = False
        self.savetracks = True
        self.trackendonly = None

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
        
    def set_params(self, nturns=None, nbatches=None, nparperbatch=None, ffile=None):
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


def tune_setup(settings):
    if settings.cose:
        line = ('qh_set_end = qh_setvalue;\n'+
                'kqf1_end = kqf1 * (1.0 + dpp_end/(1+dpp_end));\n'+
                'kqd_end = kqd * (1.0 + dpp_end/(1+dpp_end));\n\n'+

                'qh_set_start = qh_setvalue;\n'+
                'kqf1_start = kqf1 * (1.0 + dpp_start/(1+dpp_start));\n'+
                'kqd_start = kqd * (1.0 + dpp_start/(1+dpp_start));\n')
    else:
        line = ('qh = qh_end;\n'+
                "CALL, FILE='pyHOMEDIR/madx/op_matchtune_h.cmdx';\n"+
                'qh_set_end = qh_setvalue;\n'+
                'kqf1_end = kqf1;\n'+
                'kqd_end = kqd;\n\n'+

                'qh = qh_start;\n'+
                "CALL, FILE='pyHOMEDIR/madx/op_matchtune_h.cmdx';\n"+
                'qh_set_start = qh_setvalue;\n'+
                'kqf1_start = kqf1;\n'+
                'kqd_start = kqd;\n')
    return line


def track_lin(k,data,settings):
    """Creates the text to replace pyTRACKER in tracker.madx. (nominal case)"""
    line = ('c_f = (kqf1_end-kqf1_start)/('+str(settings.nturns)+'-1);\n'+
            'c_d = (kqd_end-kqd_start)/('+str(settings.nturns)+'-1);\n'+
            'c_dpp = (dpp_end-dpp_start)/('+str(settings.nturns)+'-1);\n\n'+
        
            "SYSTEM, 'mkdir "+str(k)+"';\n\n")
        
    line += ('tr$macro(turn): MACRO = {\n'+
             ' kqf1 = kqf1_start + c_f*(turn-1);\n'+
             ' kqd = kqd_start + c_d*(turn-1);\n'+
             ' dpp_turn = dpp_start + c_dpp*(turn-1);\n')

    if settings.cose:
        line += (' relerr = dpp_turn/(1+dpp_turn);\n'+ # Let's not fuss about dpp vs PT for now...
                 ' abserr = relerr*kMBA/4;\n' #kMBA=kMBB
                 ' SELECT, FLAG=ERROR, CLEAR;\n'+
                 ' SELECT, FLAG=ERROR, PATTERN="MB.*";\n'+
                 ' EFCOMP, ORDER=0, DKN={abserr};\n')
             
    if settings.dynamicbump:
        line += (" knob_x_bump = "+str(settings.dynamicbump_offx)+"+("+str(settings.dynamicbump_cx)+")*dpp_turn;\n"+
                 " knob_px_bump = "+str(settings.dynamicbump_offpx)+"+("+str(settings.dynamicbump_cpx)+")*dpp_turn;\n"+
                 " EXEC, lss2bump(knob_extr_bump, knob_x_bump, knob_px_bump);\n")
    line += "};\n\n"

    line += "EXEC, tr$macro(0);\n\n"

    line += ("OPTION, -WARN;\n"+
             "TRACK, ONEPASS, APERTURE, UPDATE, RECLOSS")
    if settings.trackendonly is None and settings.savetracks:
        line += (", FILE='"+str(k)+"/track.batch"+str(k)+"';\n\n")
    else:
        line += (";\n\n")

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
                 "MAXAPER={0.5,0.05,0.5,0.05,"+str(0.03*settings.nturns)+",0.1}, "+
                 "FFILE="+str(settings.ffile)+";\n\n"+

             "ENDTRACK;\n"+
             "OPTION, WARN;\n\n")

    if settings.trackendonly is not None and settings.savetracks:
        line+="CALL, FILE='"+settings.home+"/madx/end_only.cmdx';\n"

        for obsnum in range(len(settings.elements)+1):
            for partnum in range(settings.nparperbatch):
                temp = "obs"+str(obsnum+1).zfill(4)+".p"+str(partnum+1).zfill(4)
                line+="EXEC, end_only("+str(settings.trackendonly)+", track."+temp+", '"+str(k)+"/track.batch"+str(k)+"."+temp+"');\n"
        line+="\n"

    if(settings.saveloss):
        line+="WRITE, TABLE = trackloss, FILE = 'losses.tfs';\n\n"

    line += "SYSTEM, 'tar -czf tracks.tar.gz "+str(k)+"';"

    return line


def track_sliced(k,data,settings):
    """Creates the text to replace pyTRACKER in tracker.madx. (sliced nominal case)"""
    dpp = str(settings.slices[k/settings.nbatches])

    line = "SYSTEM, 'mkdir "+str(k)+"';\n\n"

    if settings.cose:
        line += ('relerr = ('+dpp+')/(1+('+dpp+'));\n'+
                 'abserr = relerr*kMBA/4;\n' #kMBA=kMBB
                 'SELECT, FLAG=ERROR, CLEAR;\n'+
                 'SELECT, FLAG=ERROR, PATTERN="MB.*";\n'+
                 'EFCOMP, ORDER=0, DKN={abserr};\n')

    if settings.dynamicbump:
        line += (" knob_x_bump = "+str(settings.dynamicbump_offx)+"+("+str(settings.dynamicbump_cx)+")*("+dpp+");\n"+
                 " knob_px_bump = "+str(settings.dynamicbump_offpx)+"+("+str(settings.dynamicbump_cpx)+")*("+dpp+");\n"+
                 " EXEC, lss2bump(knob_extr_bump, knob_x_bump, knob_px_bump);\n")

    line += ("dpp_matchtune = "+dpp+";\n"+
             "qh = qh_res;\n"+
             "CALL, FILE='"+settings.home+"/madx/op_matchtune_h_offmom.cmdx';\n\n")

    line += ("OPTION, -WARN;\n"+
             "TRACK, ONEPASS, APERTURE, RECLOSS")
    if settings.trackendonly is None and settings.savetracks:
        line += (", FILE='"+str(k)+"/track.batch"+str(k)+"';\n\n")
    else:
        line += (";\n\n")

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
                 "MAXAPER={0.5,0.05,0.5,0.05,"+str(0.03*settings.nturns)+",0.1}, "+
                 "FFILE="+str(settings.ffile)+";\n\n"+

             "ENDTRACK;\n"+
             "OPTION, WARN;\n\n")

    if settings.trackendonly is not None and settings.savetracks:
        line+="CALL, FILE='"+settings.home+"/madx/end_only.cmdx';\n"

        for obsnum in range(len(settings.elements)+1):
            for partnum in range(settings.nparperbatch):
                temp = "obs"+str(obsnum+1).zfill(4)+".p"+str(partnum+1).zfill(4)
                line+="EXEC, end_only("+str(settings.trackendonly)+", track."+temp+", '"+str(k)+"/track.batch"+str(k)+"."+temp+"');\n"
        line+="\n"

    if(settings.saveloss):
        line+="WRITE, TABLE = trackloss, FILE = 'losses.tfs';\n\n"

    line += "SYSTEM, 'tar -czf tracks.tar.gz "+str(k)+"';"

    return line


def submit_job(settings):
    """Creates and submits job cluster for simulation defined by settings.

    Needs some fixing:
    - We may want to pick non-random seeds for comparing loss-reduction
      methods.
    - The adapter needs to be made much more flexible, and/or different
      template options added.
    - The generated twiss table should be generated in output, not input.
    """
    starting_dir = os.getcwd()

    if settings.slices is not None:
        nslices = len(settings.slices)
        if settings.trackerrep==track_lin:
            settings.trackerrep = track_sliced
    else:
        nslices = 1
            

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

    with open("settings.info", 'w') as f:
        f.write("These files were generated using madxBatch version "+__version__+".\n\n")
        f.write("Settings used:\n")
        for item in vars(settings).items():
            f.write(str(item[0])+" = "+str(item[1])+"\n")

    #Fill in basic template data
    copyfile(settings.trackertemplate, "tracker.madx")

    if settings.thickchanges is not None:
        with open (settings.thickchanges, "r") as changefile:
            changes=changefile.read()
        replacer("tracker.madx", "/*pyTHICKCHANGES*/", changes)
    if settings.thinchanges is not None:
        with open (settings.thinchanges, "r") as changefile:
            changes=changefile.read()
        replacer("tracker.madx", "/*pyTHINCHANGES*/", changes)
    if settings.finalchanges is not None:
        with open (settings.finalchanges, "r") as changefile:
            changes=changefile.read()
        replacer("tracker.madx", "/*pyFINALCHANGES*/", changes)

    replacer("tracker.madx", 'pyTUNESETUP', tune_setup(settings))

    replacer("tracker.madx", 'pyDATADIR', settings.datadir)
    replacer("tracker.madx", 'pyHOMEDIR', settings.home)
    replacer("tracker.madx", 'pyPYCOLL', str(int(settings.pycollimate)))

    for key, replacement in settings.myreplace.iteritems():
        replacer("tracker.madx", key, replacement)

    #Generate twiss files before and after thinning, used to make initial distributions
    if settings.twissfile is None:
        log = open("twisslog.txt", 'w')
        copyfile("tracker.madx", "table.madx")
        replacer("table.madx", 'pyNTURNS', str(settings.nturns))
        replacer("table.madx", '/*pyTWISS', '')
        replacer("table.madx", 'pyTWISS*/', '')
        print "Creating Twiss table"
        if settings.pycollimate:
            subprocess.Popen(settings.pycolldir+"madxColl<table.madx", stdout=log, shell=True).wait()
        else:
            subprocess.Popen(settings.madxversion+"<table.madx", stdout=log, shell=True).wait()
        log.close()
        os.remove("table.madx")
        print 'Table creation is finished!'
    else:
        copyfile(settings.twissfile, settings.datadir+"/thin_twiss.tfs")

    if(settings.trackingbool):
        #Generate initial particle distribution
        #TODO: Non-random input
        if settings.slices is None:
            data=dis.get_gauss_distribution(output=settings.datadir+'initial_distribution.txt',
                                            twissfile=settings.datadir+"thin_twiss.tfs",
                                            n_part=settings.nbatches*settings.nparperbatch,
                                            beam_t=settings.beam_t, seed=settings.seed,
                                            dpp_d=settings.dppmax, **settings.beamkwargs)
        else:
            slicewidth = 0 if settings.slicewidth is None else settings.slicewidth
            data=dis.get_sliced_distribution(output=settings.datadir+'initial_distribution.txt',
                                             twissfile=settings.datadir+"thin_twiss.tfs",
                                             n_part=settings.nparperbatch,
                                             n_batches=settings.nbatches,
                                             beam_t=settings.beam_t, seed=settings.seed,
                                             dpps=settings.slices, dpp_d=slicewidth,
                                             **settings.beamkwargs)
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

        #Change from scattering to black in pycollimate, if needed
        if settings.pcblack:
            copyfile(settings.pycolldir+"track_inside_coll.py", "track_inside_coll.py")
            replacer("track_inside_coll.py", "black=False", "black=True")

        #Prepare job files for each batch
        copyfile(settings.home+"/other/job.sh", "jobs/start.sh")
        if settings.pycollimate:
            replacer("jobs/start.sh", 'MADXEXE', "madxColl")
        else:
            replacer("jobs/start.sh", 'MADXEXE', settings.madxversion.split('/')[-1])
        for i in range(settings.nbatches*nslices):
            outfile="jobs/"+str(i)+".madx"
            searchExp='pyTRACKER'
            replaceExp=settings.trackerrep(i,data,settings)
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
                              ("" if (settings.pcblack) else settings.pycolldir)+"track_inside_coll.py, "+
                              settings.pycolldir+"pycollimate.py, "+
                              (settings.septadb if (settings.septadbreplace is None) else "septa_DB_custom.tfs")+", "+
                              settings.home+"/other/matplotlibrc"+"\n")#, "+
                              #settings.slowexfiles+"\n")
            else:
                subfile.write("transfer_input_files = "+
                              "jobs/$(ProcId).madx, "+
                              settings.madxversion+"\n")#", "+
                              #settings.slowexfiles+"\n")
            subfile.write('arguments = "$(ProcId)"\n')
            subfile.write("initialdir = "+settings.datadir+"\n")
            if settings.saveout:
                subfile.write("output = output/$(ProcId).out\n")
            subfile.write("error = error/$(ProcId).err\n")
            subfile.write("log = log.txt\n")
            subfile.write('transfer_output_remaps = "'+
                          ('tracks.tar.gz=tracks/$(ProcId).tar.gz; ' if settings.savetracks else '')+
                          'losses.tfs=losses/$(ProcId).tfs"\n')
            if settings.flavour is None:
                subfile.write('+JobFlavour = "'+flavour(settings.nturns, settings.nparperbatch, settings.pycollimate)+'"\n')
            else:
                subfile.write('+JobFlavour = "'+settings.flavour+'"\n')
            subfile.write("\nqueue "+str(settings.nbatches*nslices))
        print 'Submit file created!'

        #Submit job
        if settings.local:
            localsub(settings.name+".sub")
            print 'Local batching completed!'
        else:
            subprocess.check_call("condor_submit "+settings.name+".sub", shell=True)
            if settings.monitor:
                subprocess.check_call("condor_wait -status "+settings.datadir+"log.txt", shell=True)

    os.chdir(starting_dir)


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
