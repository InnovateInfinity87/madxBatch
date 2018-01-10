from __future__ import division, print_function
import sys, re
import numpy as np
import pandas as pd
basePath = '/afs/cern.ch/project/sloex/'
sys.path.append(basePath+'code/madxBatch')
sys.path.append(basePath+'code/slowExtractionMADX')
from python.batching import Settings, submit_job, track_lin

## ------ Global parameters ------ ##
# TODO: gamma and magnet strengths should not be hardcoded
c = 299792458.0
gamma = 426.3167458
beta = (1.0-1.0/(gamma**2))**(.5)
v = beta*c
L = 6911.5038 # SPS length
turns_per_sec = v/L # revolution frequency: 43376 (~ kHz)

# bump dipoles
kMP = np.array([7.6765e-5, -.49e-3, -.33309e-3, -.2503e-3, .35585e-3])/4.
kMP = dict(zip(['21202', '21431', '21732', '21995', '22195'], kMP))

# bending dipoles
kMB = 8.445141542E-03 # from /afs/cern.ch/eng/sps/2014/strength/ft_noqs_ext.str
# ------------------ #

class Ripple(object):
    """
    Study the effect of power converter ripples on different magnets during slow extraction.

    @param rippleParams: 3-by-N pandas DataFrame containing the ripple parameters -- amplitude (A), frequency (f) and phase (p) -- for N types of magnets. The supported magnet types are those listed in self.elements. The dataframe columns have to be named accordingly.
    @param baseTracker: MAD-X tracking code used as template.

    """

    parentDirectory = '/'.join((basePath,'ripple')).replace('//', '/')

    def __init__(self, rippleParams, baseTracker = track_lin):
        # Tracker template
        self.baseTracker = baseTracker
        # List of magnets that can be rippled with this code - plus RF
        self.elements = ['qf', 'qd', 'qms', 'mb', 'mp', 'lse', 'lsd', 'lsf', 'rf']

        for key in self.elements:
            try:
                setattr(self, key, rippleParams[key])
            except KeyError:
                setattr(self, key, pd.Series([]))

    def getParams(self):
        return pd.DataFrame([getattr(self, e) for e in self.elements],
                            index = self.elements,
                            columns = ['a', 'f', 'p']).T

    def _addRipple(self, element):
        return '{0}*SIN(2*PI*turn/( {1} ) + ( {2} ))'.format(*getattr(self, element).values)

    def _replace(self):
        oldLine = []
        newLine = []

        if not self.qf.empty:
            oldLine.append('kqf1 = m_f * turn + n_f;')
            newLine.append('{0} + 0.5*(kqf1_start + kqf1_end)*( {1} );\n'.format(oldLine[-1], self._addRipple('qf')))

        if not self.qd.empty:
            oldLine.append('kqd = m_d * turn + n_d;')
            newLine.append('{0} + 0.5*(kqd_start + kqd_end)*( {1} );\n'.format(oldLine[-1], self._addRipple('qd')))

        if not self.qms.empty:
            oldLine.append('tr$macro(turn): MACRO = {')
            newLine.append('{0}\nkQMS = {1};\n'.format(oldLine[-1], self._addRipple('qms')))

        if not self.mb.empty:
            line = ''
            line += (' SELECT, FLAG=ERROR, CLEAR;\n'+
                     ' SELECT, FLAG=ERROR, PATTERN="MB.*";\n'+
                     ' EFCOMP, ORDER=0, DKN={};\n').format(self._addRipple('mb'))

            oldLine.append('tr$macro(turn): MACRO = {')
            newLine.append('{0}\n\n{1}'.format(oldLine[-1], line))

        if not self.mp.empty:
            line = ''
            line += (' SELECT, FLAG=ERROR, CLEAR;\n'+
                     ' SELECT, FLAG=ERROR, PATTERN="MPSH_rb\.21202.*";\n'+
                     ' EFCOMP, ORDER=0, DKN={};\n').format(self._addRipple('mp'))
            line += (' SELECT, FLAG=ERROR, CLEAR;\n'+
                     ' SELECT, FLAG=ERROR, PATTERN=""MPSH_rb\.21431.*;\n'+
                     ' EFCOMP, ORDER=0, DKN={};\n').format(self._addRipple('mp'))
            line += (' SELECT, FLAG=ERROR, CLEAR;\n'+
                     ' SELECT, FLAG=ERROR, PATTERN=""MPSH_rb\.21732.*;\n'+
                     ' EFCOMP, ORDER=0, DKN={};\n').format(self._addRipple('mp'))
            line += (' SELECT, FLAG=ERROR, CLEAR;\n'+
                     ' SELECT, FLAG=ERROR, PATTERN=""MPSH_rb\.21995.*;\n'+
                     ' EFCOMP, ORDER=0, DKN={};\n').format(self._addRipple('mp'))
            line += (' SELECT, FLAG=ERROR, CLEAR;\n'+
                     ' SELECT, FLAG=ERROR, PATTERN=""MPSH_rb\.22195.*;\n'+
                     ' EFCOMP, ORDER=0, DKN={};\n').format(self._addRipple('mp'))

            oldLine.append('tr$macro(turn): MACRO = {')
            newLine.append('{0}\n\n{1}'.format(oldLine[-1], line))

        if not self.lse.empty:
            line = ''
            line += 'klse10602_0 = klse10602;\n'
            line += 'klse22402_0 = klse22402;\n'
            line += 'klse40602_0 = klse40602;\n'
            line += 'klse52402_0 = klse52402;\n'

            oldLine.append('n_d = kqd_start - m_d;')
            newLine.append('{0}\n\n{1}'.format(oldLine[-1], line))

            line = ''
            line += ' klse10602 = klse10602_0*(1.0 + ( {} ));\n'.format(self._addRipple('lse'))
            line += ' klse22402 = klse22402_0*(1.0 + ( {} ));\n'.format(self._addRipple('lse'))
            line += ' klse40602 = klse40602_0*(1.0 + ( {} ));\n'.format(self._addRipple('lse'))
            line += ' klse52402 = klse52402_0*(1.0 + ( {} ));\n'.format(self._addRipple('lse'))

            oldLine.append('tr$macro(turn): MACRO = {')
            newLine.append('{0}\n\n{1}'.format(oldLine[-1], line))

        if not self.lsd.empty:
            line = ''
            line += 'klsda_0 = klsda;\n'
            line += 'klsdb_0 = klsdb;\n'

            oldLine.append('n_d = kqd_start - m_d;')
            newLine.append('{0}\n\n{1}'.format(oldLine[-1], line))
            line = ''
            line += ' klsda = klsda_0*(1.0 + ( {} ));\n'.format(self._addRipple('lsd'))
            line += ' klsdb = klsdb_0*(1.0 + ( {} ));\n'.format(self._addRipple('lsd'))

            oldLine.append('tr$macro(turn): MACRO = {')
            newLine.append('{0}\n\n{1}'.format(oldLine[-1], line))

        if not self.lsf.empty:
            line = ''
            line += 'klsfa_0 = klsfa;\n'
            line += 'klsfb_0 = klsfb;\n'

            oldLine.append('n_d = kqd_start - m_d;')
            newLine.append('{0}\n\n{1}'.format(oldLine[-1], line))

            line = ''
            line += ' klsfa = klsfa_0*(1.0 + ( {} ));\n'.format(self._addRipple('lsf'))
            line += ' klsfb = klsfb_0*(1.0 + ( {} ));\n'.format(self._addRipple('lsf'))

            oldLine.append('tr$macro(turn): MACRO = {')
            newLine.append('{0}\n\n{1}'.format(oldLine[-1], line))

        if not self.rf.empty:
            line = ('ACTA_new.31637: RFCAVITY, L = 17.782, VOLT = {0}, FREQ = {1}, LAG = {2};\n'+
                    'SEQEDIT, SEQUENCE = sps;\n'+
                    'FLATTEN;\n'+
                    'REPLACE, ELEMENT = ACTA.31637, BY = ACTA_new.31637;\n'+
                    'FLATTEN;\n'+
                    'ENDEDIT;\n'+
                    'USE, SEQUENCE = sps;\n\n')

            oldLine.append('/*pyTHICKCHANGES*/')
            newLine.append(line.format(*self.rf))

        return zip(oldLine, newLine)

    def tracker(self, k, data, settings):
        madxCode = self.baseTracker(k, data, settings)
        lines = self._replace()

        for oldLine, newLine in lines:
            madxCode = madxCode.replace(oldLine, newLine)

        return madxCode

    def _printTracker(self):
        """Only for checking everything works"""
        madxCode = self.tracker(0, np.random.randn(200, 200), Settings('dummy', 'dum'))

        print(madxCode)


class Align(object):

    parentDirectory = '/'.join((basePath, 'emittance/zsalign')).replace('//', '/')

    def __init__(self, alignParams, baseTracker = track_lin, trackerTemplate = 'tracker_nominal_template.madx'):
        self.baseTracker = baseTracker
        self.trackerTemplate = trackerTemplate
        self.alignParams = alignParams

        for key in alignParams.keys():
            setattr(self, key, alignParams[key])

    def _replace(self):
        line = '\n'

        for var, value in zip(self.alignParams.index, self.alignParams):
            line += '{} = {};\n'.format(var, value)

        return line

    def tracker(self, k, data, settings):
        """Change ZS alignment in tracker_nominal_template.madx.

        The function works by side-effect, changing the template in settings.trackertemplate.
        @return: None"""
        try:
            with open(settings.trackertemplate, 'r') as f:
                madxCode = f.read()
        except IOError:
            print('Error reading tracker template at {}'.format(settings.trackertemplate))
            return

        # Remove zswirethickness before adding it again
        madxCode = re.sub('\nzswirethickness.*=.*\n', '\n', madxCode)
        # Pattern to replace
        regexp = re.compile(r'\nzswireup.*=.*\nzswiredo.*=.*\n')
        # Replacement
        line = self._replace()

        if regexp.search(madxCode):
            madxCode = re.sub('\nzswireup.*\nzswiredo.*\n', line, madxCode)
        else:
            print('WARNING: no substitutions were made because the original tracker {} did not define girder alignment'.format(self.settings.trackertemplate))
            return

        newTrackerTemplate = '/'.join((self.parentDirectory, 'tracker_new_template.madx')).replace('//', '/')

        with open(newTrackerTemplate, 'w') as f:
            f.write(madxCode)

        settings.trackertemplate = newTrackerTemplate
        return self.baseTracker(k, data, settings)

    def _printTracker(self):
        """Only for checking everything works"""
        madxCode = self.tracker(0, np.random.randn(200, 200), Settings('dummy', 'dum'))

        print(madxCode)

class Scan(object):
    """Create and run several studies of the same type for different values of the parameters.

    TODO: list public methods and instance vars

    This class defines the prototype scan. It takes a settings object as input and uses it to run the simulations. The runScan method takes as input a csv file. Each row of the file must contain a certain combination of the relevant parameters (e.g. amplitude, frequency and phase for ripple studies). The method runs a simulation for each row in the file by submitting the settings object to condor manager.


    """

    def __init__(self, studyClass, paramsFile):
        """Initialise scan instance.

        @param studyClass: Study class - e.g. Ripple or Alignment.
        @param paramsFile: path to csv file containing one study (i.e. one combination of parameter values) in each row.
        """
        self.studyClass = studyClass
        self.paramsFile = paramsFile

    def paramsToFolder(self, params):
        """Create study's folder name from the parameter's values."""
        #TODO: check that if you're only rippling e.g. the QF this doesn't return QF_A_16.0_f_50.0_p_0.0_{a shitton of useless zeroes that tell you you're not rippling anything else}
        return '_'.join([str(p) for p in np.array([params.index, rippleParams]).T.flatten()])

    def listLossFiles(self):
        """Count number of loss tfs files in each study.

        @return: dict with {study: number of loss files}
        """
        d = pd.read_csv(self.paramsFile)

        folders = map(self.paramsToFolder, [row[1] for row in d.iterrows()])
        paths = ('/'.join(basePath, studyGroup, folder, 'losses').replace('//', '/') for folder in folders)

        lossCount = []

        for path in paths:
            try:
                lossCount.append(len(os.listdir(path)))
            except OSError:
                print('Folder {} does not exist'.format(path))
                return

        return dict(zip(folders, lossCount))

    def checkJobs(self):
        """[NOT FINISHED!!!] Check if there are any jobs pending. If not, check if any failed. If so, resubmit. Otherwise, start post-processing."""
        if not self.settings:
            print('Scan object has no settings attribute - you have to run the scan before checking for jobs!')
            return

        if getJobs()['jobs']:
            print('Some jobs are still running...')
            print(getJobs())
            return
        else:
            lossFiles = self.listLossFiles()

            for lossFile in lossFiles:
                if lossFiles[lossFile] < self.settings.nbatches:
                    pass
                    # resubmit
                else:
                    pass
                    # postprocess
        return

    def run(self, nturns = 50000, nbatches = 500, nparperbatch = 200, seed = 0, pycollimate = False, elements = ['AP.UP.ZS21633','AP.DO.ZS21676','AP.UP.TPST21760'], ffile = 500, savetracks = False):

        # Load parameters file
        d = pd.read_csv(self.paramsFile)

        for _, row in d.iterrows():
            # Create instance of study for these combination of parameters
            study = studyClass(row)

            # Create settings instance
            self.settings = Settings(self.paramsToFolder(row), studygroup = self.studyClass.parentDirectory, disk = 'afsproject')

            # Set folder name for the study in the studyGroup directory
            self.settings.name = self.paramsToFolder(row)

            # Set tracker to custom
            # TODO: re-write Alignment tracker so that this makes sense
            self.settings.trackerrep = lambda k, data, settings: self.study.tracker(k, data, settings)

            # Set random seed
            self.settings.seed = seed
            # Use pyCollimate scattering routine
            self.settings.pycollimate = pycollimate
            # Observe losses at these elements
            self.settings.elements = elements
            # Simulation length in number of turns
            self.settings.nturns = nturns
            self.settings.ffile = ffile # ffile MADX ?

            # Number of batches
            self.settings.nbatches = nbatches # number of batches
            # Number of particles per batch
            self.settings.nparperbatch = nparperbatch # particles per batch

            # Save track files - takes up a lot of space
            self.settings.savetracks = savetracks

            submit_job(self.settings)

    # Post-processing goes here. It could mostly be done by calling functions from the dataprocessing module and looping:
    # 1. tfs files to csv including where the particle hit
    # 2. Computing summary loss statistics
    # 3. Spill quality analysis - spill and FFT plots, duty factor, etc.
    # 4. Contour plots for losses, emittance and extraction efficiency
    # 5. Other plots ?


def parseLossFiles(scanFolder, lossFolders, align):
    """Reads loss files for every study in a scan and finds where particles hit the ZS.

    @param lossFolders: list of lists of paths to tfs files. lossFolders[i][j] contains the j-th tfs loss file for the i-th study in a scan.
    @param align: pd.DataFrame containing alignment information.
    @return: None - saves csv files (one per study) with (X,PX) coordinates at ZS upstream, at the point the particle hit the ZS and label indicating tank where it hit and from where.
    """
    hDir = '/'.join((scanFolder+ '_hits')).replace('//','/')
    if not os.path.isdir(hDir):
        os.makedirs(hDir)

    for folder in lossFolders:
        h = []
        hFile = '/'.join(hDir, folder[0].split('/')[-3], )).replace('//','/')

        if os.path.isfile(hFile):
            print('{} already exists'.format(hFile))
        else:
            print('{} needs to be computed...'.format(hFile))

            for f in folder:
                try:
                    _, table = dp.readtfs(f)
                except UnboundLocalError:
                    print('File {} was empty!'.format(f))
                    continue

                for _, particle in table.iterrows():
                    particle = fieldFreeBacktrack(particle)
                    h.append(checkIntersection(particle, align))

                pd.DataFrame(hits).to_csv(hFile)

def getJobs():
    """Call condor_q, read output and save it to a dict with number of jobs in each state (running, idle, etc.)"""

    condor = os.popen('condor_q').read()

    if condor:
        try:
            jobs = filter(lambda x: 'jobs' in x, condor.split('\n'))[0].replace(';', ',').split(',')
        except IndexError:
            print('No jobs found in condor_q output')
            return
    else:
        print('condor_q not found - are you on lxplus?')
        return

    return {i[1]: int(i[0]) for i in [j.split() for j in jobs]}


def checkIntersection(particle0, align, kick = 4.1635e-4):
    # TODO: you wrote this with the anode in mind so re-using the function for the cathode is gonna be tricky - the particle cannot come from the outside, the field is on the other side, etc. Try to fix that!
    """Calculate whether/where a particle hits a set of wires.

    @param particle0: dictionary/pd.Series containing at least X and PX coordinates of the particle at ZS upstream.
    @param align: pd.DataFrame containing girder and wires alignment data. Each row corresponds to one of the tanks or the whole girder.
                  columns: s - element's upstream end s coordinate
                           up, do - upstream and downstream (mis)alignments
                           l - element length
                           thick - wire thickness
    @param kick: ZS deflection angle
    @return: final (x,p) coordinates and label indicating which tank was hit and from where (outside or inside). Label is 'EXT.IN' if the particle was extracted.
    """

    # Slope and offset of the wires calculated from girder position and misalignment
    align['theta'] = (align['do'] - align['up'])/align['l']
    align['gamma'] = align['up']['girder'] + align['up'] - align['theta']*align['s']

    # ZS angle set by girder
    th = align['theta']['girder']
    # ZS length
    L = align['l']['girder']
    # ZS initial s
    s0 = align['s']['tank1']
    # Position of the outer and inner sides of the first wire
    uIn = align['up']['girder'] + align['up']['tank1']
    uOut = uIn - align['thick']['tank1']

    x0, p0 = map(lambda k: particle0[k], ('X', 'PX'))

    # Hits the wires before entering the ZS
    if x0 > uOut and x0 < uIn:
        return {'X0': x0, 'P0': p0, 'XF': x0, 'PF': p0, 'HIT': 'tank1.HO'}

    # Kick only if inside ZS, otherwise field free region
    a = kick if x0 > uIn else 0

    # TODO: check that this loop always goes from tank 1 to 5 - sort align df if needed
    for k, wire in align.iterrows():
        if 'tank' in wire.name:

            # TODO: Check if it's hitting head on first

            delta = (p0-th-wire['theta'])**2 - 2*a*(x0-wire['gamma'])/L

            if delta > 0:
                if a > 0:
                    sf = s0 + L*(th+wire['theta']-p0-np.sqrt(delta))/a
                else:
                    try:
                        sf = s0 + (wire['gamma']-x0)/(p0-th-wire['theta'])
                    except ZeroDivisionError:
                        print('Trajectory parallel to the ZS wires!')
                        sf = np.nan
                        break

                tankUp, tankDo = wire['s'], align['s'].shift(-1)[k]

                if sf > tankUp and sf < tankDo:
                    hit = '{tank}.{side}'.format(tank = k, side = 'OUT' if a == 0 else 'IN')
                    break
    else:
        sf = s0 + L
        hit = 'EXT.IN'

    xf = a*(sf-s0)**2/(2*L) + p0*(sf-s0) + x0
    pf = a*(sf-s0)/L + p0

    return {'X0': x0, 'P0': p0, 'XF': xf, 'PF': pf, 'HIT': hit}

def fieldFreeBacktrack(particle, s0 = 1668.9775):
    particle['X'] += -particle['PX']*(particle['S'] - s0)
    particle['Y'] += -particle['PY']*(particle['S'] - s0)

    return particle

def trackZS(particle, up = None, do = None):
    #TODO: think how to pass the alignment info and make it reasonably compatible with the class defined above
    """Tracks particle lost at (S, X, PX, Y, PY) along the ZS, given certain alignment.

    @params particle: dictionary containing the coordinates (s, x, px, y, py) at which the particle was lost by the tracking routine as well as the name of the element corresponding to this s.
    @params **kwds: dictionary containing the ZS upstream (zswireup) and downstream (zswiredo) positions, the misalignment of each cathode (up1, do1, ..., up5, do5) and the wire thickness (zswirethickness). If a parameter is not provided it defaults to its nominal value (e.g. zero for misalignments).
    @return: particle coordinates upstream of the ZS, final particle coordinates (either hitting or exiting the ZS) and labels indicating whether the particle was lost/extracted and, in case of loss, where it happened.
    """

    # Constant ZS parameters
    s0 = 1668.9776 # ZS upstream, s coordinate
    L = 18.77 # ZS length
    n = 5 # number of tanks
    l = 3.13 # tank length
    g = .78 # gap between tanks
    c = .02 # gap between cathode and anode
    apery = .023 # vertical aperture, symmetric

    particle0 = fieldFreeBacktrack(particle)

    tanks = [i*(l+g) for i in range(n)] + [L]

    align = pd.DataFrame(np.array([tanks, up, do]), columns = ['s', 'up', 'do'])
    align.index = ['tank{}'.format(i) for i in range(n)]+['girder']
    align['l'] = l
    align['l']['girder'] = L
