from __future__ import print_function
from __future__ import division
import sys, os, time, re, copy

sys.dont_write_bytecode = True

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

BASEPATH = '/afs/cern.ch/project/sloex/'

from python.batching import Settings, submit_job, track_lin
import python.dataprocessing as dp

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

    def __init__(self, rippleParams, baseTracker = track_lin):
        # Tracker template
        self.baseTracker = baseTracker
        # List of magnets that can be rippled with this code - plus RF
        self.elements = ['qf', 'qd', 'qms', 'mb', 'mp', 'lse', 'lsd', 'lsf', 'rf']

        self.folder = ''

        for key in self.elements:
            try:
                setattr(self, key, rippleParams[key])
                self.folder += '_'.join([key] + map(str, rippleParams[key].tolist()))
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

    def run(self, studyname = None, studygroup = 'ripple', nturns = 50000, nbatches = 500, nparperbatch = 200, seed = 0, pycollimate = False, elements = ['AP.UP.ZS21633','AP.DO.ZS21676','AP.UP.TPST21760'], ffile = 500, savetracks = False, **mySettings):

        # Create settings instance
        if studyname is None:
            studyname = self.folder
        self.settings = Settings(self.folder, studygroup = studygroup, disk = 'afsproject')

        # Set folder name for the study in the studyGroup directory
        self.settings.name = self.folder

        # Set tracker to custom
        self.settings.trackerrep = lambda k, data, settings: self.tracker(k, data, settings)

        # Set random seed
        self.settings.seed = seed
        # Use pyCollimate scattering routine
        self.settings.pycollimate = pycollimate
        # Observe losses at these elements
        self.settings.elements = elements
        # Simulation length in number of turns
        self.settings.nturns = nturns
        # Save tracks every ffile turns
        self.settings.ffile = ffile

        # Number of batches
        self.settings.nbatches = nbatches # number of batches
        # Number of particles per batch
        self.settings.nparperbatch = nparperbatch # particles per batch

        # Save track files - takes up a lot of space
        self.settings.savetracks = savetracks

        # Other settings
        for mySetting in mySettings:
            setattr(self.settings, mySetting, mySettings[mySetting])

        submit_job(self.settings)


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
        paths = ('/'.join(BASEPATH, studyGroup, folder, 'losses').replace('//', '/') for folder in folders)

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

    def runScan(self, nturns = 50000, nbatches = 500, nparperbatch = 200, seed = 0, pycollimate = False, elements = ['AP.UP.ZS21633','AP.DO.ZS21676','AP.UP.TPST21760'], ffile = 500, savetracks = False):

        # Load parameters file
        d = pd.read_csv(self.paramsFile)

        for _, row in d.iterrows():
            # Create instance of study for these combination of parameters
            study = studyClass(row)

            study.run(nturns = nturns, nbatches = nbatches, nparperbatch = nparperbatch, seed = seed, pycollimate = pycollimate, elements = elements, ffile = ffile, savetracks = savetracks)


class Align(object):

    parentDirectory = '/'.join((BASEPATH, 'emittance/zsalign')).replace('//', '/')

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


    # Post-processing goes here. It could mostly be done by calling functions from the dataprocessing module and looping:
    # 1. tfs files to csv including where the particle hit
    # 2. Computing summary loss statistics
    # 3. Spill quality analysis - spill and FFT plots, duty factor, etc.
    # 4. Contour plots for losses, emittance and extraction efficiency
    # 5. Other plots ?


def parseLossFiles(lossFiles, align, nParPerBatch = 200, doPrint = True):

    h = []

    for f in lossFiles:
        nBatch = int(f.split('/')[-1].split('.')[0])

        if doPrint:
            print('Reading... {}.tfs'.format(nBatch))

        try:
            _, table = dp.readtfs(f)
        except UnboundLocalError:
            print('File {} was empty!'.format(f))
            continue

        for _, particle in table.iterrows():
            # Do the thing
            hit = zsTrack(particle, align)

            # Add batch number to particle id
            hit['NUMBER'] += nParPerBatch*nBatch

            h.append(hit)

    return pd.DataFrame(h)


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











## TODO ##
# Probably remove everythin below here since I've moved it to another repo #




## Cross and checkCross where meant to be used with track files ##

def cross(d, wire, girder, kick = 4.1635e-4):
    """Check if particle hits a wire in this order: HO, IN, OUT, SLIP."""

    # TODO: Add batch number to this
    number = d.index[0]

    zsup, zsdo, th = girder.up, girder.do, girder.theta

    tank = int(wire.name[-1])

    upOut = wire.up + girder.up
    upIn = upOut + wire.thick

    doOut = wire.do + girder.up
    doIn = doOut + wire.thick

    # Head on
    ho = (d.X > upOut)*(d.PX < upIn)

    if ho.any():
        hit = d[ho].iloc[0,:]
        return {'NUMBER': number,
                'X': hit.X, 'PX': hit.PX, 'S_': hit.S,
                'HIT': 'HO', 'tank': tank}

    # Inside
    din = d[d.X > upIn]
    delta = (din.PX - th - wire.theta)**2 - 2*kick*(din.X - zsup - wire.gamma)/L
    sf = din.S + L*(th + wire.theta - din.PX - np.sqrt(delta))/kick

    ins = (sf > wire.s)*(sf < wire.s + wire.l)

    if ins.any():
        hit = din[ins].iloc[0,:]
        return {'NUMBER': number,
                'X': hit.X, 'PX': hit.P, 'S_': sf[ins][0],
                'HIT': 'IN', 'tank': tank}

    # Outside
    dout = d[d.X < upOut]
    try:
        sf = dout.S + (zsup+wire.gamma-dout.X)/(dout.PX-th-wire.theta)
    except ZeroDivisionError:
        raise

    out = (sf > wire.s)*(sf < wire.s + wire.l)
    if out.any():
        hit = dout[out].iloc[0,:]
        return {'NUMBER': number,
                'X': hit.X, 'PX': hit.PX, 'S_': sf[out].iloc[0],
                'HIT': 'OUT', 'tank': tank}

    # Enters ZS through tank gaps
    slip = (sf > wire.s + wire.l)*(sf < wire.s + wire.l + wire.gap)
    if slip.any():
        hit = dout[slip].iloc[0,:]
        return {'NUMBER': number,
                'X': hit.X, 'PX': hit.PX, 'S_': sf[slip].iloc[0],
                'HIT': 'SLIP', 'tank': tank}

    return {}

def checkCross(d, align, kick = 4.1635e-4):
    """Call cross function iteratively for every tank."""

    d = fieldFreeBacktrack(d)

    align['theta'] = (align.do - align.up)/align.l
    align['gamma'] = align.up - align.theta*align.s
    align['gap'] = np.diff(align.s[:-1]) - align.l[:-2]

    girder = align.loc['girder', :]
    h = {}

    for k, wire in align.iterrows():
        if 'tank' in wire.name and not h:
            h = cross(d, wire, girder)
    else:
        hit = d.iloc[-1,:]
        h = {'NUMBER': d.index[0],
             'X': hit.X, 'PX': hit.PX, 'S_': hit.S,
             'HIT': 'EXT', 'tank': np.nan}

    return h


### As of 19/1/18 checkIntersection seems to work fine and checkLosses is an attempt at vectorising checkIntersection that has yet to be benchmarked

def checkLosses(_d, align, kick = 4.1635e-4):
    """Compute losses
    @param _d: pd.DataFrame - loss table
    @param align: pd.DataFrame - alignment settings
    @param kick: float - ZS kick"""

    align['theta'] = (align.do - align.up)/align.l
    align['gamma'] = align.up - align.theta*align.s
    align['gap'] = np.diff(align.s[:-1]) - align.l[:-2]

    d = _d.copy()

    d = advance(d, min(d.S) - d.S, kick = 0)

    d['X0'] = d.X
    d['PX0'] = d.PX
    d['S0'] = d.S

    d['HIT'] = '_OUT'
    d['tank'] = 'none'

    for tank, wire in align.iterrows():
        if 'tank' in wire.name:
            uOut = align.up.girder + wire.up + align.theta.girder*(wire.s - align.s.girder)
            uIn = uOut + wire.thick

            # Check if the particles hit the wires head on or move into/out of the ZS
            d.pipe(moveIn, uOut)
            d.pipe(moveOut, uIn)
            d.pipe(hitHO, uIn, uOut, tank)

            # Move particles along tank, checking if they're lost or extracted
            d.loc[d.HIT == '_IN', :] = evolve(d.loc[d.HIT == '_IN', :], align, tank, kick)
            d.loc[d.HIT == '_OUT', :] = evolve(d.loc[d.HIT == '_OUT', :], align, tank, kick = 0)

    d.pipe(declareCirculating)

    return d


def evolve(d, align, tank, kick):
    """Check if particles are lost, extracted or keep circulating."""

    # Global alignment parameters
    th = align.theta.girder
    L = align.l.girder
    s0 = align.s.girder
    zsUp, zsDo = align.loc['girder', ['up', 'do']]

    # Wire parameters
    wire = align.loc[tank, :]
    theta = th + wire.theta
    sUp, sDo = wire.s, wire.s + wire.l
    xUp, xDo = zsUp + th*sUp + wire.up, zsUp + th*sUp + wire.up + wire.theta*(sDo - sUp)

    l = wire.l if tank == 'tank5' else wire.l + wire.gap

    # Coordinates
    x0, p0 = d.X, d.PX

    # Discriminant for the intersection equation
    delta = (p0 - theta)**2 - 2*kick*(x0 - xUp)/L

    if kick > 0:
        hit = 'IN'
        sHit = sUp + L*(theta - p0 + np.sqrt(delta))/kick
       # TODO: what about cathode?
    else:
        hit = 'OUT'
        try:
            sHit = sUp + (xUp - x0)/(p0 - theta)
        except ZeroDivisionError:
            pass

    # Extracted
    d.loc[(delta < 0), :] = advance(d.loc[(delta < 0), :], L - sUp, kick, 'EXT')

    # Lost
    d.loc[(delta > 0) & (sHit > sUp) & (sHit < sDo), :] = advance(d.loc[(delta > 0) & (sHit > sUp) & (sHit < sDo), :], sHit - sUp, kick, hit, tank)

    # Circulating
    d.loc[(delta > 0) & ((sHit < sUp) | (sHit > sDo)), :] = advance(d.loc[(delta > 0) & ((sHit < sUp) | (sHit > sDo)), :], l, kick)

    return d


def advance(d, l, kick, hit = None, tank = None):
    """Evolve particle trajectory by distance l and acceleration kick and update if they are lost."""
    d.S += l
    d.X += d.PX*l + kick*l**2
    d.PX += 2*kick*l

    if hit is not None:
        d.HIT = hit

    if tank is not None:
        d.tank = tank

    return d


def fieldFreeBacktrack(d, align):
    #TODO: backtrack in (T, PT) as well - you will only need compaction factor and beam's gamma if RF is off, but it gets tricky if RF is on
    # Do we really need this if we only care about the longitudinal space for extracted particles?
    """Backtrack particles lost along the ZS to the position of the first tank."""
    s0 = align.s.tank1

    d.X -= d.PX*(d.S - s0)
    d.Y -= d.PY*(d.S - s0)
    d.S = s0

    return d


def moveIn(d, uIn):
    d.loc[(d.HIT == '_OUT') & (d.X > uIn), 'HIT'] = '_IN'

def moveOut(d, uOut):
    d.loc[(d.HIT == '_IN') & (d.X < uOut), 'HIT'] = '_OUT'

def hitHO(d, uIn, uOut, tank):
    d.loc[(d.HIT.isin(('_IN', '_OUT'))) & (d.X < uIn) & (d.X > uOut), ['HIT', 'tank']] = 'HO', tank

def declareCirculating(d):
    d.loc[d.HIT.apply(lambda x: x.startswith('_')), 'HIT'] = 'CIRC'

def areExtracted(d):
    return d[d.HIT == 'EXT']


def checkIntersection(particle, align, kick = 4.1635e-4):
    # TODO: vectorise this function
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

    # Backtrack - all particles are at ZS upstream now
    particle = fieldFreeBacktrack(particle, align)

    # Slope and offset of the wires calculated from girder position and misalignment
    align['theta'] = (align.do - align.up)/align.l
    align['gamma'] = align.up - align.theta*align.s

    # ZS angle set by girder
    th = align.theta.girder
    # ZS length
    L = align.l.girder
    # ZS initial s
    s0 = align.s.tank1
    # Girder up/downstream
    zsup, zsdo = align.loc['girder', ['up', 'do']]

    # Particle coordinates at ZS.UP
    x0, p0 = map(lambda k: particle[k], ('X', 'PX'))

    # Particle coordinates along ZS
    x, p = x0, p0

    # Kick only if inside ZS, otherwise field free region
    a = kick if x0 > align.up.girder else 0

    # TODO: check that this loop always goes from tank 1 to 5 - sort align df if needed
    for k, wire in align.iterrows():
        if 'tank' in wire.name:
            uOut = align.up.girder + wire.up
            uIn = uOut + wire.thick

            if x > uOut and x < uIn:
                return {'NUMBER': particle.name, 'X': x0, 'PX': p0, 'S_': s0, 'HIT': 'HO', 'tank': 1}

            # Discriminant
            delta = (p0-th-wire.theta)**2 - 2*a*(x0-zsup-wire.gamma)/L

            # Check intersection for every wire
            if delta > 0:
                if a > 0:
                    sf = s0 + L*(th+wire.theta-p0-np.sqrt(delta))/a
                else:
                    try:
                        sf = s0 + (zsup+wire.gamma-x0)/(p0-th-wire.theta)
                    except ZeroDivisionError:
                        continue

                x = a*(sf-s0)**2/(2*L) + p0*(sf-s0) + x0
                p = a*(sf-s0)/L + p0

                #TODO: what happens when a particle slips through the gap between 2 tanks?
                tankUp, tankDo = wire.s, wire.s + wire.l

                if sf > tankUp and sf < tankDo:
                    hit = 'OUT' if a == 0 else 'IN'
                    tank = int(k[-1])
                    break
    else:
        sf = s0 + L
        hit = 'EXT'
        tank = np.nan


    return {'NUMBER': particle.name, 'X': x0, 'PX': p0, 'S_': sf, 'HIT': hit, 'tank': tank}

def zsg(align, s):
    """Get nominal x coordinate of wires at distance s given girder alignment."""
    return align.up.girder + align.theta.girder*(s - align.s.girder)

def move(x, p, s, a):
    """Track particle with initial coordinates (x, p) for distance s in field a"""
    x += p*s + .5*a*s**2
    p += a*s
    return x, p

def zsTrack(part, align, kick = 4.1635e-4):
    """Track particle along the ZS.
    @param part: particle coordinates at ZS upstream
    @param align: pandas DataFrame with ZS alignment settings. One row for each tank plus another for the girder. Columns: s = element's upstream position s coordinate; up, do = element's upstream/downstream x coordinate; thick = wire thickness; l = element's length.
    @param kick: ZS kick
    @return particles coordinates after tracking plus flags indicating if/where it hit the wires (HIT = EXT, IN, OUT, HO) and in which tank"""

    align['theta'] = (align.do - align.up)/align.l
    align['gap'] = np.nan
    align.gap.iloc[:-2] = np.diff(align.s[:-1]) - align.l.iloc[:-2]

    L = align.l.girder

    # Backtrack particles to ZS upstream if they have larger s
    x0, p0 = move(part['X'], part['PX'], align.s.girder - part['S'], a = 0)
    x, p = x0, p0

    for k, wire in align.iterrows():
        if 'tank' in wire.name:

            th = align.theta.girder + wire.theta
            l = wire.l if k == 'tank5' else wire.l + wire.gap

            uOut = zsg(align, wire.s) + wire.up
            uIn = uOut + wire.thick

            dOut = zsg(align, wire.s + wire.l) + wire.up
            dIn = dOut + wire.thick

            # Particle hits wires head on
            if x > uOut and x < uIn:
                sf = wire.s
                hit = 'HO'
                tank = k
                break
            # Particle is outside of the ZS
            elif x < uOut:
                try:
                    sf = wire.s + (x - uOut)/(th - p)
                except ZeroDivisionError:
                    x, p = move(x, p, l, a = 0)
                    continue
            # Particle is inside the ZS
            elif x > uIn:
                delta = (th - p)**2 - 2*kick*(x - uIn)/L
                if delta > 0:
                    sf = wire.s + L*(th - p - np.sqrt(delta))/(2*kick)
                else:
                    x, p = move(x, p, l, a = kick/L)
                    continue
            # Check intersection between trajectory and wires
            if sf > wire.s and sf < wire.s + wire.l:
                hit = 'OUT' if x < uOut else 'IN'
                tank = k
                x, p = move(x, p, sf - wire.s, a = kick/L if x > uIn else 0)
                break
            else:
                x, p = move(x, p, l, a = kick/L if x > uIn else 0)
    else:
        sf = align.s.girder + L
        hit = 'EXT'
        tank = 'none'

    return {
            # TODO: this should work both with pd.Series and dicts
            'NUMBER': part.name,
            'S': align.s.girder, '_S': sf,
            'X': x0, '_X': x,
            'PX': p0, '_PX': p,
            'Y': part['Y'], '_Y': part['Y'] + part['PY']*(sf-align.s.girder),
            'PY': part['PY'], '_PY': part['PY'],
            'T': part['T'], 'PT': part['PT'],
            'HIT': hit, 'tank': tank
            }
