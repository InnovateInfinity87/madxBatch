from version import __version__

import numpy as np
from scipy.stats import truncnorm as trandn
from prettytable import PrettyTable
import linecache
from datetime import datetime

# Comment: In MAD-X dispersion is w.r.t PT, not DELTAP. (see uguide p19-20)
# TODO: set default emit_x, emit_y to realistic values instead of historical?
# TODO: Truncation happensin a parallelogram, not an ellipse atm. Look into implementing truncated Rayleigh?

m_p = 0.938272081 # proton mass in GeV, from PDG website Sept.'17

def header():
    time = datetime.utcnow()
    return ('@ NAME             %20s "INITIAL DISTRIBUTION"\n'+
            '@ TYPE             %04s "USER"\n'+
            '@ TITLE            %20s "INITIAL DISTRIBUTION"\n'+
            '@ ORIGIN           %06s "PYTHON"\n'+
            '@ DATE             %08s "'+time.strftime("%d/%m/%y")+'"\n'+
            '@ TIME             %08s "'+time.strftime("%H.%M.%S")+'"\n')

def string_to_float(seq):
    for x in seq:
        try:
            yield float(x)
        except ValueError:
            yield x


def twissinit(twissfile):
    variables = linecache.getline(twissfile, 46).split()[1:]
    twiss = list(string_to_float(linecache.getline(twissfile, 48).split()))

    x0 = twiss[variables.index('X')]
    y0 = twiss[variables.index('Y')]
    px0 = twiss[variables.index('PX')]
    py0 = twiss[variables.index('PX')]


    betx = twiss[variables.index('BETX')]
    bety = twiss[variables.index('BETY')]

    alfx = twiss[variables.index('ALFX')]
    alfy = twiss[variables.index('ALFY')]

    dx = twiss[variables.index('DX')]
    dpx = twiss[variables.index('DPX')]
    dy = twiss[variables.index('DY')]
    dpy = twiss[variables.index('DPY')]

    return x0,y0,px0,py0,betx,bety,alfx,alfy,dx,dpx,dy,dpy


class Beam(object):
    def __init__(self, beam_type, **kwargs):
        if beam_type == 'FT':
            self.type = 'FT'
            self.set_params(pc=400.0, dpp_0=0.0, dpp_d=0.0015,
                            emit_xn=8E-6, emit_yn=8E-6, n_sigma=6.8,
                            pdist='unif')
        elif beam_type == 'CUSTOM':
            self.type = 'CUSTOM'
        else:
            print('Beam type '+str(beam_type)+' not recognized, '+
                  'creating CUSTOM beam. Parameters should be set '+
                  'manually.')
            self.type = 'CUSTOM'
        self.set_params(**kwargs)

    def set_params(self, pc=None, dpp_0=None, dpp_d=None,
                   emit_xn=None, emit_yn=None, n_sigma=None,
                   pdist=None):
        if pc is not None: self.pc=pc
        if dpp_0 is not None: self.dpp_0 = dpp_0
        if dpp_d is not None: self.dpp_d = dpp_d
        if emit_xn is not None: self.emit_xn = emit_xn
        if emit_yn is not None: self.emit_yn = emit_yn
        if n_sigma is not None: self.n_sigma = n_sigma
        if pdist is not None: self.pdist = pdist
        self._calc_params()

    def _calc_params(self):
        self.energy = np.sqrt(m_p**2 + self.pc**2)
        self.gamma_r = self.energy/m_p
        self.beta_r = np.sqrt(1. - 1./self.gamma_r**2)
        self.emit_x = self.emit_xn/(self.beta_r*self.gamma_r)
        self.emit_y = self.emit_yn/(self.beta_r*self.gamma_r)

    def dpp_to_pt(self, dpp):
        de = np.sqrt(m_p**2 + (self.pc*(1+dpp))**2) - self.energy
        return de/self.pc


def get_gauss_distribution(twissfile='sequence_totrack.tfs', beam_t='FT',
                           sigmas=None, seed=None, n_part=100, noseeding=False,
                           output='initial_distribution.txt', **kwargs):
    x0,y0,px0,py0,betx,bety,alfx,alfy,dx,dpx,dy,dpy = twissinit(twissfile)
    beam = Beam(beam_t, **kwargs)

    if not noseeding:
        np.random.seed(seed)

    if sigmas is None:
        sigmas = beam.n_sigma

    # Momentum distribution
    if beam.pdist == 'unif':
        dpp = np.random.uniform(beam.dpp_0-beam.dpp_d, beam.dpp_0+beam.dpp_d, n_part)
    elif beam.pdist.startswith('gauss'):
        sigp = float(beam.pdist[5:])
        dpp = beam.dpp_0+trandn(-sigp, sigp, scale=(beam.dpp_d/sigp)).rvs(n_part)
    else:
        print "Warning: unknown pdist '"+beam.pdist+"', assuming uniform."
        dpp = np.random.uniform(beam.dpp_0-beam.dpp_d, beam.dpp_0+beam.dpp_d, n_part)

    pt = np.fromiter((beam.dpp_to_pt(d) for d in dpp), float)

    # Transverse distributions
    sx = np.sqrt(beam.emit_x*betx)
    n_x = trandn(-sigmas, sigmas, scale=sx).rvs(n_part)
    n_px = trandn(-sigmas, sigmas, scale=sx).rvs(n_part)
    x = x0 + n_x + dx*pt
    px = px0 + (n_px - alfx*n_x)/betx + dpx*pt

    sy = np.sqrt(beam.emit_y*bety)
    n_y = trandn(-sigmas, sigmas, scale=sy).rvs(n_part)
    n_py = trandn(-sigmas, sigmas, scale=sy).rvs(n_part)
    y = y0 + n_y + dy*pt
    py = py0 + (n_py - alfy*n_y)/bety + dpy*pt

    # Output table
    if output is not None:
        normal = PrettyTable(['*', 'NUMBER', 'TURN', 'X', 'PX', 'Y', 'PY', 'T', 'PT', 'S', 'E'])
        normal.align = 'r'
        normal.left_padding_width = 0
        normal.right_padding_width = 8
        normal.border = False
        normal.add_row(['$', '%d', '%d', '%le', '%le', '%le', '%le', '%le', '%le', '%le', '%le'])
        for i in xrange(n_part):
            normal.add_row([' ', i + 1, 0, x[i], px[i], y[i], py[i], 0.0000000, pt[i], 0.0000000, beam.energy])
        
        with open(output, 'w') as fp:
            fp.write(header())
            fp.write(normal.get_string())

    return x, px, y, py, pt


# either n_halo=number for thin halo, or n_halo=(lower,upper) for fat halo
def get_halo_distribution(twissfile='sequence_totrack.tfs', beam_t='FT',
                          n_halo=5,seed=None, n_part=100,
                          output='initial_distribution_halo.txt', **kwargs):
    x0,y0,px0,py0,betx,bety,alfx,alfy,dx,dpx,dy,dpy = twissinit(twissfile)
    beam = Beam(beam_t, **kwargs)

    try:
        nhalo[0]
    except TypeError:
        nhalo = (nhalo,nhalo)

    np.random.seed(seed)

    # Momentum distribution
    dpp = np.random.uniform(beam.dpp_0-beam.dpp_d, beam.dpp_0+beam.dpp_d, n_part)
    pt = np.fromiter((beam.dpp_to_pt(d) for d in dpp), float)

    # Transverse distributions
    psix = np.random.uniform(0, 2*np.pi, n_part)
    widthx = np.random.uniform(n_halo[0], n_halo[1], n_part)
    x = x0 + widthx*np.sqrt(beam.emit_x*betx)*np.cos(psix) + dx*pt
    px = px0 + widthx*np.sqrt(beam.emit_x/betx)*(np.sin(psix) - alfx*np.cos(psix)) + dpx*pt

    psiy = np.random.uniform(0, 2*np.pi, n_part)
    widthy = np.random.uniform(n_halo[0], n_halo[1], n_part)
    y = y0 + widthy*np.sqrt(beam.emit_y*bety)*np.cos(psiy) + dy*pt
    py = py0 + widthy*np.sqrt(beam.emit_y/bety)*(np.sin(psiy) - alfy*np.cos(psiy)) + dpy*pt

    # Output table
    if output is not None:
        normal = PrettyTable(['*', 'NUMBER', 'TURN', 'X', 'PX', 'Y', 'PY', 'T', 'PT', 'S', 'E'])
        normal.align = 'r'
        normal.left_padding_width = 0
        normal.right_padding_width = 8
        normal.border = False
        normal.add_row(['$', '%d', '%d', '%le', '%le', '%le', '%le', '%le', '%le', '%le', '%le'])
        for i in xrange(n_part):
            normal.add_row([' ', i + 1, 0, x[i], px[i], y[i], py[i], 0.0000000, pt[i], 0.0000000, beam.energy])
        
        with open(output, 'w') as fp:
            fp.write(header())
            fp.write(normal.get_string())

    return x, px, y, py, pt


def get_sliced_distribution(twissfile='sequence_totrack.tfs', beam_t='FT',
                           sigmas=None, dpps=[None],
                           seed=None, n_batches=1, n_part=100,
                           output='initial_distribution_slices.txt', **kwargs):

    np.random.seed(seed)

    xf = np.empty(0)
    pxf = np.empty(0)
    yf = np.empty(0)
    pyf = np.empty(0)
    ptf = np.empty(0)

    energy = Beam(beam_t, **kwargs).energy

    if not 'pdist' in kwargs:
        kwargs['pdist'] = 'gauss1'

    for dpp in dpps:
        x, px, y, py, pt = get_gauss_distribution(twissfile=twissfile, beam_t=beam_t,
                                                  sigmas=sigmas, n_part=(n_batches*n_part),
                                                  output=None, noseeding=True,
                                                  dpp_0=dpp, **kwargs)
        xf = np.concatenate((xf, x))
        pxf = np.concatenate((pxf, px))
        yf = np.concatenate((yf, y))
        pyf = np.concatenate((pyf, py))
        ptf = np.concatenate((ptf, pt))

    # Output table
    if output is not None:
        normal = PrettyTable(['*', 'NUMBER', 'TURN', 'X', 'PX', 'Y', 'PY', 'T', 'PT', 'S', 'E'])
        normal.align = 'r'
        normal.left_padding_width = 0
        normal.right_padding_width = 8
        normal.border = False
        normal.add_row(['$', '%d', '%d', '%le', '%le', '%le', '%le', '%le', '%le', '%le', '%le'])
        for i in xrange(xf.size):
            normal.add_row([' ', i + 1, 0, xf[i], pxf[i], yf[i], pyf[i], 0.0000000, ptf[i], 0.0000000, energy])
        
        with open(output, 'w') as fp:
            fp.write(header())
            fp.write(normal.get_string())

    return xf, pxf, yf, pyf, ptf
