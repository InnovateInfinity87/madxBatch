from __future__ import print_function
from future.utils import iteritems, lrange

import fileinput
import matplotlib as mpl
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import os
import pandas as pd
from scipy.stats import norm
import seaborn as sns
import shutil
import sys
import tarfile

# Warn if using old matplotlib
_mplver = [int(x) for x in mpl.__version__.split('.')]
if _mplver[0]<2 or (_mplver[0]==2 and _mplver[1]==0 and _mplver[2]<2):
    print('WARNING: dataprocessing needs matplotlib version 2.0.2 or higher to function properly.')

# Set style
sns.set()
plt.style.use('seaborn-ticks')
sns.set_context("notebook", rc={"text.usetex": True})
sns.set_style({'image.cmap': 'viridis'})

#TODO python2/3 compatibility
#TODO look for plotting related stuff in todo thingy
#TODO Some more specifics down in the code
#TODO Check where to use check_lossbug=False for readtfs (tracks!)
#TODO revive old movie code?

# File contents:
# 1. Utils
# 2. Reading data
# 3. General plotting
# 4. Specialized plotting / Backward compatibility
# 5. Calculations


# 1. Utils
_units = {'X': 'm', 'PX': 'rad', 'Y': 'm', 'PY': 'rad', 'T': 'm', 'PT': '1',
          'S': 'm', 'E': 'eV', 'TURN': '1'}
_def_unit_exp = {key: (9 if key=='E' else 0) for key in _units}
_sidict = {18: 'E', 15: 'P', 12: 'T',
           9: 'G', 6: 'M', 3: 'k',
           0: '',
           -3: 'm', -6: '$\mu$', -9: 'n',
           -12: 'p', -15: 'f', -18: 'a'}
_sidict_rev = {value: key for key, value in iteritems(_sidict)}
    
# Nominal SPS apertures
nomap = {'tpstcirc': 0.03729401, # central orbit to inner blade edge (start tpst)
         'tpstblade': 0.0046, # blade thickness
         'tpstex': 0.04, # outer blade edge to outer aperture limit
         'zsupmid': 0.068, # centre of upstream wire
         'zsdomid': 0.0413, # centre of downstream wire
         'zsthick': 0.0002, # effective Zs thickness
         'zsex': 0.02} # centre wire to edge cathode

def errorcheck(errfolder):
    failed=[]
    messages=[]
    for errfile in os.listdir(errfolder):
        jobid = errfile.split(".")[0]
        if not os.path.exists(errfolder+"/../losses/"+jobid+".tfs"):
                failed += [jobid]
                if os.stat(errfolder+'/'+errfile).st_size > 0:
                    with open(errfolder+'/'+errfile, 'r') as f:
                        firstline = f.readline()
                    if firstline not in messages:
                        messages += [firstline]
                elif "Missing lossfile" not in messages:
                    messages += ["Missing lossfile"]
    return failed, messages

def fixlossfile(filename):
    linebuffer = ''
    for line in fileinput.input(filename, inplace=1):
        line = linebuffer+line
        breakdown = line.split('"')
        if len(breakdown)==2:
            linebuffer = line[:-1]
            line = ''
        elif len(breakdown)>=3:
            line = breakdown[0]+'"'+breakdown[1].split()[0]+'"\n'
        if len(line)>0:
            sys.stdout.write(line)
            linebuffer = ''
  
def getunits(data, forceunit=None):
    if forceunit is None:
        result = {}
    else:
        result = dict(forceunit)
    if not 'unit' in result:
        result['unit'] = _units[data.name]
    if not 'expon' in result:
        if 'prefix' in result:
            result['expon'] = _sidict_rev[result['prefix']]
        else:
            try:
                result['expon'] = int(round(np.log10(data.abs().mean())/3.0)*3.0)
            except ValueError:
                print('An error occured in getunits, mean was', data.abs().mean())
                result['expon'] = 0
            result['prefix'] = _sidict[result['expon']]
        try:
            result['expon'] -= _def_unit_exp[data.name]
        except:
            pass
    if not 'prefix' in result:
        try:
            result['prefix'] = _sidict[result['expon']-_def_unit_exp[data.name]]
        except:
            result['prefix'] = _sidict[result['expon']]
    if result['unit']=='1' and result['expon']!=0:
        result['unit'] = ''
        result['prefix'] = '$10^{'+str(result['expon'])+'}$'
    return result   

# Wilson score interval as described on wiki: https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval
def wilson(ns, ntot, conf=0.95):
    ntot = float(ntot)
    z = norm.ppf(1-(1-conf)/2)
    mid = (ns + z**2/2) / (ntot + z**2)
    dev = z/(ntot+z**2) * np.sqrt(ns*(ntot-ns)/ntot + z**2/4)
    return [mid-dev, mid+dev]


# 2. Reading data
def getsettings(folder):
    with open(folder+'/settings.info', 'r') as f:
        settings = {}
        first = True
        skip = True
        for line in f:
            if first:
                settings['version'] = line.strip().split()[-1][:-1]
                first = False
                continue
            if skip:
                if line.startswith('Settings used:'):
                    skip = False
                continue
            var, val = line.strip().split(" = ",1)
            settings[var.strip()] = val.strip()
    return settings

#TODO for efficiency should be more like
# frames = [ process_your_file(f) for f in files ]
# result = pd.concat(frames)
def getlosses(lossfolder, settings=None, lossloc=None, batches=None,
              filters=None, usecols=None):
    if usecols is not None:
        if 'NUMBER' not in usecols:
            usecols = usecols+['NUMBER']
    if lossloc is not None:
        if usecols is not None:
            if 'ELEMENT' not in usecols:
                usecols = usecols+['ELEMENT']
                dropel = True
            else:
                dropel = False
        else:
            dropel = False
    if settings is None:
        settings = getsettings(lossfolder+"/..")
    elif type(settings)==dict:
        pass
    else:
        settings = getsettings(settings)
        
    nppb = int(settings['nparperbatch'])
    if batches is None:
        nbatches = len(settings['slices'].split(','))*int(settings['nbatches'])
        batches = range(nbatches)
    df = pd.DataFrame()
    
    for batch in batches:
        lossfile = str(batch)+'.tfs'
        if os.path.isfile(lossfolder+'/'+lossfile):
            _, losstable = readtfs(lossfolder+'/'+lossfile, usecols=usecols)
            if lossloc is not None:
                losstable = losstable.loc[losstable['ELEMENT'] == lossloc]
                if dropel:
                    losstable.drop('ELEMENT', axis=1, inplace=True)
            if filters is not None:
                for (var, fun) in filters:
                    losstable = losstable.loc[losstable[var].apply(fun)]
            losstable.reset_index(inplace=True)
            losstable['NUMBER'] = losstable['NUMBER'] + int(lossfile[:-4])*nppb
            losstable.set_index('NUMBER', inplace=True)
            df = df.append(losstable)
        else:
            print("lossfile for batch "+str(batch)+" not found")
    df.sort_index(inplace=True)
    return df

#TODO for efficiency should be more like
# frames = [ process_your_file(f) for f in files ]
# result = pd.concat(frames)
#TODO make more efficient when getting tracks for multiple locations
def gettracks(trackfolder, settings=None, obsloc='obs0001', tpt=None,
              batches=None, filters=None, usecols=None, verbose=False):
    #Note mark_turns only works well if ffile=1 at the moment, otherwise need losstable
    """tpt=turns per track"""
    if usecols is not None:
        if 'NUMBER' not in usecols:
            usecols = usecols+['NUMBER']
        if 'TURN' not in usecols:
                usecols = usecols+['TURN']
    if settings is None:
        settings = getsettings(trackfolder+"/..")
    elif type(settings)==dict:
        pass
    else:
        settings = getsettings(settings)
    if not obsloc.startswith('obs0'):
        obsind = eval(settings['elements']).index(obsloc)
        obsloc = 'obs'+str(obsind+2).zfill(4)

    nppb = int(settings['nparperbatch'])
    if batches is None:
        nbatches = len(settings['slices'].split(','))*int(settings['nbatches'])
        batches = range(nbatches)

    def getbatch(batch):
        untarred = False
        if not os.path.isdir(trackfolder+'/'+str(batch)) and tarfile.is_tarfile(trackfolder+'/'+str(batch)+'.tar.gz'):
            if verbose:
                print('untarring tracks for batch '+str(batch))
            with tarfile.open(trackfolder+'/'+str(batch)+'.tar.gz') as tar:
                tar.extractall(path=trackfolder)
            untarred = True
        if os.path.isdir(trackfolder+'/'+str(batch)):
            if verbose:
                print('reading tracks for batch '+str(batch)+' from folder')
            dfb = pd.DataFrame()
            for part in np.arange(nppb)+1:
                trackfile = 'track.batch'+str(batch)+'.'+obsloc+'.p'+str(part).zfill(4)
                trackpath = trackfolder+'/'+str(batch)+'/'+trackfile
                if os.path.isfile(trackpath):
                    _, tracktable = readtfs(trackpath, usecols=usecols)
                    if filters is not None:
                        for (var, fun) in filters:
                            tracktable = tracktable.loc[tracktable[var].apply(fun)]
                    if tpt is not None:
                        tracktable = tracktable.tail(tpt)
                    tracktable.reset_index(inplace=True)
                    tracktable['obsnum'] = tracktable.index-len(tracktable)+1
                    tracktable['NUMBER'] = tracktable['NUMBER'] + batch*nppb
                    tracktable.set_index(['NUMBER','TURN'], inplace=True)
                    dfb = dfb.append(tracktable)
                else:
                    print('trackfile '+str(batch)+'/'+trackfile+' not found')
        else:
            print("tracks for batch "+str(batch)+" not found")
        if untarred:
            if verbose:
                print('deleting untarred files for batch '+str(batch))
            shutil.rmtree(trackfolder+'/'+str(batch))
        return dfb

    df = pd.concat([getbatch(b) for b in batches])
    df.sort_index(inplace=True)
    return df

def readtfs(filename, usecols=None, index_col=0, check_lossbug=True):
    header = {}
    nskip = 1
    closeit = False

    try:
        datafile = open(filename, 'r')
        closeit = True
    except TypeError:
        datafile = filename

    for line in datafile:
        nskip += 1
        if line.startswith('@'):
            entry = line.strip().split()
            header[entry[1]] = eval(' '.join(entry[3:]))
        elif line.startswith('*'):
            colnames = line.strip().split()[1:]
            break

    if closeit:
        datafile.close()

    table = pd.read_csv(filename, delim_whitespace = True,
                        skipinitialspace = True, skiprows = nskip,
                        names = colnames, usecols = usecols,
                        index_col = index_col)

    if check_lossbug:
        try:
            table['ELEMENT'] = table['ELEMENT'].apply(lambda x: str(x).split()[0])
        except KeyError:
            pass
        try:
            for location in table['ELEMENT'].unique():
                if not location.replace(".","").replace("_","").replace('$','').isalnum():
                    print("WARNING: some loss locations in "+filename+
                          " don't reduce to alphanumeric values. For example "+location)
                    break
                if location=="nan":
                    print("WARNING: some loss locations in "+filename+" are 'nan'.")
                    break
        except KeyError:
            pass

    return header, table


# 3. General plotting
class PlotGrid(object):
    """sns.JointGrid-like grid modified to our desires
       Not used yet...
    """
    def __init__(self, sizes=(15,5,1)):
        size = sum(sizes)/2.0
        fig = plt.figure(figsize=(size, size))
        gridsize = sum(sizes)
        gs = plt.GridSpec(gridsize, gridsize)

        js = sizes[1]+sizes[2]
        ms = sizes[2]

        ax_joint = fig.add_subplot(gs[js:, :-js])
        if sizes[1]>0:
         ax_marg_x = fig.add_subplot(gs[ms:js, :-js], sharex=ax_joint)
         ax_marg_y = fig.add_subplot(gs[js:, -js:-ms], sharey=ax_joint)
        else:
         ax_marg_x = None
         ax_marg_y = None
        if sizes[2]>0:
         ax_cbar = fig.add_subplot(gs[js:, -ms:])
        else:
         ax_cbar = None

        self.fig = fig
        self.ax_joint = ax_joint
        self.ax_marg_x = ax_marg_x
        self.ax_marg_y = ax_marg_y
        self.ax_cbar = ax_cbar

def lossmap(data, twiss, slim=None, merge=True, threshold=0.0001,
            extracted=None, save=None):
    # Kind of hacky with the twiss but ok
    mydata = lossstats(data, slim=slim, normalize=True, merge=merge)
    if extracted is not None:
        mydata[extracted[0]] -= extracted[1]
    mydata = mydata[mydata>=threshold]
    
    if not isinstance(twiss, pd.DataFrame):
        _, twiss = readtfs(twiss, usecols=['NAME', 'S'])
    
    mydata = mydata.to_frame(name='LOSS')
    mydata['S'] = twiss['S'].loc[mydata.index]
    mydata.sort_values('S', inplace=True)
    
    fig, ax = plt.subplots()
    xvals = lrange(1, len(mydata)+1)
    ax.bar(xvals, mydata['LOSS'].values, 0.05, color='k',bottom=threshold)
    ax.set_xlim(0,len(mydata)+1)
    ax.set_xticks(xvals)
    ax.set_xticklabels(mydata.index)
    fig.autofmt_xdate(bottom=0.2, rotation=30, ha='right')
    ax.set_yscale("log")
    plt.show()

    if save is not None:
        plt.savefig(save, bbox_inches='tight')
        plt.close()
        
    return fig, ax

def lossmapscan(lossdict, twiss, slim=None, merge=True, apmerge=True,
                clim=(5E-5,3E-2), extracted=None, save=None):
    if not isinstance(twiss, dict):
        onetwiss = True
        if not isinstance(twiss, pd.DataFrame):
            _, twiss = readtfs(twiss, usecols=['NAME', 'S'])
    else:
        onetwiss = False
    
    alldata = None
    for scanparam, data in iteritems(lossdict):
        mydata = lossstats(data, slim=slim, normalize=True, merge=merge, apmerge=apmerge)
        if extracted is not None:
            mydata[extracted[0]] -= extracted[1]
        mydata = mydata[mydata>=clim[0]]
        
        if onetwiss:
            mytwiss = twiss
        else:
            mytwiss = twiss[scanparam]
            if not isinstance(mytwiss, pd.DataFrame):
                _, mytwiss = readtfs(mytwiss, usecols=['NAME', 'S'])

        mydata = mydata.to_frame(name=scanparam)
        mydata['S'] = mytwiss['S'].loc[mydata.index]
        mydata.reset_index(inplace=True)
        mydata.set_index(['index','S'], inplace=True)
        
        if alldata is None:
            alldata = mydata.copy()
        else:
            alldata = alldata.merge(mydata, how='outer', left_index=True, right_index=True)
        
    alldata.reset_index(inplace=True)
    alldata.set_index('index', inplace=True)
    alldata.sort_values('S', inplace=True)
    
    labels = alldata.index.values
    scanvals = sorted(list(lossdict.keys()), key=float)
    
    alldata = alldata[scanvals]
    
    cm = plt.cm.get_cmap('RdYlGn_r')
    cm.set_under('White')
    fig, ax = plt.subplots()
    pcm = ax.imshow(alldata.transpose(), cmap=cm, norm=LogNorm(vmin=clim[0], vmax=clim[1]), interpolation='none', aspect='auto')
    
    ax.set_xlim(-0.5,len(labels)-0.5)
    ax.set_xticks(range(0, len(labels)))
    ax.set_xticklabels(labels)
    ax.set_ylim(-0.5,len(scanvals)-0.5)
    ax.set_yticks(range(0, len(scanvals)))
    ax.set_yticklabels(scanvals)
    fig.autofmt_xdate(bottom=0.2, rotation=30, ha='right')
    fig.colorbar(pcm, extend='min')

    if save is not None:
        plt.savefig(save, bbox_inches='tight')
        plt.close()
    
    return fig, ax

#TODO main plot types, with appropriate args
#TODO implement xkwargs ykwargs
#TODO colorbar? (colorbar units not ok...)
#TODO use custom grid thingy instead of the sns one?
#TODO font and font size
def plotter(data, kind='scatter', marginals=True, xax="X", yax="PX", cax=None,
            xlim=None, ylim=None, clim=[None,None], margxlim=None,
            margylim=None, xbin=None, ybin=None, log=False, mainkwargs=None,
            xkwargs=None, ykwargs=None, save=None):
    if mainkwargs is None:
        mainkwargs={}
    
    xunit = getunits(data[xax])
    yunit = getunits(data[yax])
    #if cax is not None:
    #    cunit = getunits(data[cax])
        
    if xbin is None:
        if xlim is None:
            xbin = (data[xax].max()-data[xax].min())/100.0
        else:
            xbin = (xlim[1]-xlim[0])/100.0
    if ybin is None:
        if ylim is None:
            ybin = (data[yax].max()-data[yax].min())/100.0
        else:
            ybin = (ylim[1]-ylim[0])/100.0
    if xlim is None:
        binsx = np.arange(data[xax].min()-xbin,
                          data[xax].max()+2*xbin, xbin)
    else:
        binsx = np.arange(xlim[0]-xbin, xlim[1]+xbin, xbin)
    if ylim is None:
        binsy = np.arange(data[yax].min()-ybin,
                          data[yax].max()+2*ybin, ybin)
    else:
        binsy = np.arange(ylim[0]-ybin, ylim[1]+ybin, ybin)
    xlim = (binsx[0], binsx[-1])
    ylim = (binsy[0], binsy[-1])
    
    # configure grid
    ratio = 4 if marginals else 100
    size = 6 if marginals else 4.5
    g = sns.JointGrid(xax, yax, data=data, size=size, ratio=ratio,
                      xlim=xlim, ylim=ylim)
    sns.despine(top=False, right=False)
    if marginals:
        plt.setp(g.ax_marg_x.get_yticklabels(), visible=True)
        plt.setp(g.ax_marg_y.get_xticklabels(), visible=True)
        plt.setp(g.ax_marg_x.yaxis.get_majorticklines(), visible=True)
        plt.setp(g.ax_marg_x.yaxis.get_minorticklines(), visible=True)
        plt.setp(g.ax_marg_y.xaxis.get_majorticklines(), visible=True)
        plt.setp(g.ax_marg_y.xaxis.get_minorticklines(), visible=True)
        plt.setp(g.ax_marg_x.get_yticklabels(), visible=True)
        plt.setp(g.ax_marg_y.get_xticklabels(), visible=True)
    else:
        g.ax_marg_x.set_axis_off()
        g.ax_marg_y.set_axis_off()
    
    weights = 100.0*np.ones_like(g.x)/len(g.x)
    
    # set main plot
    if kind=='scatter':
        mainkwargs.setdefault('s',4)
        if cax is not None:
            mainkwargs.setdefault('c', data[cax])
            mainkwargs.setdefault('vmin', clim[0])
            mainkwargs.setdefault('vmax', clim[1])
        g = g.plot_joint(plt.scatter, **mainkwargs)
    elif kind=='hexbin':
        if log:
            mainkwargs.setdefault('bins', 'log')
        g = g.plot_joint(plt.hexbin, extent=[xlim[0], xlim[1], ylim[0], ylim[1]],
                         gridsize=(len(binsx)-1,len(binsy)-1), **mainkwargs)
    elif kind=='hist2d':
        if log:
            mainkwargs.setdefault('norm', LogNorm())
        g.ax_joint.hist2d(g.x, g.y, bins=[binsx,binsy],
                          normed=True, weights=weights, **mainkwargs)
    else:
        print('Warning: "kind" not recognised')
        
    g.ax_joint.set_xlabel(xax+' ['+xunit['prefix']+xunit['unit']+']')
    g.ax_joint.set_ylabel(yax+' ['+yunit['prefix']+yunit['unit']+']')
    
    # set marginal plots
    if marginals:
        g.ax_marg_x.hist(g.x, bins=binsx, weights=weights, log=log)
        g.ax_marg_y.hist(g.y, bins=binsy, weights=weights, log=log,
                         orientation='horizontal')

        g.ax_marg_x.set_ylim(margxlim)
        g.ax_marg_y.set_xlim(margylim)

        g.ax_marg_x.set_ylabel('$\%$')
        g.ax_marg_y.set_xlabel('$\%$')
    
    # hack for hopefully better default ticks
    g.ax_joint.xaxis.set_minor_locator(ticker.AutoMinorLocator())
    g.ax_joint.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    if marginals:
        if log:
            g.ax_marg_x.yaxis.set_major_locator(ticker.LogLocator(base=10.0, subs=(1.0,), numticks=10))
            g.ax_marg_x.yaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs=(0.2,0.4,0.6,0.8), numticks=10))
            g.ax_marg_x.yaxis.set_minor_formatter(ticker.NullFormatter())

            g.ax_marg_y.xaxis.set_major_locator(ticker.LogLocator(base=100.0, subs=(1.0,), numticks=10))
            g.ax_marg_y.xaxis.set_minor_locator(ticker.LogLocator(base=100.0, subs=(0.2,0.4,0.6,0.8), numticks=10))
            g.ax_marg_y.xaxis.set_minor_formatter(ticker.NullFormatter())
        else:
            g.ax_marg_x.yaxis.set_minor_locator(ticker.AutoMinorLocator())
            g.ax_marg_y.xaxis.set_minor_locator(ticker.AutoMinorLocator())
            
    # fix units and formatting
    ticks = ticker.FuncFormatter(lambda x, pos: '${0:g}$'.format(x/10**xunit['expon']))
    g.ax_joint.xaxis.set_major_formatter(ticks)
    ticks = ticker.FuncFormatter(lambda x, pos: '${0:g}$'.format(x/10**yunit['expon']))
    g.ax_joint.yaxis.set_major_formatter(ticks)
    ticks = ticker.FuncFormatter(lambda x, pos: '${0:g}$'.format(x))
    g.ax_marg_x.yaxis.set_major_formatter(ticks)
    g.ax_marg_y.xaxis.set_major_formatter(ticks)

    # colorbar
    #cb = plt.colorbar(ax=[g.ax_joint, g.ax_marg_x, g.ax_marg_y], use_gridspec=False)
    #cb.set_label('$'+cax+' ['+cunit['prefix']+cunit['unit']+']$')
    #ticks = ticker.FuncFormatter(lambda x, pos: '${0:g}$'.format(x/10**cunit['expon']))
    #cb.ax.yaxis.set_major_formatter(ticks)

    # save and exit
    if save is not None:
        plt.savefig(save)
        plt.close()
    
    return g
    
    
# 4. Specialized plotting / Backward compatibility
def losshistscatter(lossfolder, lossloc="AP.UP.ZS21633", xax="X", yax="PX",
                    cax="TURN", xlim=None, ylim=None, clim=[None,None],
                    monochrom=False, xbin=None, ybin=None, log=False,
                    extra=None, save=None,
                    datalim=[[None,None],[None,None],[None,None]]):
    # datalim and extra removed
    usecols = ['NUMBER', xax, yax] if (cax is None) else ['NUMBER', xax, yax, cax]
    if xax=='S': lossloc=None
    data = getlosses(lossfolder, lossloc=lossloc, usecols=usecols)
    if monochrom:
        cax=None
    plotter(data, kind='scatter', xax=xax, yax=yax, cax=cax, xlim=xlim,
            ylim=ylim, clim=clim, xbin=xbin, ybin=ybin, log=False, save=save)
    
def trackplot(trackfolder, obsloc="obs0001", xax="X", yax="PX", cax="TURN",
              xlim=None, ylim=None, clim=[None,None], tpt=3, batches=None, save=None):
    usecols = ['NUMBER', 'TURN', xax, yax] if (cax is None) else ['NUMBER', 'TURN', xax, yax, cax]
    data = gettracks(trackfolder, obsloc=obsloc, tpt=tpt, batches=batches, usecols=usecols)
    plotter(data, kind='scatter', xax=xax, yax=yax, cax=cax, xlim=xlim,
            ylim=ylim, clim=clim, log=False, save=save)

def fullplot(folder, lossloc="AP.UP.ZS21633", obsloc="obs0002", xax="X",
             yax="PX", cax="TURN",  xlim=None, ylim=None, clim=[None,None],
             tpt=3, batches=None, save=None):
    usecols = ['NUMBER', 'TURN', xax, yax] if (cax is None) else ['NUMBER', 'TURN', xax, yax, cax]
    trackdata = gettracks(folder+"/tracks", obsloc=obsloc, tpt=tpt, batches=batches, usecols=usecols).reset_index()
    lossdata = getlosses(folder+"/losses", lossloc=lossloc, batches=batches, usecols=usecols).reset_index()
    data = pd.concat([trackdata, lossdata], ignore_index=True)
    plotter(data, kind='scatter', xax=xax, yax=yax, cax=cax, xlim=xlim,
            ylim=ylim, clim=clim, marginals=False, log=False, save=save)


# 5. Calculations
def beamstats(data, location=None, plane='X', p0=400, disp=(None,None), save=None, silent=True):
    xax = plane
    yax = 'P'+plane
    
    if save is not None:
        silent = False
    
    if not isinstance(data, pd.DataFrame):
        usecols = ['NUMBER', xax, yax]
        data = getlosses(data, lossloc=location, usecols=usecols)
    elif location is not None:
        data = data[data['ELEMENT']==location]

    result = {}
    
    result['location'] = location
    result['plane'] = plane
    result['npart'] = len(data)
    result['p0'] = p0
    
    betagamma = p0/0.938272081
    
    result[xax+'_min'] = data[xax].min()
    result[xax+'_max'] = data[xax].max()
    result[xax+'_avg'] = data[xax].mean()
    result[xax+'_std'] = data[xax].std()
    
    result[yax+'_min'] = data[yax].min()
    result[yax+'_max'] = data[yax].max()
    result[yax+'_avg'] = data[yax].mean()
    result[yax+'_std'] = data[yax].std()

    if disp[0] is None:
        fit = np.polyfit(data['PT'],data[xax],1)
        result['D'+xax] = fit[0]
        result[xax+'0'] = fit[1]
    else:
        result['D'+xax] = disp[0]
        result[xax+'0'] = result[xax+'_avg']
        
    if disp[1] is None:
        fit = np.polyfit(data['PT'],data[yax],1)
        result['D'+yax] = fit[0]
        result[yax+'0'] = fit[1]
    else:
        result['D'+yax] = disp[1]
        result[yax+'0'] = result[yax+'_avg']
    
    covmat = np.cov(data[xax]-result['D'+xax]*data['PT'], data[yax]-result['D'+yax]*data['PT'])
    
    result['emittance']  = np.sqrt(np.linalg.det(covmat))
    result['emittance_norm']  = result['emittance']*betagamma
    result['alpha'] = -1*covmat[0,1]/result['emittance']
    result['beta'] = covmat[0,0]/result['emittance']
    result['gamma'] = covmat[1,1]/result['emittance']
    
    if not silent:
        message = ('min,max '+xax+': '+str(result[xax+'_min'])+', '+str(result[xax+'_max'])+'\n'+
                   'avg,std '+xax+': '+str(result[xax+'_avg'])+', '+str(result[xax+'_std'])+'\n'+
                   'min,max '+yax+': '+str(result[yax+'_min'])+', '+str(result[yax+'_max'])+'\n'+
                   'avg,std '+yax+': '+str(result[yax+'_avg'])+', '+str(result[yax+'_std'])+'\n'+
                   'estimated orbit: '+str((result[xax+'0'], result[yax+'0']))+'\n'+
                   'estimated dispersion: '+str((result['D'+xax], result['D'+yax]))+'\n'+
                   'statistical emittance: '+str(result['emittance'])+'\n'+
                   'normalized statistical emittance at '+str(p0)+' GeV/c: '+str(result['emittance_norm'])+'\n'+
                   'stat alpha,beta,gamma: '+str([result['alpha'], result['beta'], result['gamma']])+'\n'+
                   '(Calculated with '+str(result['npart'])+' particles.)')

        if save is None:
            print(message)
        else:
            with open(save, 'w') as out:
                out.write(message)
            
    return result


def extr_tagger(data=None, pycoll=False, aperturex=[0,1], aperturey=[0,1],
                zs_len=1, zs_an=0.1, aperturex2=[0,1], aperturey2=None,
                losslocs=[]):
    if aperturey2 is None:
        aperturey2=aperturey
        
    if pycoll:
        alltags = ['extracted', 'tpst', 'tce', 'zs_wire', 'zs_app']+losslocs+['other']
        def tagger(particle):
            tag = 'other'
            if particle['ELEMENT']=="AP.UP.TPST21760":
                if (aperturex[0]<particle["X"]<aperturex[1]
                    and aperturey[0]<particle["Y"]<aperturey[1]):
                    tag = 'extracted'
                else:
                    tag = 'tpst'
            elif "TPST" in particle['ELEMENT']:
                tag = 'tpst'
            elif "TCE" in particle['ELEMENT']:
                tag = 'tce'
            elif particle['ELEMENT']==("ZS.21655"):
                tag = 'zs_wire'
            elif "ZS" in particle['ELEMENT']:
                tag = 'zs_app'
            else:
                for loc in losslocs:
                    if particle['ELEMENT']==loc:
                        tag = loc
            return tag
    else:
        alltags = ['extracted', 'zs_vert', 'zs_wire', 'zs_cath']+losslocs+['other']
        def tagger(particle):
            tag = 'other'
            if particle['ELEMENT']=="AP.UP.ZS21633":
                quadraticA = zs_an/zs_len/2
                quadraticB = (particle["PX"]
                              - (aperturex2[0]-aperturex[0])/zs_len)
                quadraticC = particle["X"]-aperturex[0]
                quadraticD = quadraticB**2 - 4*quadraticA*quadraticC
                if quadraticD >= 0:
                    hitdist = (-1*quadraticB - quadraticD**0.5) / (2*quadraticA)
                else:
                    hitdist = 0
                    
                if not (aperturey[0]<particle["Y"]<aperturey[1]):
                    tag = 'zs_vert'
                elif particle["X"]<aperturex[0]:
                    tag = 'zs_wire'
                elif particle["X"]>aperturex[1]:
                    tag = 'zs_cath'
                elif not (aperturey2[0]<(particle["Y"]+zs_len*particle["PY"])<aperturey2[1]):
                    tag = 'zs_vert'
                elif 0<hitdist<zs_len:
                    tag = 'zs_wire'
                elif particle["X"]+(particle["PX"]+zs_an/2)*zs_len>aperturex2[1]:
                    tag = 'zs_cath'
                else:
                    tag = 'extracted'
            elif "ZS" in particle['ELEMENT']:
                tag = 'zs_wire'
            else:
                for loc in losslocs:
                    if particle['ELEMENT']==loc:
                        tag = loc
            return tag
        
    if data is not None:
        data['tag'] = data.apply(tagger, axis=1)
        
    return alltags, tagger


def efficiency(data, alltags=['extracted', 'other'], tagger=None, conf=0.95,
               save=None, silent=True):
    if not isinstance(data, pd.DataFrame):
        usecols = ['NUMBER', 'ELEMENT', 'X', 'PX', 'Y', 'PY']
        data = getlosses(data, lossloc=None, usecols=usecols)
    
    if save is not None:
        silent = False
        
    if tagger is not None:
        tags = data.apply(tagger, axis=1)
    else:
        tags = data['tag']
    ntot = len(tags)
    result = tags.value_counts()
    
    for tag in alltags:
        if tag not in result.index.values:
            result[tag] = 0
            
    ci = result.apply(lambda x: wilson(x,ntot, conf=conf))
    result = result/ntot

    if not silent:
        message = str(100*result['extracted'])+' % extracted\t'+str([100*x for x in ci['extracted']])+'\n'
        for tag, frac in result.iteritems():
            if tag not in ['extracted', 'other']:
                message += str(100*frac)+' % lost on '+tag+'\t'+str([100*x for x in ci[tag]])+'\n'
        message += str(100*result['other'])+' % lost elsewhere\t'+str([100*x for x in ci['other']])+'\n'
        message += '(total: '+str(len(data))+' particles)'

        if save is None:
            print(message)
        else:
            with open(save, 'w') as out:
                out.write(message)
                
    return result, ci


def getactionangle(inits, twiss, plane='X'):
    # assumes action is given by twiss ellipses
    if not isinstance(inits, pd.DataFrame):
        usecols = ['NUMBER', plane, 'P'+plane, 'PT']
        _, inits = readtfs(inits, usecols=usecols)
    if not isinstance(twiss, pd.DataFrame):
        usecols = ['NAME', plane, 'P'+plane, 'D'+plane, 'DP'+plane,
                   'ALF'+plane, 'BET'+plane]
        _, twiss = readtfs(twiss, usecols=usecols)
    mytwiss = twiss.iloc[0]
    
    realpos = inits[plane] - mytwiss[plane] - mytwiss['D'+plane]*inits['PT']
    realang = inits['P'+plane] - mytwiss['P'+plane] - mytwiss['DP'+plane]*inits['PT']
    
    action = (mytwiss['BET'+plane]*realang**2 + 2*mytwiss['ALF'+plane]*realpos*realang
              + (1+mytwiss['ALF'+plane]**2)/mytwiss['BET'+plane]*realpos**2)/2.0
    
    cosang = realpos / np.sqrt(2*action*mytwiss['BET'+plane])
    sinang = -realang / np.sqrt(2*action/mytwiss['BET'+plane]) - mytwiss['ALF'+plane]*cosang
    
    angle = np.arctan2(sinang,cosang)
    
    result = pd.DataFrame({'J'+plane: action, 'PHI'+plane: angle})
    
    return result

def lossstats(data, slim=None, merge=True, apmerge=False, normalize=False,
              save=None, silent=True):
    # Note: if s limits are given and normalization is on, losses will be normalized to the loss within that region.
    # Note: if some slices are inside the s limits and some out, only losses on the slices inside will be taken into account, even is merge=True
    if not isinstance(data, pd.DataFrame):
        usecols = ['NUMBER', 'S', 'ELEMENT']
        kwargs = {}
        if slim is not None:
            kwargs['filters'] = [('S', lambda s: s>=slim[0] and s<=slim[1])]
        data = getlosses(data, lossloc=None, usecols=usecols, **kwargs)
    elif slim is not None:
        data = data[data['S'].between(*slim)]
        
    if 'tag' in data:
        data = data[data['tag']!='extracted']
        
    losslocs = data['ELEMENT']
    if merge:
        losslocs = losslocs.apply(lambda x: x.split('..')[0])
    if apmerge:
        f = lambda x: x[:-5]+'.'+x[-5:] # Dodgy but works...
        g = lambda y: y if not y.startswith('AP') else f(y.split('.')[-1])
        losslocs = losslocs.apply(g)

    lossstats = losslocs.value_counts()
        
    lossstats.name = 'LOSSSTATS'
    if normalize:
        lossstats /= len(data)

    if save is not None:
        lossstats.to_csv(save, sep=' ')
    elif not silent:
        print(lossstats.to_csv(sep=' '))

    return lossstats
