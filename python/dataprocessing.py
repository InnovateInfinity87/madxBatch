import os
import sys
import numpy as np
import pandas as pd
import re
import operator
import fileinput
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.lines import Line2D
from matplotlib.ticker import NullFormatter
from matplotlib.colors import LogNorm

_units = {'X': 'm', 'PX': 'rad', 'Y': 'm', 'PY': 'rad', 'T': 'm', 'PT': '1', 'S': 'm', 'E': 'eV', 'TURN': '1'}
_def_unit_exp = {'X': 0, 'PX': 0, 'Y': 0, 'PY': 0, 'T': 0, 'PT': 0, 'S': 0, 'E': 9, 'TURN':0}

def readtfs(filename, usecols=None, index_col=0):
    header = {}
    nskip=1

    with open(filename, 'r') as datafile:
        for line in datafile:
            nskip += 1
            if line.startswith('@'):
                entry = line.strip().split()
                header[entry[1]] = ' '.join(entry[3:]).replace('"','')
            elif line.startswith('*'):
                colnames = line.strip().split()[1:]
                break

    if usecols is not None:
        colnames = [colnames[i] for i in usecols]

    table = pd.read_csv(filename, delim_whitespace = True,
                        skipinitialspace = True, skiprows = nskip,
                        names = colnames, usecols = usecols,
                        index_col = index_col)

    # Take care of lossloc bug in madxColl...
    try:
        table['ELEMENT'] = table['ELEMENT'].apply(lambda x: str(x).split()[0])
    except KeyError:
        pass

    # Check for the extra buggy things in lossfiles
    try:
        for location in table['ELEMENT']:
            if not location.replace(".","").replace("_","").isalnum():
                print "WARNING: some loss locations in "+filename+" don't reduce to alphanumeric values. For example "+location
                break
            if location=="nan":
                print "WARNING: some loss locations in "+filename+" are 'nan'."
                break
    except KeyError:
        pass

    return header, table

def fixlossfile(filename):
    question = 'This function is dangerous, only use on lossfile with accidental additional "-symbols or newlines.\n Do you want to proceed?'
    reply = str(raw_input(question+' (y/n): ')).lower().strip()
    if reply[0] == 'y':
        linebuffer = ''
        for line in fileinput.input(filename, inplace=1):
            line = linebuffer+line
            breakdown = line.split('"')
            if len(breakdown)==2:
                linebuffer = line[:-1]
                line = ''
            elif len(breakdown)>3:
                line = breakdown[0]+'"'+breakdown[1]+'"\n'
            if len(line)>0:
                sys.stdout.write(line)
                linebuffer = ''
    else:
        print 'Aborting.'


def readsingletrack(filename):
    header, table = readtfs(filename, usecols=[1,2,3,4,5,6,7,8,9])
    return header, table


# TODO: change scale
def lossplot(lossfolder, lossloc="AP.UP.ZS21633", xax="X", yax="PX", cax="TURN", xlim=None, ylim=None, clim=[None,None], save=None):
    if xax=='S':
        lossloc = None
    xdata = []
    ydata = []
    cdata = []

    empty = 0

    for lossfile in os.listdir(lossfolder):
        if os.stat(lossfolder+'/'+lossfile).st_size > 0:
            _, losstable = readtfs(lossfolder+'/'+lossfile)
            if lossloc is not None:
                losstable = losstable.loc[losstable['ELEMENT'] == lossloc]
            xdata += losstable[xax].tolist()
            ydata += losstable[yax].tolist()
            cdata += losstable[cax].tolist()
        else:
            empty += 1

    if empty>0:
        print "warning: "+str(empty)+" empty loss files found!"

    xunit = _units[xax]
    yunit = _units[yax]
    
    cm = plt.cm.get_cmap('viridis')
    fig, ax = plt.subplots()
    plt.autoscale(enable=True, axis='both', tight=True)
    plot = ax.scatter(xdata, ydata, c=cdata, cmap=cm, vmin=clim[0], vmax=clim[1], edgecolor='')
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)
    plt.colorbar(plot)

    ax.set_xlabel(xax+' ['+xunit+']')
    ax.set_ylabel(yax+' ['+yunit+']')

    if save is None:
        plt.show()
    else:
        plt.savefig(save)
        plt.close()


# TODO: normalize y, change scale x, label y
def losschart(lossfolder, lossloc="AP.UP.ZS21633", xax="X", binwidth=5E-4, startx=None, endx=None, save=None):
    if xax=='S':
        lossloc = None
    data = []

    for lossfile in os.listdir(lossfolder):
        _, losstable = readtfs(lossfolder+'/'+lossfile)
        if lossloc is not None:
            losstable = losstable.loc[losstable['ELEMENT'] == lossloc]
        data += losstable[xax].tolist()

    xunit = _units[xax]

    if startx is None:
        startx=min(data)
    if endx is None:
        endx=max(data)+binwidth
    bins = np.arange(startx, endx + binwidth, binwidth)
    
    fig, ax = plt.subplots()
    plot = ax.hist(data, bins=bins, facecolor='c', alpha=0.75)

    ax.set_xlabel(xax+' ['+xunit+']')

    if save is None:
        plt.show()
    else:
        plt.savefig(save)
        plt.close()


# TODO: normalize y, change scale x, label y
def lossstats(lossfolder, region=None, printstats=True):
    data = {}

    for lossfile in os.listdir(lossfolder):
        if os.stat(lossfolder+'/'+lossfile).st_size > 0:
            _, losstable = readtfs(lossfolder+'/'+lossfile)
            if region is not None:
                losstable = losstable.loc[losstable['S']>=region[0] and losstable['S']<=region[1]]
            for location in losstable['ELEMENT'].tolist():
                try:
                    data[location] += 1
                except KeyError:
                    data[location] = 1

    if printstats:
        print "Start loss stats:"
        for location in sorted(data, key=data.get, reverse=True):
            print location, data[location]
        print "End loss stats."

    return data


def lossmap(lossfolder, region=None, save=None, threshold=0.0001, extracted=None):
    rawdata = lossstats(lossfolder)
    total = float(sum(rawdata.values()))
    for x in rawdata: rawdata[x] /= total
    if extracted is not None:
        rawdata[extracted[0]] -= extracted[1]
    data = {}
    positions = {}
    _, twissfile = readtfs(lossfolder+"/../thin_twiss.tfs")

    for location in rawdata:
            # Reduce sliced elements to one
            if re.search(r'\.\.[0-9]*$', location) is not None:
                newloc = "..".join(location.split("..")[:-1])
                try:
                    data[newloc] += rawdata[location]
                except KeyError:
                    data[newloc] = rawdata[location]
            else:
                try:
                    data[location] += rawdata[location]
                except KeyError:
                    data[location] = rawdata[location]

    for location in data:
        # Keep only significant loss locations
        if data[location] >= threshold:
            positions[location] = twissfile.loc[location,"S"]

    labels = [x[0] for x in sorted(positions.items(), key=operator.itemgetter(1))]
    lossvals = [data[location] for location in labels]
    labels =[x[:-2] if x.endswith("_S") else x for x in labels]
    xvals = range(1, len(labels)+1)

    fig, ax = plt.subplots()
    ax.bar(xvals, lossvals, 0.05, color='k')
    ax.set_xlim(0,len(labels)+1)
    ax.set_xticks(xvals)
    ax.set_xticklabels(labels)
    fig.autofmt_xdate(bottom=0.2, rotation=30, ha='right')
    ax.set_yscale("log")

    if save is None:
        plt.show()
    else:
        plt.savefig(save)
        plt.close()


def lossmapscan(lossfolders, region=None, save=None, threshold=0.0001, extracted=None, clim=[5E-5,3E-2]):
    rawdata = {}
    data = {}
    positions = {}
    notwiss = True
    for scanparam, lossfolder in lossfolders.iteritems():
        if notwiss:
            _, twissfile = readtfs(lossfolder+"/../thin_twiss.tfs")
            masterkey = scanparam
            notwiss = False

        rawdata[scanparam] = lossstats(lossfolder, printstats=False)
        total = float(sum(rawdata[scanparam].values()))
        for x in rawdata[scanparam]: rawdata[scanparam][x] /= total
        if extracted is not None:
            rawdata[scanparam][extracted[0]] -= extracted[1][scanparam]
        data[scanparam] = {}

        for location in rawdata[scanparam]:
                # Reduce sliced elements to one
                if re.search(r'\.\.[0-9]*$', location) is not None:
                    newloc = "..".join(location.split("..")[:-1])
                    try:
                        data[scanparam][newloc] += rawdata[scanparam][location]
                    except KeyError:
                        data[scanparam][newloc] = rawdata[scanparam][location]
                else:
                    try:
                        data[scanparam][location] += rawdata[scanparam][location]
                    except KeyError:
                        data[scanparam][location] = rawdata[scanparam][location]

        for location in data[scanparam]:
            # Keep only significant loss locations
            if data[scanparam][location] >= threshold:
                positions[location] = twissfile.loc[location,"S"]

    labels = [x[0] for x in sorted(positions.items(), key=operator.itemgetter(1))]
    scanvals = sorted(list(data.keys()))
    plotvals = np.zeros((len(scanvals), len(labels)))

    for i,scanval in enumerate(scanvals):
        for j,label in enumerate(labels):
            try:
                plotvals[i,j] = data[scanval][label]
            except KeyError:
                plotvals[i,j] = clim[0]/2.0

    labels =[x[:-2] if x.endswith("_S") else x for x in labels]

    cm = plt.cm.get_cmap('RdYlGn_r')
    cm.set_under('White')
    fig, ax = plt.subplots()
    pcm = ax.imshow(plotvals, cmap=cm, norm=LogNorm(vmin=clim[0], vmax=clim[1]), interpolation='none', aspect='auto')

    ax.set_xlim(-0.5,len(labels)-0.5)
    ax.set_xticks(range(0, len(labels)))
    ax.set_xticklabels(labels)
    ax.set_ylim(-0.5,len(scanvals)-0.5)
    ax.set_yticks(range(0, len(scanvals)))
    ax.set_yticklabels(scanvals)

    fig.autofmt_xdate(bottom=0.2, rotation=30, ha='right')
    fig.colorbar(pcm, extend='min')

    if save is None:
        plt.show()
    else:
        plt.savefig(save)
        plt.close()


# TODO: include Francesco's loss calculation?
def beamstats(lossfolder, lossloc="AP.UP.ZS21633", plane="X", save=None):
    if plane=='X':
        xax='X'
        pax='PX'
    elif plane=='Y':
        xax='Y'
        pax='PY'
    else:
        print "plane should be 'X' or'Y'"
        return

    xdata = []
    pdata = []
    totloss = 0

    for lossfile in os.listdir(lossfolder):
        _, losstable = readtfs(lossfolder+'/'+lossfile)
        totloss += len(losstable.index)
        if lossloc is not None:
            losstable = losstable.loc[losstable['ELEMENT'] == lossloc]
        xdata += losstable[xax].tolist()
        pdata += losstable[pax].tolist()

    minx = str(min(xdata))
    avgx = str(np.mean(xdata))
    maxx = str(max(xdata))
    minp = str(min(pdata))
    avgp = str(np.mean(pdata))
    maxp = str(max(pdata))
    stdx = str(np.std(xdata))
    stdp = str(np.std(pdata))
    covmat = np.cov(xdata,pdata)
    statemit = np.sqrt(np.linalg.det(covmat))

    if lossloc is None:
        message = "Globally:\n"
    else:
        message = "At "+lossloc+":\n"

    message += ("min,max "+xax+": "+minx+", "+maxx+"\n"+
                "avg,std "+xax+": "+avgx+", "+stdx+"\n"+
                "min,max "+pax+": "+minp+", "+maxp+"\n"+
                "avg,std "+pax+": "+avgp+", "+stdp+"\n"+
                "statistical emittance: "+str(statemit)+"\n"+
                "stat alpha,beta,gamma: "+str([-1*covmat[0,1]/statemit, covmat[0,0]/statemit, covmat[1,1]/statemit])+"\n"+
                "(# of particles local/global: "+str(len(xdata))+"/"+str(totloss)+")")

    if save is None:
        print message
    else:
        with open(save, 'w') as out:
            out.write(message)


def efficiency(lossfolder, pycoll=False, aperturex=[0,1], aperturey=[0,1],
               zs_len=1, zs_an=0.1, aperturex2=[0,1], aperturey2=None, report=True, errorbin=None, losslocs=[]):
    xax='X'
    yax='Y'
    losses={}
    binnum=0

    if aperturey2 is None:
        aperturey2=aperturey

    empty = 0

    if pycoll:
        losses_t = ['extracted', 'zs_wire', 'zs_app', 'tce', 'tpst', 'other']+losslocs
        losses[0]={}
        for key in losses_t:
            losses[0][key] = 0

        for lossfile in os.listdir(lossfolder):
            if os.stat(lossfolder+'/'+lossfile).st_size > 0:
                _, losstable = readtfs(lossfolder+'/'+lossfile)
                for pid, particle in losstable.iterrows():
                    if errorbin is not None:
                        binnum = int(pid//errorbin)
                        if binnum not in losses:
                            losses[binnum] = {}
                            for key in losses_t:
                                losses[binnum][key] = 0
                    if particle['ELEMENT']=="AP.UP.TPST21760":
                        if (particle["X"]>aperturex[0] and particle["X"]<aperturex[1]
                            and particle["Y"]>aperturey[0] and particle["Y"]<aperturey[1]):
                            losses[binnum]['extracted'] += 1
                        else:
                            losses[binnum]['tpst'] += 1
                    elif "TPST" in particle['ELEMENT']:
                        losses[binnum]['tpst'] += 1
                    elif "TCE" in particle['ELEMENT']:
                        losses[binnum]['tce'] += 1
                    elif particle['ELEMENT']==("ZS.21655"):
                        losses[binnum]['zs_wire'] += 1
                    elif "ZS" in particle['ELEMENT']:
                        losses[binnum]['zs_app'] += 1
                    else:
                        other = True
                        for loc in losslocs:
                            if particle['ELEMENT']==loc:
                                losses[binnum][loc] += 1
                                other = False
                        if other:
                            losses[binnum]['other'] += 1
            else:
                empty += 1

    else:
        losses_t = ['extracted', 'zs_wire', 'zs_cath', 'zs_vert', 'other']+losslocs
        losses[0]={}
        for key in losses_t:
            losses[0][key] = 0

        for lossfile in os.listdir(lossfolder):
            if os.stat(lossfolder+'/'+lossfile).st_size > 0:
                _, losstable = readtfs(lossfolder+'/'+lossfile)
                for pid, particle in losstable.iterrows():
                    if errorbin is not None:
                        binnum = int(pid//errorbin)
                        if binnum not in losses:
                            losses[binnum] = {}
                            for key in losses_t:
                                losses[0][key] = 0
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

                        if particle["Y"]<aperturey[0] or particle["Y"]>aperturey[1]:
                            losses[binnum]['zs_vert'] += 1
                        elif particle["X"]<aperturex[0]:
                            losses[binnum]['zs_wire'] += 1
                        elif particle["X"]>aperturex[1]:
                            losses[binnum]['zs_cath'] +=1
                        elif (particle["Y"]+zs_len*particle["PY"]<aperturey2[0]
                              or particle["Y"]+zs_len*particle["PY"]>aperturey2[1]):
                            losses[binnum]['zs_vert'] += 1
                        elif hitdist>0 and hitdist<zs_len:
                            losses[binnum]['zs_wire'] += 1
                        elif particle["X"]+(particle["PX"]+zs_an/2)*zs_len>aperturex2[1]:
                            losses[binnum]['zs_cath'] += 1
                        else:
                            losses[binnum]['extracted'] += 1
                    elif "ZS" in particle['ELEMENT']:
                        losses[binnum]['zs_wire'] += 1
                    else:
                        other = True
                        for loc in losslocs:
                            if particle['ELEMENT']==loc:
                                losses[binnum][loc] += 1
                                other = False
                        if other:
                            losses[binnum]['other'] += 1
            else:
                empty += 1

    if empty > 0:
        print str(empty)+" empty loss files found"

    # Add total, convert to %
    for binnum in losses:
        losses[binnum]['total'] = sum(losses[binnum].values())
        for key in losses_t:
            try:
                losses[binnum][key] = 100.0*losses[binnum][key]/losses[binnum]['total']
            except ZeroDivisionError:
                losses[binnum][key] = 0

    # Determine means and errorbars
    if errorbin is None:
        reslosses = losses[0]
        reslosses['truetotal'] = reslosses['total']
    else:
        reslosses = {}
        for key in losses_t+['total']:
            temp = [losses[i][key] for i in losses]
            reslosses[key] = np.mean(temp)
            reslosses['std_'+key] = np.std(temp)
        reslosses['truetotal'] = 0
        for binnum in losses:
            reslosses['truetotal'] += losses[binnum]['total']

    if report:
        if pycoll:
            if errorbin is None:
                print (str(reslosses['extracted'])+" % extracted\n"
                       +str(reslosses['zs_wire'])+" % lost on the ZS wires\n"
                       +str(reslosses['zs_app'])+" % lost on the ZS aperture\n"
                       +str(reslosses['tce'])+" % lost on the TCE\n"
                       +str(reslosses['tpst'])+" % lost on the TPST\n"
                       +str(reslosses['other'])+" % lost elsewhere.\n"
                       +"(Total: "+str(reslosses['truetotal'])+" particles.)")
            else:
                print (str(reslosses['extracted'])+'+-'+str(reslosses['std_extracted'])+" % extracted\n"
                       +str(reslosses['zs_wire'])+'+-'+str(reslosses['std_zs_wire'])+" % lost on the ZS wires\n"
                       +str(reslosses['zs_app'])+'+-'+str(reslosses['std_zs_app'])+" % lost on the ZS aperture\n"
                       +str(reslosses['tce'])+'+-'+str(reslosses['std_tce'])+" % lost on the TCE\n"
                       +str(reslosses['tpst'])+'+-'+str(reslosses['std_tpst'])+" % lost on the TPST\n"
                       +str(reslosses['other'])+'+-'+str(reslosses['std_other'])+" % lost elsewhere.\n"
                       +"(Total: "+str(reslosses['truetotal'])+" particles.)")

        else:
            if errorbin is None:
                print (str(reslosses['extracted'])+" % extracted\n"
                       +str(reslosses['zs_vert'])+" % lost on the ZS vertical aperture\n"
                       +str(reslosses['zs_wire'])+" % lost on the ZS wires\n"
                       +str(reslosses['zs_cath'])+" % lost on the ZS cathode\n"
                       +str(reslosses['other'])+" % lost elsewhere.\n"
                       +"(Total: "+str(reslosses['truetotal'])+" particles.)")
            else:
                print (str(reslosses['extracted'])+'+-'+str(reslosses['std_extracted'])+" % extracted\n"
                       +str(reslosses['zs_vert'])+'+-'+str(reslosses['std_zs_vert'])+" % lost on the ZS vertical aperture\n"
                       +str(reslosses['zs_wire'])+'+-'+str(reslosses['std_zs_wire'])+" % lost on the ZS wires\n"
                       +str(reslosses['zs_cath'])+'+-'+str(reslosses['std_zs_cath'])+" % lost on the ZS cathode\n"
                       +str(reslosses['other'])+'+-'+str(reslosses['std_other'])+" % lost elsewhere.\n"
                       +"(Total: "+str(reslosses['truetotal'])+" particles.)")


    return reslosses


def wireangle(lossfolder, lossloc="AP.UP.ZS21633", wiremax=0.69, save=None):
    xax='X'
    pax='PX'

    xdata = []
    pdata = []

    for lossfile in os.listdir(lossfolder):
        if os.stat(lossfolder+'/'+lossfile).st_size > 0:
            _, losstable = readtfs(lossfolder+'/'+lossfile)
            losstable = losstable.loc[losstable['ELEMENT'] == lossloc]
            xdata += losstable[xax].tolist()
            pdata += losstable[pax].tolist()

    wireang = [pdata[i] for i in range(len(pdata)) if xdata[i]<wiremax]

    minp = str(min(wireang))
    avgp = str(np.mean(wireang))
    maxp = str(max(wireang))
    stdp = str(np.std(wireang))

    message = ("Angular spread at "+lossloc+
               ", for lost particles with x<"+str(wiremax)+":\n"+
               "min,max "+pax+": "+minp+", "+maxp+"\n"+
               "avg,std "+pax+": "+avgp+", "+stdp+"\n"+
               "("+str(len(wireang))+" particles)\n")

    if save is None:
        print message
    else:
        with open(save, 'w') as out:
            out.write(message)

def emittance(lossfolder, lossloc="AP.UP.ZS21633", ap=[0.06815, 0.08815], save=None):
    xax='X'
    pax='PX'

    xdata = []
    pdata = []

    for lossfile in os.listdir(lossfolder):
        if os.stat(lossfolder+'/'+lossfile).st_size > 0:
            _, losstable = readtfs(lossfolder+'/'+lossfile)
            losstable = losstable.loc[losstable['ELEMENT'] == lossloc]
            xdata += losstable[xax].tolist()
            pdata += losstable[pax].tolist()

    exbeam = [i for i in range(len(xdata)) if (xdata[i]>ap[0] and xdata[i]<ap[1])]
    xdata = [xdata[i] for i in exbeam]
    pdata = [pdata[i] for i in exbeam]

    statemit = str(np.sqrt(np.linalg.det(np.cov(xdata,pdata))))

    message = ("For particles within the aperture at "+lossloc+":\n"+
               "the statistical emittance is "+statemit+"\n")

    if save is None:
        print message
    else:
        with open(save, 'w') as out:
            out.write(message)

def errorcheck(errfolder, checkerr=True, checkloss=True):
    failed=[]
    messages=[]

    for errfile in os.listdir(errfolder):
        jobid = errfile.split(".")[0]
        if checkloss:
            if not os.path.exists(errfolder+"/../losses/"+jobid+".tfs"):
                failed += [jobid]
                if "Missing lossfile" not in messages:
                    messages += ["Missing lossfile"]
                continue
        if checkerr:
            if os.stat(errfolder+'/'+errfile).st_size > 0:
                failed += [jobid]
                with open(errfolder+'/'+errfile, 'r') as f:
                    firstline = f.readline()
                if firstline not in messages:
                    messages += [firstline]

    return failed, messages


def trackplot(trackfolder, obsloc="obs0001", xax="X", yax="PX", cax="TURN", xlim=None, ylim=None, clim=[None,None], tpt=3, batches=None, save=None):
    if xax=='S':
        lossloc = None
    xdata = []
    ydata = []
    cdata = []

    empty = 0

    if batches is None:
        for trackbatch in os.listdir(trackfolder):
            if os.path.isdir(trackfolder+'/'+trackbatch):
                for trackfile in os.listdir(trackfolder+'/'+trackbatch):
                    if obsloc in trackfile.split('.'):
                        trackpath = trackfolder+'/'+trackbatch+'/'+trackfile
                        if os.stat(trackpath).st_size > 0:
                            _, tracktable = readtfs(trackpath)
                            if tpt is not None:
                                xdata += tracktable[xax].tolist()[-tpt:]
                                ydata += tracktable[yax].tolist()[-tpt:]
                                cdata += tracktable[cax].tolist()[-tpt:]
                            else:
                                xdata += tracktable[xax].tolist()
                                ydata += tracktable[yax].tolist()
                                cdata += tracktable[cax].tolist()
                        else:
                            empty += 1
    else:
        for batchnum in batches:
            if os.path.isdir(trackfolder+'/'+str(batchnum)):
                for trackfile in os.listdir(trackfolder+'/'+str(batchnum)):
                    if obsloc in trackfile.split('.'):
                        trackpath = trackfolder+'/'+str(batchnum)+'/'+trackfile
                        if os.stat(trackpath).st_size > 0:
                            _, tracktable = readtfs(trackpath)
                            if tpt is not None:
                                xdata += tracktable[xax].tolist()[-tpt:]
                                ydata += tracktable[yax].tolist()[-tpt:]
                                cdata += tracktable[cax].tolist()[-tpt:]
                            else:
                                xdata += tracktable[xax].tolist()
                                ydata += tracktable[yax].tolist()
                                cdata += tracktable[cax].tolist()
                        else:
                            empty += 1
            else:
                print "tracks for batch",str(batchnum),"not found"

    if empty>0:
        print "warning: "+str(empty)+" empty track files found!"

    xunit = _units[xax]
    yunit = _units[yax]
    
    cm = plt.cm.get_cmap('viridis')
    fig, ax = plt.subplots()
    plt.autoscale(enable=True, axis='both', tight=True)
    plot = ax.scatter(xdata, ydata, c=cdata, cmap=cm, vmin=clim[0], vmax=clim[1], edgecolor='')
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)
    plt.colorbar(plot)

    ax.set_xlabel(xax+' ['+xunit+']')
    ax.set_ylabel(yax+' ['+yunit+']')

    if save is None:
        plt.show()
    else:
        plt.savefig(save)
        plt.close()


def fullplot(folder, lossloc="AP.UP.ZS21633", obsloc="obs0002", xax="X", yax="PX", cax="TURN", xlim=None, ylim=None, clim=[None,None], tpt=3, save=None):
    trackfolder = folder+"/tracks"
    lossfolder = folder+"/losses"
    if xax=='S':
        lossloc = None
    xdata = []
    ydata = []
    cdata = []

    empty = 0

    for lossfile in os.listdir(lossfolder):
        if os.stat(lossfolder+'/'+lossfile).st_size > 0:
            _, losstable = readtfs(lossfolder+'/'+lossfile)
            if lossloc is not None:
                losstable = losstable.loc[losstable['ELEMENT'] == lossloc]
            xdata += losstable[xax].tolist()
            ydata += losstable[yax].tolist()
            cdata += losstable[cax].tolist()
        else:
            empty += 1

    if empty>0:
        print "warning: "+str(empty)+" empty loss files found!"

    empty = 0

    for trackbatch in os.listdir(trackfolder):
        if os.path.isdir(trackfolder+'/'+trackbatch):
            for trackfile in os.listdir(trackfolder+'/'+trackbatch):
                if obsloc in trackfile.split('.'):
                    trackpath = trackfolder+'/'+trackbatch+'/'+trackfile
                    if os.stat(trackpath).st_size > 0:
                        _, tracktable = readtfs(trackpath)
                        if tpt is not None:
                            xdata += tracktable[xax].tolist()[-tpt:]
                            ydata += tracktable[yax].tolist()[-tpt:]
                            cdata += tracktable[cax].tolist()[-tpt:]
                        else:
                            xdata += tracktable[xax].tolist()
                            ydata += tracktable[yax].tolist()
                            cdata += tracktable[cax].tolist()
                    else:
                        empty += 1

    if empty>0:
        print "warning: "+str(empty)+" empty track files found!"

    xunit = _units[xax]
    yunit = _units[yax]
    
    cm = plt.cm.get_cmap('viridis')
    fig, ax = plt.subplots()
    plt.autoscale(enable=True, axis='both', tight=True)
    plot = ax.scatter(xdata, ydata, c=cdata, cmap=cm, vmin=clim[0], vmax=clim[1], edgecolor='')
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)
    plt.colorbar(plot)

    ax.set_xlabel(xax+' ['+xunit+']')
    ax.set_ylabel(yax+' ['+yunit+']')

    if save is None:
        plt.show()
    else:
        plt.savefig(save)
        plt.close()


def steinbachplot(folder, xlim=None, ylim=None, save=None):
    trackfolder = folder+"/tracks"
    lossfolder = folder+"/losses"
    xdata_in = []
    ydata_in = []
    xdata_out = []
    ydata_out = []

    empty = 0

    _, twisstable = readtfs(folder+"/thin_twiss.tfs")
    xco = twisstable.loc["BEGI.10010", "X"]
    pxco = twisstable.loc["BEGI.10010", "PX"]
    alpha = twisstable.loc["BEGI.10010", "ALFX"]
    beta = twisstable.loc["BEGI.10010", "BETX"]
    gamma = (1+alpha**2)/beta

    for trackbatch in os.listdir(trackfolder):
        if os.path.isdir(trackfolder+'/'+trackbatch):
            _, losstable = readtfs(lossfolder+'/'+trackbatch+".tfs")
            lostparticles = losstable.index.values.tolist()
            for trackfile in os.listdir(trackfolder+'/'+trackbatch):
                if trackfile.split('.')[2]=="obs0001":
                    trackpath = trackfolder+'/'+trackbatch+'/'+trackfile
                    if os.stat(trackpath).st_size > 0:
                        _, tracktable = readtfs(trackpath)
                        x = tracktable["X"].iloc[0]-xco
                        px = tracktable["PX"].iloc[0]-pxco
                        if int(trackfile.split('.')[-1][-4:]) in lostparticles:
                            # TODO: USE A COLORMAP BASED ON TURN FOR THESE
                            xdata_out += [tracktable["PT"].iloc[0]]
                            ydata_out += [gamma*x**2 + beta*px**2 + alpha*x*px]
                        else:
                            xdata_in += [tracktable["PT"].iloc[0]]
                            ydata_in += [gamma*x**2 + beta*px**2 + alpha*x*px]
                    else:
                        empty += 1

    if empty>0:
        print "warning: "+str(empty)+" empty track files found!"

    fig, ax = plt.subplots()
    plt.autoscale(enable=True, axis='both', tight=True)
    # TODO: Add histograms of extracted/non-extracted?
    plot = ax.scatter(xdata_in, ydata_in, color="k")
    plot = ax.scatter(xdata_out, ydata_out, color="r")
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)

    ax.set_xlabel("PT [1]")
    ax.set_ylabel("amplitude (disregarding dispersion...) [1]")

    if save is None:
        plt.show()
    else:
        plt.savefig(save)
        plt.close()


def losshistscatter(lossfolder, lossloc="AP.UP.ZS21633",
                    xax="X", yax="PX", cax="TURN",
                    xlim=None, ylim=None, clim=[None,None],
                    monochrom=False,
                    xbin=None, ybin=None,
                    log=False, extra=None, save=None,
                    datalim=[[None,None],[None,None],[None,None]]):
    if xax=='S':
        lossloc = None
    xdata = []
    ydata = []
    cdata = []

    empty = 0

    for lossfile in os.listdir(lossfolder):
        if os.stat(lossfolder+'/'+lossfile).st_size > 0:
            _, losstable = readtfs(lossfolder+'/'+lossfile)
            if lossloc is not None:
                losstable = losstable.loc[losstable['ELEMENT'] == lossloc]
            xdata += losstable[xax].tolist()
            ydata += losstable[yax].tolist()
            cdata += losstable[cax].tolist()
        else:
            empty += 1

    if empty>0:
        print "warning: "+str(empty)+" empty loss files found!"

    selector = [((True if (datalim[0][0] is None) else (datalim[0][0]<=xdata[i])) and
                 (True if (datalim[0][1] is None) else (xdata[i]<=datalim[0][1])) and
                 (True if (datalim[1][0] is None) else (datalim[1][0]<=ydata[i])) and
                 (True if (datalim[1][1] is None) else (ydata[i]<=datalim[1][1])) and
                 (True if (datalim[2][0] is None) else (datalim[2][0]<=cdata[i])) and
                 (True if (datalim[2][1] is None) else (cdata[i]<=datalim[2][1])))
                 for i in range(len(xdata))]
    xdata = [xdata[i] for i in range(len(selector)) if selector[i]]
    ydata = [ydata[i] for i in range(len(selector)) if selector[i]]
    cdata = [cdata[i] for i in range(len(selector)) if selector[i]]


    xunit = _units[xax]
    yunit = _units[yax]

    left, width = 0.1, 0.65
    left_h = left + width + 0.02
    bottom, height = 0.1, 0.65
    bottom_h = bottom + height + 0.02

    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom_h, width, 0.2]
    rect_histy = [left_h, bottom, 0.2, height]

    plt.figure(1, figsize=(8, 8))
    cm = plt.cm.get_cmap('viridis')

    axScatter = plt.axes(rect_scatter)
    #plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
    #plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))

    axHistx = plt.axes(rect_histx)
    axHisty = plt.axes(rect_histy)

    axHistx.xaxis.set_major_formatter(NullFormatter())
    axHisty.yaxis.set_major_formatter(NullFormatter())
    if monochrom:
        axScatter.scatter(xdata, ydata, edgecolor='')
    else:
        axScatter.scatter(xdata, ydata, c=cdata, cmap=cm, vmin=clim[0], vmax=clim[1], edgecolor='')

    if extra is not None:
        for line in extra:
            linex = [v[0] for v in line]
            liney = [v[1] for v in line]
            axScatter.plot(linex, liney, 'k--')

    weights = 100*np.ones_like(xdata)/len(xdata)

    if xbin is None:
        xbin = (max(xdata)-min(xdata))/100
    if ybin is None:
        ybin = (max(xdata)-min(xdata))/100

    if xlim is None:
        binsx = np.arange(min(xdata)-xbin, max(xdata)+2*xbin, xbin)
    else:
        binsx = np.arange(xlim[0]-xbin, xlim[1]+xbin, xbin)

    if ylim is None:
        binsy = np.arange(min(ydata)-ybin, max(ydata)+2*ybin, ybin)
    else:
        binsy = np.arange(ylim[0]-ybin, ylim[1]+ybin, ybin)
        
    axScatter.set_xlim((binsx[0], binsx[-1]))
    axScatter.set_ylim((binsy[0], binsy[-1]))

    axHistx.hist(xdata, bins=binsx, weights=weights, log=log)
    axHisty.hist(ydata, bins=binsy, weights=weights,
                 orientation='horizontal', log=log)

    axHistx.set_xlim(axScatter.get_xlim())
    axHisty.set_ylim(axScatter.get_ylim())

    axScatter.set_xlabel(xax+' ['+xunit+']')
    axScatter.set_ylabel(yax+' ['+yunit+']')

    axHistx.set_ylabel(r'$\%$')
    axHisty.set_xlabel(r'$\%$')

    #plt.colorbar(axScatter) # Does this work??

    if save is None:
        plt.show()
    else:
        plt.savefig(save)
        plt.close()














# HEAVILY CUSTOM STUFF THAT SHOULDN'T REALLY BE IN THIS LIBRARY
def zsbacktrack(trackfolder, cax="TURN", save_d=None, save_u=None):
    xax="X"
    yax="PX"
    xunit = _units[xax]
    yunit = _units[yax]
    cm = plt.cm.get_cmap('viridis')

    xdata = []
    ydata = []
    cdata = []
    iddata = []

    xdata_u = []
    ydata_u = []
    cdata_u = []

    # Read in downstream ZS data
    empty = 0
    obsloc = 'obs0003'

    for trackbatch in os.listdir(trackfolder):
        if os.path.isdir(trackfolder+'/'+trackbatch):
            for trackfile in os.listdir(trackfolder+'/'+trackbatch):
                if obsloc in trackfile.split('.'):
                    trackpath = trackfolder+'/'+trackbatch+'/'+trackfile
                    if os.stat(trackpath).st_size > 0:
                        _, tracktable = readtfs(trackpath)
                        xdata += tracktable[xax].tolist()[-1:]
                        ydata += tracktable[yax].tolist()[-1:]
                        cdata += tracktable[cax].tolist()[-1:]
                        iddata += [(trackfile.split('.')[1],trackfile.split('.')[3],tracktable["TURN"].tolist()[-1:])]
                    else:
                        empty += 1
    if empty>0:
        print "warning: "+str(empty)+" empty track files found!"

    # Find the weird region and plot downstream
    weird = [i for i in range(len(xdata)) if (xdata[i]>0.06 and ydata[i]<-0.001)]

    xdata = [xdata[i] for i in weird]
    ydata = [ydata[i] for i in weird]
    cdata = [cdata[i] for i in weird]
    iddata = [iddata[i] for i in weird]

    fig, ax = plt.subplots()
    plt.autoscale(enable=True, axis='both', tight=True)
    plot = ax.scatter(xdata, ydata, c=cdata, cmap=cm, edgecolor='')
    plt.colorbar(plot)

    ax.set_xlabel(xax+' ['+xunit+']')
    ax.set_ylabel(yax+' ['+yunit+']')

    if save_d is None:
        plt.show()
    else:
        plt.savefig(save_d)
        plt.close()

    # Read the corresponding upstream locations and plot
    for tag in iddata:
        trackbatch = tag[0][5:]
        trackpath = trackfolder+'/'+trackbatch+'/'+'.'.join(("track",tag[0],"obs0002",tag[1]))
        if os.stat(trackpath).st_size > 0:
            _, tracktable = readtfs(trackpath)
            xdata_u += tracktable[xax].tolist()[-1:]
            ydata_u += tracktable[yax].tolist()[-1:]
            cdata_u += tracktable[cax].tolist()[-1:]
        else:
            print "file "+trackpath+" not found"

    fig, ax = plt.subplots()
    plt.autoscale(enable=True, axis='both', tight=True)
    plot = ax.scatter(xdata_u, ydata_u, c=cdata_u, cmap=cm, edgecolor='')
    plt.colorbar(plot)

    ax.set_xlabel(xax+' ['+xunit+']')
    ax.set_ylabel(yax+' ['+yunit+']')

    if save_u is None:
        plt.show()
    else:
        plt.savefig(save_u)
        plt.close()

    # Sanity check
    if cdata==cdata_u:
        print "Sanity check ok"
    else:
        print "NOT OK..."
    



























# BELOW IS SOME OLD CODE I MIGHT RE-USE, DO NOT RELY ON IT...


class TrackOut:
    """Python representation of tracking output

    Expect MAD-X data generated without the "ONEFILE" option
    """
    def __init__(self, madxoutput, obspoint="0001", particles=[5,4,3,2,1]):
        self.units = {'X': 'm', 'PX': 'rad', 'Y': 'm', 'PY': 'rad', 'T': 'm', 'PT': '1', 'S': 'm', 'E': '???'}
        with open(madxoutput, 'r') as file:
            line = file.readline()
            while not line.startswith('*'):
                line = file.readline()
            keys = line.split()[3:]
            line = file.readline()
            line = file.readline()
            self.nturns = int(line.split()[2])-1
            self.npart = int(line.split()[3])

            coords = np.empty([self.npart, self.nturns, len(keys)])
            coords[:] = np.nan

            obs = False

            for line in file:
                if line.startswith('#'):
                    if line.split()[5] == obspoint:
                        obs = True
                    else:
                        obs = False
                elif obs == True :
                    partid = int(line.split()[0])
                    turnid = int(line.split()[1])
                    coords[partid-1,turnid-1,:] = line.split()[2:]

            self.coords = pd.Panel(coords,
                items = range(self.npart),
                major_axis = range(self.nturns),
                minor_axis = keys)

            self.min = {key: self.coords.minor_xs(key).min().min() for key in keys}
            self.max = {key: self.coords.minor_xs(key).max().max() for key in keys}














    def plot(self, xax, yax, pids="All", turnstart=0, turnend=None, save=None):
        """Make a plot of (xax, yax) tracking data (e.g. ('X','PX'))

        pids is a list of particle IDs to be plotted, or the string
        "All" to plot all particles. Particle IDs start at 1, not 0, in
        order to conform to the MAD-X standard.
        """
        (xexp, xunit) = siprefix(self.min[xax], self.max[xax], self.units[xax])
        (yexp, yunit) = siprefix(self.min[yax], self.max[yax], self.units[yax])

        if turnend is None:
            turnend = self.nturns

        if pids == "All":
            pids = range(self.npart)
        else:
            pids = [pid-1 for pid in pids]
        colors = ['r', 'y', 'g', 'c', 'b', 'm', 'k']
        ncolors = len(colors)
        markers = ['v','s','H','d','^','o','*','<','8','D','>','p','h']
        nmarkers = len(markers)            

        fig, ax = plt.subplots()
        for pid in pids:
            ax.plot(self.coords[pid][xax][turnstart:turnend]/10**xexp,
                    self.coords[pid][yax][turnstart:turnend]/10**yexp,
                    color = colors[pid%ncolors],
                    label = "particle "+str(pid+1),
                    linestyle = 'None', markeredgecolor = 'none',
                    marker = markers[pid%nmarkers])

        ax.set_xlabel(xax+' ['+xunit+']')
        ax.set_ylabel(yax+' ['+yunit+']')
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        ax.legend(numpoints=1, loc='center left', bbox_to_anchor=(1, 0.5))
        plt.show() #save it first, add option to show

    def movie(self, xax, yax, pids="All", nhist="All", save="test", fps=None, dpi=None):
        """Make a plot of (xax, yax) tracking data (e.g. ('X','PX'))

        Like plot, but it moves :)
        """
        (xexp, xunit) = siprefix(self.min[xax], self.max[xax], self.units[xax])
        (yexp, yunit) = siprefix(self.min[yax], self.max[yax], self.units[yax])

        if nhist == "All":
            nhist = self.nturns

        if pids == "All":
            pids = range(self.npart)
        else:
            pids = [pid-1 for pid in pids]
        colors = ['r', 'y', 'g', 'c', 'b', 'm', 'k']
        ncolors = len(colors)
        markers = ['v','s','H','d','^','o','*','<','8','D','>','p','h']
        nmarkers = len(markers)

        fig, ax = plt.subplots()

        lines = [plt.plot([], [], color = colors[pid%ncolors],
                    label = "particle "+str(pid+1), linestyle = 'None',
                    markeredgecolor = 'none',
                    marker = markers[pid%nmarkers])[0] for pid in pids]
        ax.set_xlim(self.min[xax]/10**xexp, self.max[xax]/10**xexp)
        ax.set_ylim(self.min[yax]/10**yexp, self.max[yax]/10**yexp)
        ax.set_xlabel(xax+' ['+xunit+']')
        ax.set_ylabel(yax+' ['+yunit+']')
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width*0.8, box.height])
        ax.legend(numpoints=1, loc='center left', bbox_to_anchor=(1, 0.5))
        ax.set_title('Turn = ')

        def init():    
            for line in lines:
                line.set_data([], [])
            return lines

        def animate(i):
            for j,line in enumerate(lines):
                line.set_data(self.coords[
                              pids[j]][xax][max([0,i-nhist]):i]/10**xexp,
                              self.coords[
                              pids[j]][yax][max([0,i-nhist]):i]/10**yexp)
            ax.set_title('Turn = '+str(i))
            return lines

        anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=self.nturns+1, interval=1, blit=False)
        anim.save(save+'.mp4', fps=fps, dpi=dpi)

        plt.show() #save it first, add option to show



# Put into seperate library?
def siprefix(datamin, datamax, unit_in):
    sismall = ['', 'm', '$\mu$', 'n']
    exponent = 0
    unitout = unit_in
    for i,_ in enumerate(sismall):
        if (datamax<10*10**(-3*i) and
              datamin>-10*10**(-3*i)):
            exponent = -3*i
            if unit_in=='1':
                unit_out = '1E'+str(-1*exponent)
            else:
                unit_out = str(sismall[i])+unit_in
    return(exponent, unit_out)
