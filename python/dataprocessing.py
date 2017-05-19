import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.lines import Line2D
from matplotlib.ticker import NullFormatter

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

    return header, table


def readsingletrack(filename):
    header, table = readtfs(filename, usecols=[1,2,3,4,5,6,7,8,9])
    return header, table


# TODO: change scale
def lossplot(lossfolder, lossloc="AP.UP.ZS21633", xax="X", yax="PX", cax="TURN", xlim=None, ylim=None, clim=None, save=None):
    if xax=='S':
        lossloc = None
    if clim is None:
        clim=(None,None)
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
def lossstats(lossfolder):
    data = {}

    for lossfile in os.listdir(lossfolder):
        if os.stat(lossfolder+'/'+lossfile).st_size > 0:
            _, losstable = readtfs(lossfolder+'/'+lossfile)
            for location in losstable['ELEMENT'].tolist():
                try:
                    data[location]+=1
                except KeyError:
                    data[location]=1

    print "Start loss stats:"
    for location in sorted(data, key=data.get, reverse=True):
        print location, data[location]
    print "End loss stats."


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
    statemit = str(np.sqrt(np.linalg.det(np.cov(xdata,pdata))))

    if lossloc is None:
        message = "Globally:\n"
    else:
        message = "At "+lossloc+":\n"

    message += ("min,max "+xax+": "+minx+", "+maxx+"\n"+
                "avg,std "+xax+": "+avgx+", "+stdx+"\n"+
                "min,max "+pax+": "+minp+", "+maxp+"\n"+
                "avg,std "+pax+": "+avgp+", "+stdp+"\n"+
                "statistical emittance: "+statemit+"\n"+
                "(# of particles local/global: "+str(len(xdata))+"/"+str(totloss)+")")

    if save is None:
        print message
    else:
        with open(save, 'w') as out:
            out.write(message)


def efficiency(lossfolder, pycoll=False, aperturex=[0,1], aperturey=[0,1]):
    xax='X'
    yax='Y'

    empty = 0

    if pycoll:
        extracted = 0
        tpst = 0
        tce = 0
        zs_wire = 0
        zs_app = 0
        other = 0

        for lossfile in os.listdir(lossfolder):
            if os.stat(lossfolder+'/'+lossfile).st_size > 0:
                _, losstable = readtfs(lossfolder+'/'+lossfile)
                for _, particle in losstable.iterrows():
                    if particle['ELEMENT']=="AP.UP.TPST21760":
                        if (particle["X"]>aperturex[0] and particle["X"]<aperturex[1]
                            and particle["Y"]>aperturey[0] and particle["Y"]<aperturey[1]):
                            extracted += 1
                        else:
                            tpst +=1
                    elif "TPST" in particle['ELEMENT']:
                        tpst += 1
                    elif "TCE" in particle['ELEMENT']:
                        tce += 1
                    elif particle['ELEMENT']==("ZS.21655"):
                        zs_wire += 1
                    elif "ZS" in particle['ELEMENT']:
                        zs_app += 1
                    else:
                        other += 1
            else:
                empty += 1

        if empty > 0:
            print str(empty)+" empty loss files found"
                    
        print (str(extracted)+" extracted\n"
               +str(zs_wire)+" lost on the ZS wires\n"
               +str(zs_app)+" lost on the ZS aperture\n"
               +str(tce)+" lost on the TCE\n"
               +str(tpst)+" lost on the TPST\n"
               +str(other)+" lost elsewhere.")
    else:
        extracted = 0
        zs = 0
        other = 0

        for lossfile in os.listdir(lossfolder):
            if os.stat(lossfolder+'/'+lossfile).st_size > 0:
                _, losstable = readtfs(lossfolder+'/'+lossfile)
                for _, particle in losstable.iterrows():
                    if particle['ELEMENT']=="AP.UP.ZS21633":
                        if (particle["X"]>aperturex[0] and particle["X"]<aperturex[1]
                            and particle["Y"]>aperturey[0] and particle["Y"]<aperturey[1]):
                            extracted += 1
                        else:
                            zs +=1
                    elif "ZS" in particle['ELEMENT']:
                        zs += 1
                    else:
                        other += 1
            else:
                empty += 1

        if empty > 0:
            print str(empty)+" empty loss files found"
                    
        print (str(extracted)+" extracted\n"
               +str(zs)+" lost on the ZS\n"
               +str(other)+" lost elsewhere.")


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
               "avg,std "+pax+": "+avgp+", "+stdp+"\n")

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

def errorcheck(errfolder):
    failed=[]
    messages=[]

    for errfile in os.listdir(errfolder):
        jobid = errfile.split(".")[0]
        if os.stat(errfolder+'/'+errfile).st_size > 0:
            failed += [jobid]
            with open(errfolder+'/'+errfile, 'r') as f:
                firstline = f.readline()
            if firstline not in messages:
                messages += [firstline]

    return failed, messages


def trackplot(trackfolder, obsloc="obs0001", xax="X", yax="PX", cax="TURN", xlim=None, ylim=None, clim=None, tpt=3, batches=None, save=None):
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
    plot = ax.scatter(xdata, ydata, c=cdata, cmap=cm, edgecolor='')
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


def fullplot(folder, lossloc="AP.UP.ZS21633", obsloc="obs0002", xax="X", yax="PX", cax="TURN", xlim=None, ylim=None, clim=None, tpt=3, save=None):
    trackfolder = folder+"/tracks"
    lossfolder = folder+"/losses"
    if clim is None:
        clim=(None,None)
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


def losshistscatter(lossfolder, lossloc="AP.UP.ZS21633",
                    xax="X", yax="PX", cax="TURN",
                    xlim=None, ylim=None, clim=None,
                    monochrom=False,
                    xbin=None, ybin=None,
                    log=False, extra=None, save=None):
    if xax=='S':
        lossloc = None
    if clim is None:
        clim=(None,None)
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

    if log:
        axHistx.set_ylabel(r'$\log(\%)$')
        axHisty.set_xlabel(r'$\log(\%)$')
    else:
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
