import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.lines import Line2D

_units = {'X': 'm', 'PX': 'rad', 'Y': 'm', 'PY': 'rad', 'T': 'm', 'PT': '1', 'S': 'm', 'E': 'eV'}
_def_unit_exp = {'X': 0, 'PX': 0, 'Y': 0, 'PY': 0, 'T': 0, 'PT': 0, 'S': 0, 'E': 9}

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

    return header, table


def readsingletrack(filename):
    header, table = readtfs(filename, usecols=[1,2,3,4,5,6,7,8,9])
    return header, table


# TODO: change scale
def lossplot(lossfolder, lossloc="AP.UP.ZS21633", xax="X", yax="PX", cax="TURN", xlim=None, ylim=None, save=None):
    if xax=='S':
        lossloc = None
    xdata = []
    ydata = []
    cdata = []

    for lossfile in os.listdir(lossfolder):
        _, losstable = readtfs(lossfolder+'/'+lossfile)
        if lossloc is not None:
            losstable = losstable.loc[losstable['ELEMENT'] == lossloc]
        xdata += losstable[xax].tolist()
        ydata += losstable[yax].tolist()
        cdata += losstable[cax].tolist()

    xunit = _units[xax]
    yunit = _units[yax]
    
    cm = plt.cm.get_cmap('viridis')
    fig, ax = plt.subplots()
    plt.autoscale(enable=True, axis='both', tight=True)
    plot = ax.scatter(xdata, ydata, c=cdata, cmap = cm)
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
