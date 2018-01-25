# -*- coding: utf-8 -*-
"""
Ugly script to make plots for the amplitude extraction.
"""
import os
import sys
import subprocess
import matplotlib as mpl
mpl.use('Agg') #To use plotting routines on machines without display
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
import numpy as np
sys.path.insert(1, '/afs/cern.ch/project/sloex/code/madxBatch')
import python.dataprocessing as datproc

def losshistscatter_alt(lossfolder, lossloc="AP.UP.ZS21633",
                        xax="X", yax="PX", cax="TURN", otherax="PT",
                        xlim=None, ylim=None, clim=[None,None],
                        monochrom=False,
                        xbin=None, ybin=None,
                        log=False, extra=None, save=None,
                        datalim=[[None,None],[None,None],[None,None]],
                        ofilter=(lambda x: True)):
    if xax=='S':
        lossloc = None
    xdata = []
    ydata = []
    cdata = []
    odata = []

    empty = 0

    for lossfile in os.listdir(lossfolder):
        if os.stat(lossfolder+'/'+lossfile).st_size > 0:
            _, losstable = datproc.readtfs(lossfolder+'/'+lossfile)
            if lossloc is not None:
                losstable = losstable.loc[losstable['ELEMENT'] == lossloc]
            xdata += losstable[xax].tolist()
            ydata += losstable[yax].tolist()
            cdata += losstable[cax].tolist()
            odata += losstable[otherax].tolist()
        else:
            empty += 1

    if empty>0:
        print "warning: "+str(empty)+" empty loss files found!"

    selector = [((True if (datalim[0][0] is None) else (datalim[0][0]<=xdata[i])) and
                 (True if (datalim[0][1] is None) else (xdata[i]<=datalim[0][1])) and
                 (True if (datalim[1][0] is None) else (datalim[1][0]<=ydata[i])) and
                 (True if (datalim[1][1] is None) else (ydata[i]<=datalim[1][1])) and
                 (True if (datalim[2][0] is None) else (datalim[2][0]<=cdata[i])) and
                 (True if (datalim[2][1] is None) else (cdata[i]<=datalim[2][1])) and
                 ofilter(odata[i]))
                 for i in range(len(xdata))]
    xdata = [xdata[i] for i in range(len(selector)) if selector[i]]
    ydata = [ydata[i] for i in range(len(selector)) if selector[i]]
    cdata = [cdata[i] for i in range(len(selector)) if selector[i]]


    xunit = datproc._units[xax]
    yunit = datproc._units[yax]

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






folder = sys.argv[1]
lossfolder = folder+'/losses'
trackfolder = folder+'/tracks'
errfolder = folder+'/error'
name = folder.split("/")[-1]
if name=="":
    name = folder.split("/")[-2]

lossloc = "AP.UP.ZS21633"

failed, messages = datproc.errorcheck(errfolder)
print (str(len(failed))+" processes failed")
print "errors"
for message in messages:
    print message
print ""

plotfolder = folder+"/../plots/"+name

if not os.path.exists(plotfolder):
    os.makedirs(plotfolder)

if not os.path.exists(trackfolder+'/0'):
    subprocess.check_call('./unpacker.sh', shell=True, cwd=trackfolder)

# Print efficiency and loss stats
stdout = sys.stdout
with open(plotfolder+"/stats.txt", 'w') as sys.stdout:
    datproc.efficiency(lossfolder, aperturex=[0.06815,0.08815], aperturey=[-0.023,0.023], zs_len=18.77, zs_an=4.1635E-4, aperturex2=[0.04140,0.06140])
    print ""
    datproc.wireangle(lossfolder)
    print ""
    datproc.emittance(lossfolder, lossloc=lossloc, ap=[0.04219,0.08219], betagamma=426.3167)
sys.stdout = stdout

# Make plots
datproc.lossplot(lossfolder, lossloc=lossloc, xax='TURN', yax='PT', cax='PT', ylim=[-0.0025, 0.0020], clim=[-0.0025, 0.0020], save=plotfolder+"/sweep_"+name+".png")
datproc.trackplot(trackfolder, obsloc="obs0002", cax="PT", clim=[-0.0025, 0.0020], tpt=1, save=plotfolder+"/zs_up_"+name+".png")
datproc.trackplot(trackfolder, obsloc="obs0003", cax="PT", clim=[-0.0025, 0.0020], tpt=1, save=plotfolder+"/zs_down_"+name+".png")
datproc.trackplot(trackfolder, obsloc="obs0004", cax="PT", clim=[-0.0025, 0.0020], tpt=1, save=plotfolder+"/tpst_circ_"+name+".png")

datproc.losshistscatter(lossfolder, lossloc=lossloc, xlim=[0.06815,0.08515], ylim=[-0.00185, -0.00135],cax="PT", clim=[-0.0025, 0.0020], xbin=0.0002, ybin=0.00001, log=False, save=plotfolder+"/losshist_"+name+".png")
datproc.losshistscatter(lossfolder, lossloc=lossloc, xlim=[0.06815,0.08515], ylim=[-0.00185, -0.00135],cax="TURN", clim=[0,50000], xbin=0.0002, ybin=0.00001, log=False, save=plotfolder+"/losshist_turn_"+name+".png")
datproc.losshistcombo(lossfolder, lossloc=lossloc, xlim=[0.06815,0.08515], ylim=[-0.00185, -0.00135], xbin=0.0002, ybin=0.000005, cm='viridis', log=False, save=plotfolder+"/losshist2_"+name+".png")
datproc.losshistcombo(lossfolder, lossloc=lossloc, xlim=[0.06815,0.07215], ylim=[-0.00155, -0.0014], xbin=0.00008, ybin=0.000004, cm='viridis', log=False, save=plotfolder+"/losshist2_zoom_"+name+".png")

res = 5000

lim1 = 0
lim2 = res
while lim1 < 50000:
    datproc.losshistscatter(lossfolder, lossloc=lossloc, xlim=[0.06815,0.08515], ylim=[-0.00185, -0.00135],cax="TURN", clim=[0,50000], xbin=0.0002, ybin=0.00001, log=False, save=plotfolder+"/losshist_turn_"+str(lim1)+"_"+str(lim2)+"_"+name+".png", datalim=[[None,None],[None,None],[lim1,lim2]])
    losshistscatter_alt(lossfolder, lossloc=lossloc, xlim=[0.06815,0.08515], ylim=[-0.00185, -0.00135],cax="PT", otherax="TURN" , clim=[-0.0025, 0.0020], xbin=0.0002, ybin=0.00001, log=False, save=plotfolder+"/losshist_pt_turn_"+str(lim1)+"_"+str(lim2)+"_"+name+".png",
                        datalim=[[None,None],[None,None],[None,None]], ofilter=(lambda x: (lim1<=x and x<=lim2)))
    lim1 += res
    lim2 += res

losshistscatter_alt(lossfolder, lossloc=lossloc, xlim=[0.06815,0.08515], ylim=[-0.00185, -0.00135],cax="TURN", otherax="TURN" , clim=[0,50000], xbin=0.0002, ybin=0.00001, log=False, save=plotfolder+"/losshist_movement_"+name+".png",
                    datalim=[[None,None],[None,None],[None,None]], ofilter=(lambda x: ((15000<=x and x<=20000) or (x>=45000))))



