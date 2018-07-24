# -*- coding: utf-8 -*-
"""
Script to make nominal plots

@author: Linda Stoel
"""
import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import os
import sys

import python.dataprocessing as datproc

ap = datproc.nomap

#TODO python 2/3 compatibility

# Get command line args
folder = sys.argv[1]
try:
    kind = sys.argv[2]
except IndexError:
    kind = '01'
    

# Check for errors
failed, messages = datproc.errorcheck(folder+'/error')
print (str(len(failed))+" processes failed")
if not len(failed)==0:
    print "Errors:"
    for message in messages:
        print message
    print ""

# Get study name
name = folder.split("/")[-1]
if name=="":
    name = folder.split("/")[-2]

# Make plot folder
plotfolder = folder+"/../plots/"+name
if not os.path.exists(plotfolder):
    os.makedirs(plotfolder)

# Get settings
settings = datproc.getsettings(folder)
sliced = False if settings['slices']=='None' else True
pycoll = eval(settings['pycollimate'])
myloc = "AP.UP.TPST21760" if pycoll else "AP.UP.ZS21633"

if '0' in kind or '1' in kind:
    # Get losses
    losses = datproc.getlosses(folder+'/losses', settings=settings)
    myloss = losses[(losses['ELEMENT']==myloc)]

if '0' in kind:
    # Print efficiency and loss stats
    if pycoll:
        myap = {'aperturex': ap['tpstcirc']+ap['tpstblade']+np.array([0, ap['tpstex']]),
                'aperturey': [-0.0100,0.0100]}
    else:
        myap = {'aperturex': ap['zsupmid']+ap['zsthick']/2+np.array([0, ap['zsex']]),
                'aperturex2': ap['zsdomid']+ap['zsthick']/2+np.array([0, ap['zsex']]),
                'aperturey': [-0.023,0.023],
                'zs_len': 18.77,
                'zs_an': 4.1635E-4}
    stdout = sys.stdout
    with open(plotfolder+"/stats.txt", 'w') as sys.stdout:
        print 'Probabilities with 95% confidence intervals:'
        alltags,_ = datproc.extr_tagger(data=losses, pycoll=pycoll, **myap)
        datproc.efficiency(losses, alltags=alltags, silent=False)
        print '\nLoss stats:'
        datproc.lossstats(losses, merge=False, silent=False)
        print 'Beam stats:'
        datproc.beamstats(losses[losses['tag']=='extracted'], silent=False)
    sys.stdout = stdout

if '1' in kind:
    # Make loss plots
    datproc.plotter(myloss, xax="TURN", yax="PT", cax='PT', kind='hexbin',
                    ylim=[-0.0025, 0.0020], clim=[-0.0025, 0.0020],
                    mainkwargs={'bins': 'log'},
                    save=plotfolder+"/sweep_"+name+'.png')
    if pycoll:
        for plotkind in ['scatter', 'hist2d']:
            datproc.plotter(myloss, xax='X', yax='PX', cax='PT', kind=plotkind,
                            xlim=[ap['tpstcirc'], 0.095], ylim=[0.0003,0.0027],
                            clim=[-0.0025, 0.0020], xbin=0.0005, ybin=0.000025, log=True,
                            save=plotfolder+'/tpst_loss_'+plotkind+'_h_'+name+'.png')
            datproc.plotter(myloss, xax='Y', yax='PY', cax='PT', kind=plotkind,
                            xlim=[-0.02, 0.02], ylim=[-0.0006,0.0006],
                            clim=[-0.0025, 0.0020], xbin=0.0003, ybin=0.000015, log=True,
                            save=plotfolder+'/tpst_loss_'+plotkind+'_v_'+name+'.png')
            datproc.plotter(myloss, xax='X', yax='Y', cax='PT', kind=plotkind,
                            xlim=[ap['tpstcirc'], 0.095], ylim=[-0.02, 0.02],
                            clim=[-0.0025, 0.0020], xbin=0.0005, ybin=0.0003, log=True,
                            save=plotfolder+'/tpst_loss_'+plotkind+'_s_'+name+'.png')
        myrange = [1668,1688]
        datproc.plotter(losses[losses['S'].between(*myrange)], xax='S', yax='X',
                        cax="PT", xlim=myrange, clim=[-0.0025, 0.0020], 
                        ylim=[ap['zsdomid']-ap['zsthick']/2-ap['zsex'],
                              ap['zsupmid']+ap['zsthick']/2+ap['zsex']],
                        log=True, save=plotfolder+"/zs_loss_"+name+".png")
    else:
        lim_pt = [-0.0025, 0.0020]
        lim_x = ap['zsupmid']+np.array([-ap['zsthick']/2, ap['zsthick']+ap['zsex']])
        lim_px = [-0.00195, -0.00135]
        datproc.plotter(myloss, xax='X', yax='PX', cax='PT',
                        xlim=lim_x, ylim=lim_px, clim=lim_pt,
                        xbin=0.0002, ybin=0.000005, log=False,
                        save=plotfolder+"/losshist_"+name+".png")
        datproc.plotter(myloss, xax='X', yax='PX', cax='PT', kind='hist2d',
                        xlim=lim_x, ylim=lim_px, clim=lim_pt,
                        xbin=0.0002, ybin=0.000005, log=False,
                        save=plotfolder+"/losshist2_"+name+".png")
        datproc.plotter(myloss, xax='X', yax='PX', cax='PT', kind='hist2d',
                        xlim=[lim_x[0],0.07215], ylim=[-0.00155, -0.0014], clim=lim_pt,
                        xbin=0.00008, ybin=0.000004, log=False,
                        save=plotfolder+"/losshist2_zoom_"+name+".png")

if '2' in kind:
    # Get track data and make track plots
    tracklocs = eval(settings['elements'])
    tpt = None if sliced else 9
    for trackloc in tracklocs:
        tracks = datproc.gettracks(folder+'/tracks', settings=settings, obsloc=trackloc)
        for obsnum in [1,3,6,9]:
            datproc.plotter(tracks[tracks['obsnum']>-obsnum], cax="PT",
                            clim=[-0.0025, 0.0020],
                            save=plotfolder+'/'+trackloc+'_'+str(obsnum)+'obs_'+name+'.png')
        if tpt is None:
            datproc.plotter(tracks, cax="PT", clim=[-0.0025, 0.0020],
                            save=plotfolder+'/'+trackloc+'_full_'+name+'.png')
