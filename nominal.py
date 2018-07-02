# -*- coding: utf-8 -*-
"""
Nominal SPS slow extraction with and without scattering

@author: Linda Stoel
"""
import os
from python.batching import Settings, submit_job
import sys

bl = [True, False]
modes = ['qsweep', 'cose', 'ampex']

# Simulations with sweep
for mode in modes:
    for pc, db, fast in [(x,y,z) for x in bl for y in bl for z in bl]:
        if pc and not fast:
            continue

        if db and mode=='cose':
            continue

        name = mode
        name += "_db" if db else "_nodb"
        name += "_fast" if fast else ""
        name += "_pc" if pc else "_nopc"

        settings=Settings(name, studygroup="nominal", disk='afsproject')
        settings.seed = 0
        settings.pycollimate = pc
        settings.dynamicbump = db

        if mode=='cose':
            settings.cose = True
        elif mode=='ampex':
            here = os.path.dirname(os.path.realpath(sys.argv[0]))
            settings.finalchanges = here+'/madx/ampex_finalchanges.cmdx'
            settings.ampex = True

        if fast:
            settings.nturns = 50000
            settings.ffile = 500
            if pc:
                settings.flavour = "testmatch"
        else:
            settings.nturns = 204565
            settings.ffile = 2045

        if pc:
            settings.elements = ['AP.UP.ZS21633_M','AP.DO.ZS21676_M','AP.UP.TPST21760']
            settings.nbatches = 100
            settings.nparperbatch = 2500
        else:
            settings.elements = ['AP.UP.ZS21633','AP.DO.ZS21676','AP.UP.TPST21760']
            settings.nbatches = 500
            settings.nparperbatch = 500

        submit_job(settings)

# Simulations of momentum slices to be fixed in future release
