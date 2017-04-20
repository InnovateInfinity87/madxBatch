# -*- coding: utf-8 -*-
"""
Nominal SPS slow extraction with and without scattering

@author: Linda Stoel
"""
from python.batching import Settings, submit_job

#Simulate without pycollimate
settings=Settings('nominal_nopycoll', disk='afspublic')

settings.trackertemplate=settings.home+"/madx/tracker_nominal_template3b.madx"
settings.pycollimate=False
settings.elements=['AP.UP.ZS21633','TPST.21760','AP.UP.MST21774']

settings.seed = 0

settings.nturns=204565
settings.nbatches=2
settings.nparperbatch=10
settings.ffile=5

submit_job(settings)

#And with faster sweep
settings.set_name('nominal_nopycoll_small')
settings.nturns=10000

submit_job(settings)

# And with pycollimate
settings.set_name('nominal_pycoll')
settings.pycollimate=True
settings.elements=['AP.UP.ZS21633_M','TPST.21760','AP.UP.MST21774']
settings.nturns=204565

submit_job(settings)

#And with faster sweep
settings.set_name('nominal_pycoll_small')
settings.nturns=10000

submit_job(settings)
