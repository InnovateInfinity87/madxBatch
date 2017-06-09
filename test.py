# -*- coding: utf-8 -*-
"""
Example SPS slow extraction study: decapole folding

@author: Linda Stoel
"""
from python.batching import Settings, submit_job

settings=Settings('thickslice_bump_alt', studygroup='test', disk='afsproject')

settings.trackingbool=True
settings.trackertemplate=settings.home+"/madx/tracker_nominal_template.madx"
settings.local=False
settings.monitor=False
settings.seed = 0

settings.elements=['AP.UP.ZS21633','AP.DO.ZS21676','AP.UP.TPST21760']
#settings.elements=['AP.UP.ZS21633_M','AP.DO.ZS21676_M','TPST.21760']

settings.nturns=300
settings.nbatches=10
settings.nparperbatch=100
settings.ffile=1

settings.slices=[-0.0015,0.0,0.0015]
#settings.slicewidth=0.0
settings.dynamicbump=True
settings.pycollimate=False

settings.flavour=None

settings.saveout = False

submit_job(settings)
