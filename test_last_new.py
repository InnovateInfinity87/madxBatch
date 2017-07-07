# -*- coding: utf-8 -*-
"""
Example SPS slow extraction study: test.

@author: Linda Stoel
"""
from python.batching import Settings, submit_job


settings=Settings('endonly_pc', studygroup='test', disk='afsproject')

settings.saveout = False
settings.seed = 0

#settings.elements=['AP.UP.ZS21633','AP.DO.ZS21676','TPST.21760']
settings.elements=['AP.UP.ZS21633_M','AP.DO.ZS21676_M','TPST.21760']

settings.nturns=1000
settings.nbatches=2
settings.nparperbatch=2
settings.ffile=1
settings.trackendonly=30

settings.dynamicbump=False
settings.pycollimate=True

settings.flavour="tomorrow"

submit_job(settings)
