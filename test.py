# -*- coding: utf-8 -*-
"""
Example SPS slow extraction study: decapole folding

@author: Linda Stoel
"""
from python.batching import Settings, submit_job

settings=Settings('twiss', studygroup='multitest', disk='afspublic')

settings.trackingbool=False
settings.trackertemplate=settings.home+"/madx/tracker_multipole_template.madx"
settings.local=True
settings.monitor=False
settings.seed = 0

settings.elements=['AP.UP.ZS21633','AP.DO.ZS21676','AP.UP.TPST21760']

settings.nturns=10000
settings.nbatches=1
settings.nparperbatch=100
settings.ffile=1

settings.dynamicbump=False
settings.pycollimate=False

settings.flavour=None

submit_job(settings)
