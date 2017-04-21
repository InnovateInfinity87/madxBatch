# -*- coding: utf-8 -*-
"""
Example SPS slow extraction study: decapole folding

@author: Linda Stoel
"""
import math
from python.batching import Settings, submit_job

for nturns in [10000, 204565, 5000, 50000, 100000]:
    settings=Settings("sweepspeed_nominal_n"+str(nturns), studygroup="benchmarks", disk='afspublic')
    settings.pycollimate=False
    settings.elements=['AP.UP.ZS21633']

    settings.seed = 0

    settings.nturns=nturns
    settings.nbatches=10000
    settings.nparperbatch=10
    settings.ffile=math.floor(nturns/100)

    submit_job(settings)
