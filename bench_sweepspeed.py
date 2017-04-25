# -*- coding: utf-8 -*-
"""
Example SPS slow extraction study: decapole folding

@author: Linda Stoel
"""
import math
from python.batching import Settings, submit_job

for nturns in [1000, 5000, 10000, 50000, 100000, 150000, 204565]:
    for db in [False, True]:
        settings=Settings("sweepspeed_db"+str(int(db))+"_n"+str(nturns), studygroup="benchmarks", disk='afspublic')
        settings.dynamicbump=db
        settings.pycollimate=False
        settings.elements=['AP.UP.ZS21633']
        settings.trackertemplate=settings.home+"/madx/tracker_nominal_template_v2.madx"

        settings.seed = 0

        settings.nturns=nturns
        settings.nbatches=10000
        settings.nparperbatch=10
        settings.ffile=math.floor(nturns/100)

        submit_job(settings)
