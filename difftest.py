# -*- coding: utf-8 -*-
"""
Example SPS slow extraction study: diffuser positioning

@author: Linda Stoel
"""
from python.batching import Settings, submit_job

for diffpos in [6900, 6895, 6890, 6905]:
    settings=Settings(str(diffpos), studygroup='difftest', disk='afspublic')

    settings.trackingbool=True
    settings.trackertemplate=settings.home+"/madx/tracker_diffuser_template.madx"
    settings.local=False
    settings.monitor=False
    settings.seed = 0

    settings.elements=['AP.UP.ZS21633_M','AP.UP.TPST21760']

    settings.septadb = settings.home+"/input/septa_DB_template.tfs"
    settings.septadbreplace = {'diffpos': str(diffpos/100000.0)}

    settings.nturns=10000
    settings.nbatches=1000
    settings.nparperbatch=25
    settings.ffile=1

    settings.dynamicbump=True
    settings.pycollimate=True

    settings.flavour='workday'

    submit_job(settings)
