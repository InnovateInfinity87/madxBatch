# -*- coding: utf-8 -*-
"""
Example SPS slow extraction study: downstream zs scan

@author: Linda Stoel
"""
from python.batching import Settings, submit_job

for zsdown in [4150, 4120, 4180]:
    settings=Settings(str(zsdown), studygroup='alignzs', disk='afspublic')

    settings.seed = 0

    settings.elements=['AP.UP.ZS21633_M','AP.UP.TPST21760']

    settings.septadb = settings.home+"/input/septa_DB_scan.tfs"
    settings.septadbreplace = {'zsdown': str(zsdown/100000.0)}

    settings.nturns=50000
    settings.nbatches=100
    settings.nparperbatch=100
    settings.ffile=1

    settings.dynamicbump=False
    settings.pycollimate=True

    settings.flavour='workday'

    submit_job(settings)
