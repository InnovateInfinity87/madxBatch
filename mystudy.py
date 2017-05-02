# -*- coding: utf-8 -*-
"""
Example SPS slow extraction study: decapole folding

@author: Linda Stoel
"""
from python.batching import Settings, submit_job

replacements = [{'pyK2': '-0.16', 'pyK4': '0.0'},
                {'pyK2': '-0.8', 'pyK4': '5222.6'},
                {'pyK2': '-1.05', 'pyK4': '7991.0'}]

for seed in range(1):
    for i, replace in enumerate(replacements):
        name = "case"+str(i)+"_seed"+str(seed)
        settings=Settings(name, studygroup='multitest', disk='afspublic')

        settings.trackertemplate=settings.home+"/madx/tracker_multipole_template.madx"
        settings.myreplace=replace
        settings.pycollimate=False
        settings.elements=['AP.UP.ZS21633','AP.DO.ZS21676','AP.UP.TPST21760']

        settings.seed = seed

        settings.nturns=50000
        settings.nbatches=100
        settings.nparperbatch=20
        settings.ffile=500

        submit_job(settings)
