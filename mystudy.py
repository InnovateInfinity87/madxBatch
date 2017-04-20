# -*- coding: utf-8 -*-
"""
Example SPS slow extraction study: decapole folding

@author: Linda Stoel
"""
from python.batching import Settings, submit_job

replacements = [{'pyK2': '0.16', 'pyK4': '0.0'},
                {'pyK2': '0.8', 'pyK4': '5222.6'},
                {'pyK2': '1.05', 'pyK4': '7991.0'}]

for seed in range(1):
    for i, replace in enumerate(replacements):
        name = "mystudy_case"+str(i)+"_seed"+str(seed)
        settings=Settings(name, disk='afsprivate')

        settings.trackertemplate=settings.home+"/madx/tracker_multipole_template.madx"
        settings.myreplace=replace
        settings.pycollimate=True
        settings.elements=['AP.UP.ZS21633_M','TPST.21760','AP.UP.MST21774']

        settings.seed = seed

        settings.nturns=10
        settings.nbatches=1
        settings.nparperbatch=10
        settings.ffile=1

        submit_job(settings)
