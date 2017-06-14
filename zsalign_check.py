# -*- coding: utf-8 -*-
"""
Example SPS slow extraction study: downstream zs scan

@author: Linda Stoel
"""
from python.batching import Settings, submit_job

#pos = [41250, 41300, 41350, 41400]
#pos = [41000, 41050, 41100]
pos = [40950, 41150]
pc = False
db = True

for zsdown in pos:
    name = "scandown_"+str(zsdown)+"_sweep"
    name += "_pc" if pc else "_nopc"
    name += "_db" if db else "_nodb"

    settings=Settings(name, studygroup='zsalign', disk='afsproject')

    settings.seed = 0
    settings.savetracks = False

    if pc:
        settings.pycollimate=True
        settings.septadb = settings.home+"/input/septa_DB_scan.tfs"
        settings.septadbreplace = {'zsdown': str(zsdown/1E6)}
        settings.elements=['AP.UP.ZS21633_M','AP.DO.ZS21676_M','TPST.21760']
        settings.nbatches = 100
        settings.nparperbatch = 1000
    else:
        settings.pycollimate=False
        settings.elements=['AP.UP.ZS21633','AP.DO.ZS21676','AP.UP.TPST21760']
        settings.nbatches = 500
        settings.nparperbatch = 200

    settings.myreplace = {"zswiredo = 0.04245;": "zswiredo = "+str(zsdown/1.0E6)+";"}

    settings.nturns=50000
    settings.ffile=50000

    if db:
        settings.dynamicbump=True
    else:
        settings.dynamicbump=False

    if pc and db:
        settings.flavour = "tomorrow"

    submit_job(settings)
