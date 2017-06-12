# -*- coding: utf-8 -*-
"""
Example SPS slow extraction study: downstream zs scan

@author: Linda Stoel
"""
from python.batching import Settings, submit_job

pos = [41550, 41050, 40550, 40050]
pc = False
db = True

for zsdown in pos:
    nameadd=""
    if pc:
        nameadd+="_pc"
    else:
        nameadd+="_nopc"
    if db:
        nameadd+="_db"
    else:
        nameadd+="_nodb"

    settings=Settings("do_"+str(zsdown)+"_nopc", studygroup='zsalign', disk='afsproject')

    settings.seed = 0
    settings.flavour = "espresso"
    settings.savetracks = False

    if pc:
        settings.pycollimate=True
        settings.septadb = settings.home+"/input/septa_DB_scan.tfs"
        settings.septadbreplace = {'zsdown': str(zsdown/1E6)}
        settings.elements=['AP.UP.ZS21633_M','AP.DO.ZS21676_M','TPST.21760']
    else:
        settings.pycollimate=False
        settings.elements=['AP.UP.ZS21633','AP.DO.ZS21676','AP.UP.TPST21760']

    settings.myreplace = {"zswirepos = 0.04245;": "zswirepos = "+str(zsdown/1.0E6)+";"}

    settings.nturns=300
    settings.nbatches=10
    settings.nparperbatch=1000
    settings.ffile=1

    if db:
        settings.dynamicbump=True
    else:
        settings.dynamicbump=False
    settings.slices=[-0.0015,-0.0010,-0.0005,0.0,0.0005,0.0010,0.0015]
    settings.slicewidth=0.0

    submit_job(settings)
