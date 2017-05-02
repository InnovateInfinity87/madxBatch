# -*- coding: utf-8 -*-
"""
Nominal SPS slow extraction with and without scattering

@author: Linda Stoel
"""
from python.batching import Settings, submit_job

#Simulate without pycollimate, without dynamic bump
settings=Settings('nominal_nopc_nodb', studygroup="nominal", disk='afspublic')

settings.pycollimate=False
settings.elements=['AP.UP.ZS21633','AP.DO.ZS21676','TPST.21760']

settings.seed = 0

settings.nturns=204565
settings.nbatches=1000
settings.nparperbatch=100
settings.ffile=2045

#submit_job(settings)

#And with faster sweep
settings.set_name('nominal_nopc_nodb_fast')
settings.nturns=50000
settings.ffile=500

#submit_job(settings)

#With pycollimate, without dynamic bump
settings.set_name('nominal_pc_nodb')

settings.pycollimate=True
settings.elements=['AP.UP.ZS21633_M','AP.DO.ZS21676_M','TPST.21760']

settings.nturns=204565
settings.nbatches=1000
settings.nparperbatch=100
settings.ffile=2045

#submit_job(settings)

#And with faster sweep
settings.set_name('nominal_pc_nodb_fast')
settings.nturns=50000
settings.ffile=500

submit_job(settings)

#Without pycollimate, with dynamic bump
settings.set_name('nominal_nopc_db')

settings.dynamicbump=True
settings.pycollimate=False
settings.elements=['AP.UP.ZS21633','AP.DO.ZS21676','TPST.21760']

settings.seed = 0

settings.nturns=204565
settings.nbatches=1000
settings.nparperbatch=100
settings.ffile=2045

#submit_job(settings)

#And with faster sweep
settings.set_name('nominal_nopc_db_fast')
settings.nturns=50000
settings.ffile=500
settings.flavour="tomorrow"

#submit_job(settings)

#With pycollimate, with dynamic bump
settings.set_name('nominal_pc_db')

settings.pycollimate=True
settings.elements=['AP.UP.ZS21633_M','AP.DO.ZS21676_M','TPST.21760']

settings.nturns=204565
settings.nbatches=1000
settings.nparperbatch=100
settings.ffile=2045

#submit_job(settings)

#And with faster sweep
settings.set_name('nominal_pc_db_fast')
settings.nturns=50000
settings.ffile=500
settings.flavour="tomorrow"

submit_job(settings)
