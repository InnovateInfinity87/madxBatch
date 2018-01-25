"""
Simulation of amplitude based extraction

Will work for pc/nopc and db/nodb, but only with sweep, not with slices.
Is a bit more hacky than momentum based...
"""
import sys
import numpy as np
sys.path.insert(1, '/afs/cern.ch/project/sloex/code/madxBatch')
from python.batching import Settings, submit_job
from python.batching import __version__ as currversion

testversion = '0.02.003'
if not currversion==testversion:
    print ('This example was tested with version '+testversion+
           ' and is not guaranteed to work properly in the current '+
           'version (v'+currversion+'). Please check.')


settings=Settings('ampex', studygroup='test', disk='afsproject')

settings.madxversion = '/afs/cern.ch/user/m/mad/bin/rel/5.03.06/madx-linux64-gnu'

settings.finalchanges = "/afs/cern.ch/user/l/listoel/Desktop/SloExCode/madxBatch/madx/ampex_finalchanges.cmdx"

settings.dynamicbump=True
settings.pycollimate=False
settings.seed = 0

settings.elements=['AP.UP.ZS21633','AP.DO.ZS21676','AP.UP.TPST21760']

settings.nturns=50000
settings.nbatches = 100
settings.nparperbatch = 50
settings.ffile=50000

if settings.dynamicbump:
    settings.flavour="workday"

    settings.dynamicbump_cx = -2.0
    settings.dynamicbump_cpx = 70.0

    settings.myreplace['pyDB'] = '1'
    settings.myreplace['pyCX'] = str(settings.dynamicbump_cx)
    settings.myreplace['pyCPX'] = str(settings.dynamicbump_cpx)
else:
    settings.flavour="longlunch"

submit_job(settings)
