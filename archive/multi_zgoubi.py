# -*- coding: utf-8 -*-
"""
Example SPS slow extraction study: decapole folding

@author: Linda Stoel
"""
from python.batching import Settings, submit_job

replacements = [{'pyK2': '-0.37', 'pyK4': '0.0', 'pyEXTRBUMP': None},
                {'pyK2': '-0.73', 'pyK4': '5018.3', 'pyEXTRBUMP': None},
                {'pyK2': '-1.02', 'pyK4': '7640.7', 'pyEXTRBUMP': 0.032}]

for seed in range(1):
    for i, replace in enumerate(replacements):
        if replace['pyEXTRBUMP']==None:
            replace['pyEXTRBUMP'] = 'knob_extr_bump'
        else:
            replace['pyEXTRBUMP'] = str((0.068-replace['pyEXTRBUMP'])/0.0436825)+"*knob_extr_bump"

        

        name = "case"+str(i)+"_seed"+str(seed)
        settings=Settings(name, studygroup='multi_zgoubi', disk='afsprivate')

        settings.local=True
        settings.trackingbool=False

        settings.nturns=50000
        settings.nbatches=100
        settings.nparperbatch=25
        settings.ffile=500
        settings.dppmax=0.0000

        replace['pyNTURNS']=str(settings.nturns)
        replace['pyDPPMAX']=str(settings.dppmax)
        if i==0:
            replace['pyBET1']='103.2759543';
            replace['pyBET2']='105.3714172';
            replace['pyBET3']='95.23825079';
        elif i==1:
            replace['pyBET1']='103.2370701';
            replace['pyBET2']='103.8088785';
            replace['pyBET3']='95.85321099';
        else:
            replace['pyBET1']='102.4682338';
            replace['pyBET2']='103.2312348';
            replace['pyBET3']='97.21188502';

        settings.trackertemplate=settings.home+"/madx/tracker_multipole_zgoubi_template.madx"
        settings.myreplace=replace
        settings.pycollimate=False
        settings.elements=['AP.UP.ZS21633','AP.DO.ZS21676','AP.UP.TPST21760']

        settings.seed = seed

        submit_job(settings)
