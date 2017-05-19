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
        settings=Settings(name, studygroup='multitest', disk='afsprivate')

        #settings.local=True
        #settings.trackingbool=False

        settings.nturns=50000
        settings.nbatches=1000
        settings.nparperbatch=25
        settings.ffile=500
        settings.dppmax=0.0000

        replace['pyNTURNS']=str(settings.nturns)
        replace['pyDPPMAX']=str(settings.dppmax)
        if i==0:
            replace['pyBET1']='104.3784939';
            replace['pyBET2']='108.7895038';
            replace['pyBET3']='99.16930054';
        elif i==1:
            replace['pyBET1']='104.380386';
            replace['pyBET2']='107.1748819';
            replace['pyBET3']='99.8022297';
        else:
            replace['pyBET1']='103.6068183';
            replace['pyBET2']='106.5956281';
            replace['pyBET3']='101.2237158';

        settings.trackertemplate=settings.home+"/madx/tracker_multipole_template.madx"
        settings.myreplace=replace
        settings.pycollimate=False
        settings.elements=['AP.UP.ZS21633','AP.DO.ZS21676','AP.UP.TPST21760']

        settings.seed = seed

        submit_job(settings)
