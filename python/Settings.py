# -*- coding: utf-8 -*-
"""
Settings for the batched simulation

@author: Wouter van de Pontseele, Linda Stoel
"""

import os
import sys

class Settings:
    """Settings for the batched simulations"""
    def __init__(self, name):
        home = sys.path[0]
        user = os.environ["USER"]
        outputdir = '/afs/cern.ch/work/'+user[0]+'/'+user+'/private/madxBatchData/'

        self.name=name
        self.datadir=outputdir+name+"/"
        self.pycolldir = home+'/../pycollimate/'
        self.slowexdir = home+'/../slowExtractionMADX/'

        self.trackertemplate = home+"/madx/tracker_template.madx"

        self.Nturns=34095
        self.Nbatches=1000
        self.Nparperbatch=50
        self.turnmultiplicity=50 #ffile MADX

        self.elements=['AP.UP.ZS21633','ZS.21633','AP.DO.ZS21633']

        self.createTwiss=False
        self.trackingBool=True
        self.writetrack=True
        self.pycollimate=True

        self.monitorcondor=False
