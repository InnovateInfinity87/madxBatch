# -*- coding: utf-8 -*-
"""
Example SPS slow extraction study: test.

@author: Linda Stoel
"""
from python.batching import Settings, submit_job

def track_test(k,data,settings):
    """Creates the text to replace pyTRACKER in tracker.madx. (nominal case)"""
    line = ("m_f = (kqf1_start - kqf1_end) / (1 - "+str(settings.nturns)+");\n"+
            "m_d = (kqd_start -  kqd_end ) / (1 - "+str(settings.nturns)+");\n\n"+

            "n_f = kqf1_start - m_f;\n"+
            "n_d = kqd_start - m_d;\n\n"+
        
            "SYSTEM, 'mkdir "+str(k)+"';\n\n")

    if settings.dynamicbump:
        line += (orthogonal_bumps()+"\n\n"+

                "x_knob_start = 406*dpp_start + (-0.0765*dpp_start*1e3);\n"+
                "px_knob_start = (-20.667*dpp_start) + (-0.0022*start*1e3);\n"+
                "EXEC, obump(x_knob_start, px_knob_start);\n"+
                "dpp_inc = (dpp_end-dpp_start)/"+str(settings.nturns)+";\n"
                "x_knob_inc = 406*dpp_inc + (-0.0765*dpp_inc*1e3);\n"+
                "px_knob_inc = (-20.667*dpp_inc) + (-0.0022*dpp_inc*1e3);\n\n")
        
    line += ("tr$macro(turn): MACRO = {\n"+
             " kqf1 = m_f * turn + n_f;\n"+
             " kqd = m_d * turn + n_d;\n")
    if settings.dynamicbump:
        line += (" EXEC, obump(x_knob_inc, px_knob_inc);\n")
    line += "};\n\n"

    line += ("OPTION, -WARN;\n"+
             "TRACK, ONEPASS, APERTURE, RECLOSS, UPDATE;\n\n")

    for i in range(k*settings.nparperbatch,(k+1)*settings.nparperbatch):
        line += (" START, X ="+str(data[0][i])+", "+
                       "PX = "+str(data[1][i])+", "+
                       "Y = "+str(data[2][i])+", "+
                       "PY = "+str(data[3][i])+", "+
                       "T = "+str(0)+", "+
                       "PT = "+str(data[4][i])+";\n")
    for place in settings.elements:
        line += " OBSERVE, PLACE = "+place+';\n';
    line += (" RUN, TURNS="+str(settings.nturns)+", "+
                 "MAXAPER={0.5,0.05,0.5,0.05,"+str(0.03*settings.nturns)+",0.1}, "+
                 "FFILE="+str(settings.ffile)+";\n\n"+

             "ENDTRACK;\n"+
             "OPTION, WARN;\n\n")

    line+="CALL, FILE='"+settings.home+"/madx/savetrack.cmdx';\n"

    for obsnum in range(len(settings.elements)+1):
        for partnum in range(settings.nparperbatch):
            temp = "obs"+str(obsnum+1).zfill(4)+".p"+str(partnum+1).zfill(4)
            line+="EXEC, savetrack(50, track."+temp+", '"+str(k)+"/"+temp+"');\n"
    line+="\n"

    if(settings.saveloss):
        line+="WRITE, TABLE = trackloss, FILE = 'losses.tfs';\n\n"

    line += "SYSTEM, 'tar -czf tracks.tar.gz "+str(k)+"';"

    return line




settings=Settings('pc_tracks', studygroup='test', disk='afsproject')

settings.trackerrep = track_test
settings.saveout = False
settings.seed = 0

settings.elements=['AP.UP.ZS21633_M','AP.DO.ZS21676_M','TPST.21760']

settings.nturns=50000
settings.nbatches=100
settings.nparperbatch=1000
settings.ffile=1

settings.dynamicbump=False
settings.pycollimate=True

settings.flavour="tomorrow"

submit_job(settings)
