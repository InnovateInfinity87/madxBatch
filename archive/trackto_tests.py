import sys
from python.trackfurther import tracklossto

if __name__ == "__main__":
    folder = "/afs/cern.ch/project/sloex/zsalign/scandown_41550_sweep_nopc_nodb/test_losses"
    madxexe = "/afs/cern.ch/user/m/mad/bin/madx_dev"
#    tracklossto("tpstup", folder, "zsup", sloexcodedir="/afs/cern.ch/user/l/listoel/Desktop/SloExCode", madxexe=madxexe)
#    print "\n\n *** NEXT TRACKTO *** \n\n"
#    sys.stdout.flush()
#    tracklossto("handover", folder, "zsup", sloexcodedir="/afs/cern.ch/user/l/listoel/Desktop/SloExCode", madxexe=madxexe)
#    print "\n\n *** NEXT TRACKTO *** \n\n"
#    sys.stdout.flush()
#    tracklossto("QFA.21610", folder, "zsup", backtrack=True, sloexcodedir="/afs/cern.ch/user/l/listoel/Desktop/SloExCode", madxexe=madxexe)
#    print "\n\n *** NOW THE TWISSTRACKS *** \n\n"
    sys.stdout.flush()
    tracklossto("tpstup", folder, "zsup", sloexcodedir="/afs/cern.ch/user/l/listoel/Desktop/SloExCode", twisstrack=True, madxexe=madxexe)
    print "\n\n *** NEXT TRACKTO *** \n\n"
    sys.stdout.flush()
    tracklossto("handover", folder, "zsup", sloexcodedir="/afs/cern.ch/user/l/listoel/Desktop/SloExCode", twisstrack=True, madxexe=madxexe)
    print "\n\n *** NEXT TRACKTO *** \n\n"
    sys.stdout.flush()
    tracklossto("QFA.21610", folder, "zsup", backtrack=True, sloexcodedir="/afs/cern.ch/user/l/listoel/Desktop/SloExCode", twisstrack=True, madxexe=madxexe)
