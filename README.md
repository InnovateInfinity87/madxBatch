# madxBatch

Conventions:
- The "python" command on your machine should be python 2.7
- By default output is written to /afs/cern.ch/work/USER_INITIAL/USERNAME/private/madxBatchData/PROJECTNAME This can be changed by setting c.outputdir=DESIRED_OUTPUT_DIR in your script.
- By default pycollimate is on, it can be turned off by setting c.pycollimate=False in your script
- By default the madxBatch folder and the pycollimate folder are assumed to be located in the same folder. This can be changed by setting c.pycoldir=PYCOLL_DIR in your script.

Other reminders:
- If pycollimate is turned off, madx_dev is used. This means that in case of a large MAD-X update the code could be broken until it is updated.
