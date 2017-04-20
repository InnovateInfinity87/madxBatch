# madxBatch

Conventions:
* The "python" command on your lxplus should be python 2.7
* An example script is given in mystudy.py the important points are:
    * from python.batching import Settings, submit_job
    * create a Settings object, by calling Settings with your project name as an argument. Adjust the settings from default as needed.
    * submit your study to the cluster with submit_job(settings)
* By default output is written to /eos/USER_INITIAL/USERNAME/madxBatchData/PROJECTNAME, but a different disk can be chosen when creating the Settings object with the optional argument disk=...:
    * 'afsprivate' : /afs/cern.ch/work/USER_INITIAL/USERNAME/private/madxBatchData/PROJECTNAME
    * 'afspublic' : /afs/cern.ch/work/USER_INITIAL/USERNAME/public/madxBatchData/PROJECTNAME
* Aditionally you can override this directory by using settings.datadir=DESIRED_OUTPUT_DIR in your script after the creation of the settings object, but be carefull that the directory will no longer automatically change with set_name if you do that!
* By default pycollimate is on, it can be turned off by setting settings.pycollimate=False in your script
* By default the madxBatch folder, the pycollimate folder and the slowExtractionMADX folder are assumed to be located in the same folder. This can be changed by setting settings.pycolldir=PYCOLL_DIR and settings.slowexdir=SLOWEX_DIR in your script.

Other reminders:
* If pycollimate is turned off, madx_dev is used. This means that in case of a large MAD-X update the code could be broken until it is updated.
