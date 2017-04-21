# -*- coding: utf-8 -*-
"""
Fake batching by local execution of an HTCondor submit file

Functions to read in a simple HTCondor submit file and execute the tasks
sequentially on the local machine instead, in order to test batched code
without the need to run on the cluster.

Not guaranteed to work well at all for most cases, just for my simple ones. E.g.:
* Only a single "queue" command per submit file is supported for the moment,
and only with it's most basic use cases: "queue" and "queue #".
* The only variables handled are $(ClusterId) and $(ProcId)

@author: Linda Stoel
"""
import os
import shutil
import subprocess

def readsub(filename, here):
    subfile = {}
    with open(filename) as f:
        for line in f:
            if "=" in line:
                command = line.split('=', 1)
                subfile[command[0].strip()] = command[1].strip()
            elif "queue" in line:
                command = line.split()
                if len(command)==1:
                    subfile["queue"] = 1
                else:
                    subfile["queue"] = int(command[1])
    if "initialdir" not in subfile:
        subfile["initialdir"] = here

    return subfile


def runsub(subfile):
    if os.path.exists("fakeremote"):
        print("ERROR: Directory fakeremote already exists, please delete before trying again.")
        return

    queue = subfile["queue"]
    del subfile["queue"]
    if "requirements" in subfile:
        del subfile["requirements"]
    for key in subfile.keys():
        if key.startswith("+"):
            del subfile[key]
        else:
            subfile[key] = subfile[key].replace("$(ClusterId)","fake")
            if key not in ["queue", "initialdir", "input", "output", "error",
                           "log", "transfer_input_files", "transfer_output_files",
                           "transfer_output_remaps", "executable", "arguments"]:
                print("WARNING: keyword "+key+" is currently not supported in local mode.")

    for i in range(queue):
        #Make job-dependent version of the submit file
        locsub = subfile.copy()
        for key in locsub:
            locsub[key] = locsub[key].replace("$(ProcId)",str(i))
        if not locsub["initialdir"].endswith("/"):
            locsub["initialdir"] = locsub["initialdir"]+"/"
        for key in ["input", "output", "error", "log"]:
            if key in locsub:
                if not locsub[key].startswith("/"):
                    locsub[key] = locsub["initialdir"]+locsub[key]
        for key in ["transfer_input_files"]:
            if key in locsub:
                temp = ""
                for fname in [x.strip() for x in locsub[key].split(",")]:
                    if not fname.startswith("/"):
                        temp = temp+locsub["initialdir"]+fname+", "
                    else:
                        temp = temp+fname.strip()+", "
                locsub[key] = temp[:-2]
        if "transfer_output_remaps" in locsub:
            temp = ""
            for fnames in locsub["transfer_output_remaps"][1:-1].split(";"):
                remote,local = fnames.split("=")
                remote = remote.strip()
                local = local.strip()
                if not local.startswith("/"):
                    local = locsub["initialdir"]+local
                temp = temp+remote+"="+local+'; '
            locsub["transfer_output_remaps"] = temp[:-2]

        #Create a fake remote and start logging
        os.mkdir("fakeremote")
        log = open(locsub["log"],"a")
        log.write("submitting job "+str(i)+" to fake remote\n")

        #Copy needed input to fake remote
        shutil.copy(locsub["executable"],"fakeremote/executable")
        if "input" in locsub:
            shutil.copy(locsub["input"],"fakeremote/input")
        if "transfer_input_files" in locsub:
            for fname in locsub["transfer_input_files"].split(","):
                shutil.copy(fname.strip(),"fakeremote/"+fname.split("/")[-1].strip())

        #Go to the fake remote and execute the job
        os.chdir("fakeremote")
        with open("out.out", 'w') as out, open("err.err", 'w') as err:
            subprocess.check_call("chmod 755 ./executable", shell=True)
            if "arguments" in locsub:
                arguments = locsub["arguments"][1:-1]
            else:
                arguments = ""
            if "input" in locsub:
                try:
                    subprocess.check_call("./executable "+arguments+" < input", stdout=out, shell=True)
                except subprocess.CalledProcessError as e:
                    err.write(str(e.output))
                    err.write("\n")
                    return
            else:
                try:
                    subprocess.check_call("./executable "+arguments, stdout=out, shell=True)
                except subprocess.CalledProcessError as e:
                    err.write(str(e.output))
                    err.write("\n")
                    return

        #Exit the fake remote, copy output back to 'local' and destroy the remote
        os.chdir("..")
        if "output" in locsub:
            shutil.copyfile("fakeremote/out.out",locsub["output"])
        if "transfer_output_files" in locsub:
            for fname in locsub["transfer_output_files"].split(","):
                if os.path.isfile(fname):
                    try:
                        shutil.copyfile("fakeremote/"+fname.strip(), fname.strip().split("/")[-1])
                    except IOError as e:
                        log.write("Error occured in job "+str(i)+", copying "+fname.strip()+"\n")
                        log.write("I/O error({0}): {1}".format(e.errno, e.strerror))
                        log.write("\n")
                else:
                    try:
                        shutil.copytree("fakeremote/"+fname.strip(), fname.strip().split("/")[-1])
                    except IOError as e:
                        log.write("Error occured in job "+str(i)+", copying "+fname.strip()+"\n")
                        log.write("I/O error({0}): {1}".format(e.errno, e.strerror))
                        log.write("\n")
        if "transfer_output_remaps" in locsub:
            for fnames in locsub["transfer_output_remaps"].split(";"):
                remote,local = fnames.split("=")
                remote = remote.strip()
                local = local.strip()
                try:
                    shutil.copyfile("fakeremote/"+remote, local)
                except IOError as e:
                    log.write("Error occured in job "+str(i)+", copying "+remote+" to "+local+"\n")
                    log.write("I/O error({0}): {1}".format(e.errno, e.strerror))
                    log.write("\n")
        if "error" in locsub:
            shutil.copyfile("fakeremote/err.err",locsub["error"])
        shutil.rmtree("fakeremote")
        log.close()


def localsub(filename):
    runsub(readsub(filename, os.getcwd()))
