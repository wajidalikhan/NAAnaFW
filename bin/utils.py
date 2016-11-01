import os
import optparse 
import os.path
import optparse
import subprocess
import sys
import glob
import utils
import commands

def formSamples(proclist):
    samp=[]
    procs= proclist.split(",")
    for proc in procs:
        isAllProcesses = proc=="All"
        if proc=="ST" or isAllProcesses:
            from samplesST import samples as stemp
            samp.extend(stemp)
        if proc=="TT" or isAllProcesses:
            from samplesTT import samples as stemp
            samp.extend(stemp)
        if proc=="VJ" or isAllProcesses:
            from samplesVJ import samples as stemp
            samp.extend(stemp)
        if proc=="VV" or isAllProcesses:
            from samplesVV import samples as stemp
            samp.extend(stemp)
        if proc=="QCDMu" or isAllProcesses:
            from samplesQCDMu import samples as stemp
            samp.extend(stemp)
        if proc=="QCDEle" or isAllProcesses:
            from samplesQCDEle import samples as stemp
            samp.extend(stemp)
        if proc.startswith("_"):
            samp.append(proc[1:])
    return samp

def mymkdir(path):
    if commands.getstatusoutput("ls "+path)[0]!=0:
        os.system("mkdir " +path )
