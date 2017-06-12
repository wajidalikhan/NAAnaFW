import os
import commands
import os.path
import optparse
import subprocess
import sys
import glob
import commands
import optparse
from subprocess import call


os.system("rm -r listjobs.txt")
os.system("bjobs -rp -o jobid | grep -v ^JOBID | cut -d ':' -f2 > listjobs.txt")
alljobs =open("listjobs.txt","r")
allines= alljobs.readlines()
resubmitlines=[]
for l in allines:
    cmdlsbase= "bpeek "+ l.split("\n")[0]
    print cmdlsbase
    linesout = commands.getstatusoutput(cmdlsbase)[1].splitlines()
    #linesout = ['hello','Error in <TNetXNGFile::Open>: [ERROR] Operation expired']
#    print linesout
    haserror=False
    for lerr in linesout:
        #print lerr
        #if "output from stderr" in lerr:
        if haserror: continue
        if ( ("[FATAL]" in lerr) and ("Redirect" in lerr) ) or ( ("[ERROR]" in lerr) and ("Invalid" in lerr)):# or ("No servers" in lerr)) or ("Operation expired")):
            print "job",l," has fatal redirect error: ",lerr
            resubmitlines.append(l)
            haserror=True
        elif (("[ERROR]" in lerr)) and (("Operation expired" in lerr) or ("No servers are available to read the file" in lerr)):
            print "Job", l, "has Operation expired or No server are available : ", lerr
            call("./fancy.sh") 
            resubmitlines.append(l)
            haserror=True
for r in resubmitlines:
    print "Job with error is",r,"Resubmitting ..."
    os.system("brequeue "+r.split("\n")[0])

