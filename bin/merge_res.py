import os
import optparse 
import os.path
import optparse
import subprocess
import sys
import glob
import commands
from utils import *
#import colorama
#from colorama import Fore, Back, Style
#text = "The quick brown fox jumps over the lazy dog"
#print(Fore.RED + text)

#More complex working example: this will run only the ST_T_tch and the V+Jets samples and split them into batches of 10 files, taking them from the remote folder on Orso's public:
#Usage: python merge_res.py -l ./res -P ST,TT,VJ,VV --rm True

from os.path import join,exists
print 'Python version', sys.version_info
print '--------------------'
if sys.version_info < (2, 7):
    raise "Must use python 2.7 or greater. Have you forgotten to do cmsenv?"

workdir = 'work'
fileListDir = join(workdir,'files')

usage = 'usage: python merge_res.py -l ./res -P ST,TT,VJ,VV --rm True'
parser = optparse.OptionParser(usage)
#Input files:
parser.add_option('-P', '--process',        dest='process',  type='string',     default = 'all', help="samples to add, according to the convention in 'script_rename.py'. Options are: 'All','ST','VJ','VV','QCDMu','QCDEle', or '_'+specific process, e.g. _ST_T_tch to add only the _ST_T_tch. Accepts also multiple processes, separated by comma, e.g. -P ST,_TT,VJ will select the V+jets samples and single top sample sets, as well as the one sample named TT.")

#Details of the analysis step:
parser.add_option('-n', '--dryrun',   dest='dryrun',  action='store_true', default=False)
#Splitting options:
parser.add_option('-l', '--localpath',        dest='localpath',  type='string',     default = 'res/', )
parser.add_option('--merge',   dest='domerge',  action='store_true', default=True)
parser.add_option('--rm',   dest='doremove',  action='store_true', default=False)
parser.add_option('--rmf',   dest='doforceremove',  action='store_true', default=False)

(opt, args) = parser.parse_args()



#define samples, one folder for each mass value
samples = []

samples = formSamples(opt.process)

print "samples are:",samples


path = opt.localpath
if not path.endswith("/"): path = path+"/"
channels=["muon","electron"]
channels=["muon"]
systs = ["noSyst"]
for s in samples:
    for c in channels:
        systs=[]
        cmdlsbase = "ls "+path+s+"*_part0*"+c+"*.root"
        
        lines = commands.getstatusoutput(cmdlsbase)[1].splitlines()
        for l in lines:
            systappend=(l.split("_")[-1])[:-5]
            systs.append(systappend) if systappend!=c else systs.append("noSyst")
            print "appended syst",systs[-1]
        for sys in systs:
            print "check syst",sys
            cmdls = "ls "+path+s+"*_part*"+c+"_"+sys+".root"
            if "noSys" in sys:
                cmdls = "ls "+path+s+"*_part*"+c+".root"

            if (commands.getstatusoutput(cmdls)[0]==0):
                linessys = commands.getstatusoutput(cmdls)[1].splitlines()
                print "looking for files with the command: ", cmdls
                #print linessys[0]
                halves= linessys[0].split("_part0")
                #print halves
                mergeoutput=False
                if(len(halves)==2 and opt.domerge):
                    cmdmerge="hadd -f "+halves[0]+halves[1]+" "+path+s+"*_part*"+c+"_"+sys+".root"
                    if "noSys" in sys:
                        cmdmerge="hadd -f "+halves[0]+halves[1]+" "+path+s+"*_part*"+c+".root"
                    print "merging with the command: ",cmdmerge
                    mergeoutput=True
                    mergeoutput=(commands.getstatusoutput(cmdmerge)[0]==0)
                    if (mergeoutput): print "merge successful!"
                    
                if (mergeoutput and opt.doremove):
                    cmdrm = "rm "+path+s+"*_part*"+c+"_"+sys+".root"
                    if "noSys" in sys:
                        cmdrm = "rm "+path+s+"*_part*"+c+".root"
                    print " now removing with the command: ", cmdrm
                    os.system(cmdrm)

                if opt.doforceremove :
                    cmdrm = "rm "+path+s+"*_part*"+c+"_"+sys+".root"
                    if "noSys" in sys:
                        cmdrm = "rm "+path+s+"*_part*"+c+".root"
                    print " forcing removal with with the command: ", cmdrm
                    os.system(cmdrm)
                    
                
                else: 
                    if len(halves)!=2: print "careful! didn't find part 0 for channel ",c," sys ",sys
