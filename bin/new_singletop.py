import os
import optparse 
import os.path
import optparse
import subprocess
import sys
import glob

#Usage: python new_singletop.py -c fullhadronic -s noSys --t3batch 
#Usage: python new_singletop.py -c fullhadronic -s noSys -m t3se 
#Usage: python new_singletop.py -c fullhadronic -s noSys -m local 

from os.path import join,exists
print 'Python version', sys.version_info
print '--------------------'
if sys.version_info < (2, 7):
    raise "Must use python 2.7 or greater. Have you forgotten to do cmsenv?"

workdir = 'work'
fileListDir = join(workdir,'files')

#define samples paths
#pathlocal = "/afs/cern.ch/work/w/wajid/NapoliFW/CMSSW_8_0_16/src/Analysis/NAAnaFW/test/crab_projects/crab_st_top/results/ST/" 
pathlocal = "/afs/cern.ch/work/w/wajid/NapoliFW/CMSSW_8_0_20/src/Analysis/NAAnaFW/bin/ST/"
filepath='/afs/cern.ch/work/w/wajid/NapoliFW/CMSSW_8_0_20/src/Analysis/NAAnaFW/bin/files/trees/'

#define samples, one folder for each mass value
samples = []
samples.append("ST")
#samples.append("STbar")

usage = 'usage: %prog -l lumi'
parser = optparse.OptionParser(usage)
parser.add_option('-c', '--channel',  dest='channel', type='string',     default = 'singleH', help='Channel to analyze: singleH or singleZ')
parser.add_option('-C', '--cat',      dest='cat',     type='string',     default = 'cat2',    help='Category to analyze: cat0 or cat1 or cat2')
parser.add_option('-s', '--sys',      dest='sys',     type='string',     default = 'noSys',   help='Systematics: jesUp, jesDown, jerUp, jerDown')
parser.add_option('',   '--sync',     dest='sync',    type='string',     default = 'noSync',  help='Synchro exercise')
parser.add_option('-d', '--isData',   dest='isData',  type='string',     default = 'MC',      help='is Data or MC?')
parser.add_option('-g', '--gdb',      dest='gdb',     action='store_true', default=False)
parser.add_option('-n', '--dryrun',   dest='dryrun',  action='store_true', default=False)
parser.add_option('-m', '--mode',     dest='mode',    default='t3se', choices=['local','t3se'])
parser.add_option('--t3batch',        dest='t3batch', action='store_true', default=False)

isData="MC"
(opt, args) = parser.parse_args()

if opt.sys not in ["noSys", "jesUp", "jesDown", "jerUp", "jerDown", "metUnclUp", "metUnclDown"]:
    parser.error('Please choose an allowed value for sys: "noSys", "jesUp", "jesDown", "jerUp", "jerDown","metUnclUp", "metUnclDown"')

# Create working area if it doesn't exist
if not exists(fileListDir):
    os.makedirs(fileListDir)

# Create test/ directory if it doesn't exist
if not exists('test'):
    os.makedirs('test')

for s in samples:
    if (s.startswith("JetHT") or s.startswith("SingleMu") or s.startswith("SingleEl") or  s.startswith("MET")): isData="DATA"
    
    if opt.mode == 'local':
        print 'Info: Running in local mode ...'
        sPath = join(pathlocal,'*.root')
        print 'Info: Looking for the *.root files at: ',sPath
        
        # Get the complete list of files
        # listing = subprocess.check_output(lLs.split()+[sPath])
        files = glob.glob(sPath)
        print 'Info: Sample',s,'Files found',len(files)

    elif opt.mode == 't3se':
        print 'Info: Running on t3se in interactive mode ...'
        sPath = join(pathlocal,'*.root')
        sT2Path = join(filepath,s+'.txt')
        print sT2Path 
        f = open(sT2Path,'r')
        listing = f.read()
        files = listing.split()
        print 'Info: Sample',s,'Files found',len(files) 
        f.close() 

    # Save it to a semi-temp file
    sampleFileList = join(fileListDir,s+'.txt')
    
    print 'Info :',sampleFileList
    with open(sampleFileList,'w') as sl:
        sl.write('\n'.join(files))
    
    outDirs = ['res','trees']
    
    for d in outDirs:
        if exists(d): continue
        os.makedirs(d)
    
    cmd = 'SingleTopAnalysis '+ s + ' ' + sampleFileList  + ' ' + opt.channel + ' ' + opt.cat + ' ' + opt.sys + ' ' + opt.sync + ' ' + isData
    print cmd

    if opt.gdb:
        cmd = 'gdb --args '+cmd
    
    elif opt.t3batch:
        #usage: qexe.py [-h] [-w WORKDIR] [-q QUEUE] [-n] [-e {qsub,bsub}]
        print 'Info: Running on t3se in batch mode ...'
        sPath = join(pathlocal,'*.root')
        que = '-q 8nh'
        dry = '-n'
        batch = '-e bsub'
        jid = '%s_%s_%s_%s' % (s,opt.channel,opt.cat,opt.sys)
        #cmd = './run.py -w ' + workdir + ' ' + que + ' ' + dry + ' ' + batch + ' ' + jid + ' '+cmd
        cmd = './run.py -w ' + workdir + ' ' + que + ' ' + batch + ' ' + jid + ' '+cmd
     
    print 'Info:',cmd
 
    if opt.dryrun:
        print 'Dry Run (command will not be executed)'
        continue
    
    print '--------------------'
    subprocess.call(cmd,shell=True)

