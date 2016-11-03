import os
import optparse 
import os.path
import optparse
import subprocess
import sys
import glob
from utils import *
import commands
#Usage:python script_replacexrd.py -f ./files/renamed/ -o files/final -x xrootd.ba.infn.it -P ST,VV,VJ,TT,SingleMuon 

usage = ''
parser = optparse.OptionParser(usage)

parser.add_option('-f', '--filepath',        dest='filepath',  type='string',     default = './files/trees/mc/')
parser.add_option('-o', '--outputpath',        dest='outputpath',  type='string',     default = './files/trees/mc/newxrd/')
parser.add_option('-x', '--xrd',        dest='xrd',  type='string',     default = 'xrootd.ba.infn.it', help="new xrd redirector" )
parser.add_option('-P', '--process',        dest='process',  type='string',     default = 'all', help="samples to add, according to the convention in 'script_rename.py'. Options are: 'All','ST','VJ','VV','QCDMu','QCDEle', or '_'+specific process, e.g. _ST_T_tch to add only the _ST_T_tch. Accepts also multiple processes, separated by comma, e.g. -P ST,_TT,VJ will select the V+jets samples and single top sample sets, as well as the one sample named TT.")

(opt, args) = parser.parse_args()

samples = []
samples = formSamples(opt.process)

print "changing"

mymkdir(opt.outputpath)

for s in samples:
    fpath=opt.filepath+"/"+s+".txt"
    foutpath=opt.outputpath+"/"+s+".txt"
    f = open(fpath,'r')
    fout = open(foutpath,'w')
    listing = f.read()
    files = listing.split()
    talktome = True
    for l in files:
        xrdorig=l.split("//")[1]
        fout.write(l.replace(xrdorig,opt.xrd)+"\n")
        if talktome: 
            print 'substituting: ',xrdorig,' with ', opt.xrd, " for sample ",s
            talktome=False
    f.close() 
    fout.close()
