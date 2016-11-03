import os, sys
import shutil
import errno
os.seteuid(os.geteuid())
import optparse 
from utils import *
#usage: python script_rename.py -s ~oiorio/public/xWajid/files/trees/oct31/ -d files/renamed 

usage = 'usage: %prog '
parser = optparse.OptionParser(usage)
parser.add_option('-s', '--src',        dest='src',  type='string',     default = './', help="directory to rename files")
parser.add_option('-d', '--dest',       dest='dest',  type='string',    default = './renamed/', help="directory to rename files")

(opt, args) = parser.parse_args()
#src='./'
#dest='./renamed'

src = opt.src
dest= opt.dest
mymkdir(dest)
#remove_folder(dest)

print '*** Copying the list of files from : ',src, 'to : ', dest 
#copytree(src,dest)
copyfiles(src,dest)
os.chdir(dest)
print "Current directory is: %s" %os.getcwd()

# listing directories
list=os.listdir(os.getcwd())

for i in list: 
  if 'listST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1.txt' in i:
    print '*** Renaming file', i 
    os.renames(i,"ST_T_tch.txt")
    
  if 'listST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1.txt' in i:
    print '*** Renaming file', i 
    os.renames(i,"ST_Tbar_tch.txt")

  if 'listST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1.txt' in i:
    print '*** Renaming file', i 
    os.renames(i,"ST_T_sch.txt")
  
  if 'listST_tW_top_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1.txt' in i:
    print '*** Renaming file', i 
    os.renames(i,"ST_T_tW.txt")
  
  if 'listST_tW_antitop_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1.txt' in i:
    print '*** Renaming file', i 
    os.renames(i,"ST_Tbar_tW.txt")

  if 'listTT_TuneCUETP8M1_13TeV-powheg-pythia8.txt' in i:
    print '*** Renaming file', i 
    os.renames(i,"TT.txt")
  
  if 'listWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.txt' in i:
    print '*** Renaming file', i 
    os.renames(i,"WJetsToLNu.txt")
  
  if 'listDYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.txt' in i:
    print '*** Renaming file', i 
    os.renames(i,"DYJetsToLL.txt")

  if 'listWWTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8.txt' in i:
    print '*** Renaming file', i 
    os.renames(i,"WWTo1L1Nu2Q.txt")
  
  if 'listWWTo2L2Nu_13TeV-powheg.txt' in i:
    print '*** Renaming file', i 
    os.renames(i,"WWTo2L2Nu.txt")

  if 'listWZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8.txt' in i:
    print '*** Renaming file', i 
    os.renames(i,"WZTo1L1Nu2Q.txt")
   
  if 'listWZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8.txt' in i:
    print '*** Renaming file', i 
    os.renames(i,"WZTo2L2Q.txt")
    
  if 'listZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8.txt' in i:
    print '*** Renaming file', i 
    os.renames(i,"ZZTo2L2Q.txt")

  if 'listQCD_Pt-20toInf_MuEnrichedPt15_TuneCUETP8M1_13TeV_pythia8.txt' in i:
    print '*** Renaming file', i 
    os.renames(i,"QCDMuEPt20toInf.txt")

  if 'listSingleMuon' in i:
    print '*** Renaming file', i 
    os.renames(i,i.replace("listSingleMuon","SingleMuon"))

print '*** Renaming Done ***'
