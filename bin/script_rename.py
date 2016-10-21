import os, sys
import shutil
import errno

src='/afs/cern.ch/user/o/oiorio/public/xWajid/files/trees/mc/'
dest='/afs/cern.ch/work/w/wajid/NapoliFW/CMSSW_8_0_20/src/Analysis/NAAnaFW/bin/files/trees/mc'


def remove_folder(path):
    # check if folder exists
    if os.path.exists(path):
         # remove if exists
         print '*** Removing the folder : ',dest
         shutil.rmtree(path)
remove_folder(dest)
 
def copy(src, dest):
    try:
        shutil.copytree(src, dest)
    except OSError as e:
        # If the error was caused because the src wasn't a directory
        if e.errno == errno.ENOTDIR:
            shutil.copy(src, dest)
        else:
            print('Directory not copied. Error: %s' % e)

print '*** Copying the list of files from : ',src, 'to : ', dest 
copy(src,dest)
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
  
  if 'listST_s-channel_antitop_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1.txt' in i:
    print '*** Renaming file', i 
    os.renames(i,"ST_Tbar_sch.txt")
  
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
    os.renames(i,"WZTo1L1Nu2Q")
   
  if 'listWZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8.txt' in i:
    print '*** Renaming file', i 
    os.renames(i,"WZTo2L2Q")
    
  if 'listZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8.txt' in i:
    print '*** Renaming file', i 
    os.renames(i,"ZZTo2L2Q")

  if 'listQCD_Pt-20toInf_MuEnrichedPt15_TuneCUETP8M1_13TeV_pythia8.txt' in i:
    print '*** Renaming file', i 
    os.renames(i,"QCDMuEPt20toInf.txt")

print '*** Renaming Done ***'
