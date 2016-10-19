############################################################
######## part 1: Environment setup and compilation  ########

setenv SCRAM_ARCH slc6_amd64_gcc530

cd CMSSW_8_0_20/src

cmsenv

source /cvmfs/cms.cern.ch/crab3/crab.csh

cmsenv

mkdir Analysis/NAAnaFW

git clone https://github.com/oiorio/NAAnaFW.git Analysis/NAAnaFW

scram b -j 10 > & compilation.log &

#######################################
######## Part 2: Tree running  ########

cd test 

##MC:
nohup cmsRun topplusdmTrees_skim_cfg.py
 #runs a skimmed sample, just the nominal trees (no jes/jer), preselection of >=1 tight lepton, >=1 pfLoose jet with pt >40 GeV

##Data:
nohup cmsRun topplusdmTrees_skim_cfg.py isData=True

##Finding files for trees:
cd macros

python treesFinder.py -f samples_HLTv1_muons.txt -s storage01.lcg.cscs.ch:8443/srm/managerv2 -o samples/mc/ -p /pnfs/lcg.cscs.ch/cms/trivcat/store/user/oiorio/ ttDM/samples/2016/Oct/ -V 2 
 # This one fetches recursively all directories in a path [option -p] in the storage element [option -s] looking for samples indicated in 
 # a txt file [option -f], with some parsing options (veto  a word, look for a particular date etc. ). 
 # The output are python files including the list of samples, as list of strings with the xrootd paths to all.
 # Syntax:
  # python treesFinder.py [-f crabFilesToFetch] [-s storageElement] [-p basePathOnTheSE] [-o outputDirToStoreFiles] [-d dateInTheYYMMDD_HHMMSSCrabFormat] [-v VetoWord] [-V verboseLevelFrom0To2]
 # Can also use the built-in help function:
  #python treesFinder.py --help 

######################################################################
######## part 3: Running the detailed selection and analysis  ########


################################################
######## part 4: Statistical inference  ########
