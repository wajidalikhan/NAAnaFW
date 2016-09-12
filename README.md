Instructions to use the Framework. Assuming the user has CERN computing account.
--------------------------------------------------------------------------------
Setting up the computing enviroment: The framework is tested on CMSSW_8_0_16. 

--  ssh -Y username@lxplus.cern.ch
--  export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
--  source $VO_CMS_SW_DIR/cmsset_default.sh  
--  export SCRAM_ARCH=slc6_amd64_gcc530  
--  mkdir NapoliFW
--  cd NapoliFW
--  cmsrel CMSSW_8_0_16 
--  cd src
--  cmsenv

Or copy and past the following lines to your ~/.bashrc, Just type NAanaFW and press enter:

NAanaFW(){
    export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
    source $VO_CMS_SW_DIR/cmsset_default.sh
    export SCRAM_ARCH=slc6_amd64_gcc530
    echo "VO_CMS_SW_DIR :" $VO_CMS_SW_DIR
    echo "Setting the SCRAM_ARCH ..."
    echo "SCRAM_ARCH :" $SCRAM_ARCH
    cd /afs/cern.ch/work/w/username/NapoliFW/CMSSW_8_0_16/src
    cmsenv
    echo "CMSSW BASE VERSION: "$CMSSW_VERSION 
    echo "PWD: "$PWD
    }

Step 1: Check out the local working copy of the framework:
--  export CMSSW_GIT_REFERENCE=/cvmfs/cms.cern.ch/cmssw.git.daily
--  git cms-init
--  git clone git@github.com:cmsb2g/B2GAnaFW.git Analysis/B2GAnaFW -b CMSSW_8_0_X_V2
--  git clone git@github.com:cms-jet/JetToolbox.git JMEAnalysis/JetToolbox -b jetToolbox_763
--  scram b -j 10

Running the framework:
The python configuration file for cmsRun is B2GAnaFW/test/b2gedmntuples_cfg.py. It runs on the 
miniAOD data tier and produces an EDM-ntuple. 

--  cd Analysis/B2GAnaFW/test
--  cmsRun b2gedmntuples_cfg.py maxEvents=1000 DataProcessing='MC_MiniAODv2_80X_reHLT'

Step 2:  
--  cd /afs/cern.ch/work/w/username/NapoliFW/CMSSW_8_0_16/src
--  mkdir Analysis/NAAnaFW
--  git clone https://github.com/oiorio/NAAnaFW.git Analysis/NAAnaFW
--  scram b -j 10
--  cd Analysis/NAAnaFW/test
--  cmsRun topplusdmTrees_cfg.py maxEvts=10 sample="file:/afs/cern.ch/work/w/username/NapoliFW/CMSSW_8_0_16/src/Analysis/B2GAnaFW/test/B2GEDMNtuple.root" outputLabel='analysisTTDM.root'

For the CRAB submission
--  python submit_all.py -c b2gedmntuples_cfg.py -f <your_dataset_file> -s "T3_US_FNALLPC" -p "DataProcessing=MC_MiniAODv2_80X" -d B2GEDMNTuples_80x_V1p0 -o "/store/group/lpctlbsm/B2GAnaFW_80X_V1p0" -v RunIISpring16MiniAODv2_B2GAnaFW_80x_V1p0 -i 'Summer15_25nsV7*'

Step 3: 
Next steps in progress ...

 
