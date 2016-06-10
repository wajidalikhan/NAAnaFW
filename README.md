# NAAnaFW
setenv SCRAM_ARCH slc6_amd64_gcc530
#cmsrel CMSSW_8_0_10_patch2
cd CMSSW_8_0_10_patch2/src
cmsenv
source /cvmfs/cms.cern.ch/crab3/crab.csh
cmsenv

mkdir Analysis/NAAnaFW
git clone https://github.com/oiorio/NAAnaFW.git Analysis/NAAnaFW