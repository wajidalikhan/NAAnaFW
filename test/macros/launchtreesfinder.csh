#cd $HOME
cd ../../../../../../CMSSW_8_1_0_pre12/src/
cmsenv
cd -
#python treesFinder.py -g gfal -F txt -f samples_HLTv1_muons.txt -o ../../bin/trees/mc/ -p /pnfs/lcg.cscs.ch/cms/trivcat/store/user/oiorio/ttDM/trees/2016/Oct/31Oct/ -d 161031 -V 2 > & mcfinder.log &
python treesFinder.py -g gfal -F txt -f samples_data_mu_edmndtuple.txt -D True -o ../../bin/trees/mc/ -p /pnfs/lcg.cscs.ch/cms/trivcat/store/user/oiorio/ttDM/trees/2016/Oct/12Oct/SingleMuon/ -V 2 > & datafinder.log &
#cmsenv
