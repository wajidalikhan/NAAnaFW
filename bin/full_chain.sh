#!/bin/sh
echo "Renaming the files ..."
#python script_rename.py -s ~oiorio/public/xWajid/files/trees/nov09/ -d files/renamed
python script_rename.py -s ~oiorio/public/xWajid/files/trees/dec19/ -d files/final

echo "Renaming the xrd ..."
python script_replacexrd.py -f ./files/renamed/ -o files/final -x xrootd.ba.infn.it -P ST,VV,VJ,TT,SingleMuon

echo "Submitting jobs to cluster ... Remember to merge the output with the script:" 
python merge_res.py -l ./res/ -P ST,TT,VJ,VV,QCDMuEPt20toInf --rm True

python merge_res.py -l ./res -P SingleMuon -c "muon,muonantiiso"
python merge_res.py -l ./res/ -P ST,TT,VJ,VV,QCDMuEPt20toInf -c "muon,muonantiiso"

python merge_res.py -l ./res/ -P _QCDMuEPt20toInf -c "muon,muonantiiso"

nohup python merge_res.py -l /afs/cern.ch/work/w/wajid/public/xWajid/res -P SingleMuon,ST,TT,VJ,VV,QCDMu -c "muon" -o /afs/cern.ch/work/w/wajid/public/xWajid/mergeres >& mergemuon.log &

