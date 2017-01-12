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


#Running over muon and muonantiiso channel for MC
python new_singletop.py -c muon -s noSys --t3batch -f files/final/ -P ST,TT,VJ,VV,QCDMuEPt20toInf -S 10
python new_singletop.py -c muonantiiso -s noSys --t3batch -f files/final/ -P ST,TT,VJ,VV,QCDMuEPt20toInf -S 10

#Running if for Data and DDriven QCD
python new_singletop.py -c muon -s noSys --t3batch -f files/final/ -P SingleMuon -d DATA -S 10
python new_singletop.py -c muonantiiso -s noSys --t3batch -f files/final/ -P SingleMuon -d DATA -S 10



#Running over the MC 
#python new_singletop.py --t3batch -f files/final/ -P ST,TT,VJ,VV -S 10

#Running over muonantiiso
#python new_singletop.py -c muonantiiso -s noSys --t3batch -f files/final/ -P ST,TT,VJ,VV,QCDMuEPt20toInf -S 10

#Running over the data
#python new_singletop.py --t3batch -f files/final/ -P SingleMuon -d DATA -S 10

#Running if for Data Driven QCD
#python new_singletop.py -c muonantiiso -s noSys --t3batch -f files/final/ -P SingleMuon -d DATA -S 10
