#!/bin/sh
#Running over muon and muonantiiso channel for MC
python new_singletop.py -c muon -s noSys --t3batch -f files/final/ -P ST,TT,VJ,VV,QCDMuEPt20toInf -S 10
python new_singletop.py -c muonantiiso -s noSys --t3batch -f files/final/ -P ST,TT,VJ,VV,QCDMuEPt20toInf -S 10

#Running if for Data and DDriven QCD
python new_singletop.py -c muon -s noSys --t3batch -f files/final/ -P SingleMuon -d DATA -S 10
python new_singletop.py -c muonantiiso -s noSys --t3batch -f files/final/ -P SingleMuon -d DATA -S 10

