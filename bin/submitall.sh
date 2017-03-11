#!/bin/sh
python new_singletop.py -C CMVAT -c muon -s noSys --t3batch -f files/final/ -P ST,TT,VJ,VV,QCDMu -S 9 
python new_singletop.py -C CMVAT -c muonantiiso -s noSys --t3batch -f files/final/ -P ST,TT,VJ,VV,QCDMu -S 9

python new_singletop.py -C CMVAT -c muon -s noSys --t3batch -f files/final/ -P SingleMuon -d DATA -S 9
python new_singletop.py -C CMVAT -c muonantiiso -s noSys --t3batch -f files/final/ -P SingleMuon -d DATA -S 9

