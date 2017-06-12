#!/bin/sh
#rm -rf res/* trees/* work/* core.*

#Running over muon and muonantiiso channel for MC
#python new_singletop.py -c muon -t trees --mva mva -s noSys --t3batch -f files/final/ -P ST,TT,VJ,VV,QCDMu -S 4 
#python new_singletop.py -c muonantiiso -t trees --mva mva -s noSys --t3batch -f files/final/ -P ST,TT,VJ,VV,QCDMu -S 8
#nohup python new_singletop.py -c muon -t trees --mva mva -s noSys --t3batch -f files/final/ -P ST,TT,VJ,VV,QCDMu -S 5 &> mc.log &

nohup python new_singletop.py -c muon -t trees --mva mva -s noSys --t3batch -f files/final/ -P ST,TT,VJ,VV,QCDMu -S 3  -o /afs/cern.ch/work/w/wajid/public/xWajid &> muon_mc.log &
nohup python new_singletop.py -c muonantiiso -t trees --mva mva -s noSys --t3batch -f files/final/ -P ST,TT,VJ,VV,QCDMu -S 3  -o /afs/cern.ch/work/w/wajid/public/xWajid &> muonantisoi_mc.log &


#Running if for Data and DDriven QCD
#python new_singletop.py -c muon -t trees --mva mva -s noSys --t3batch -f files/final/ -P SingleMuon -d DATA -S 4
#python new_singletop.py -c muonantiiso -t trees --mva mva -s noSys --t3batch -f files/final/ -P SingleMuon -d DATA -S 8
#nohup python new_singletop.py -c muon -t trees --mva mva -s noSys --t3batch -f files/final/ -P SingleMuon -d DATA -S 3 &> data.log &
nohup python new_singletop.py -c muon -t trees --mva mva -s noSys --t3batch -f files/final/ -P SingleMuon -d DATA -S 3 -o /afs/cern.ch/work/w/wajid/public/xWajid &> muon_data.log &
nohup python new_singletop.py -c muonantiiso -t trees --mva mva -s noSys --t3batch -f files/final/ -P SingleMuon -d DATA -S 3 -o /afs/cern.ch/work/w/wajid/public/xWajid &> muonantiiso_data.log &

