#!/bin/sh
#python merge_res.py -l ./res -P SingleMuon -c "muon,muonantiiso"
#python merge_res.py -l ./res -P ST,TT,VJ,VV,QCDMu -c "muon,muonantiiso"

nohup python merge_res.py -l ./res -P SingleMuon,ST,TT,VJ,VV,QCDMu -c "muon" >& merge.log &

#python merge_res.py -l trees -P ST,TT,VJ,VV,QCDMu -o treesmerge -t 1
nohup python merge_res.py -l trees -P SingleMuon,ST,TT,VJ,VV,QCDMu -o treesmerge -t 1 >& trees.log &
