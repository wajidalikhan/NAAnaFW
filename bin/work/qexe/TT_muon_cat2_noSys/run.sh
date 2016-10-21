#!/bin/bash
START_TIME=`date`
cd /afs/cern.ch/work/w/wajid/NapoliFW/CMSSW_8_0_20/src/Analysis/NAAnaFW/bin/work/qexe/TT_muon_cat2_noSys
source environment.sh
cd /afs/cern.ch/work/w/wajid/NapoliFW/CMSSW_8_0_20/src/Analysis/NAAnaFW/bin
CMD="SingleTopAnalysis TT work/files/TT.txt muon cat2 noSys noSync MC"
echo $CMD
SingleTopAnalysis TT work/files/TT.txt muon cat2 noSys noSync MC
res=$?
echo =========================================================
echo Exit code: $res
echo =========================================================
echo Done on `date` \| $res -  TT_muon_cat2_noSys [started: $START_TIME ]>> qexe.log
