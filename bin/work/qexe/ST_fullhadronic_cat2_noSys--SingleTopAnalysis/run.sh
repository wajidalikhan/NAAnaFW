#!/bin/bash
START_TIME=`date`
cd /afs/cern.ch/work/w/wajid/NapoliFW/CMSSW_8_0_20/src/Analysis/NAAnaFW/bin/work/qexe/ST_fullhadronic_cat2_noSys--SingleTopAnalysis
source environment.sh
cd /afs/cern.ch/work/w/wajid/NapoliFW/CMSSW_8_0_20/src/Analysis/NAAnaFW/bin
CMD="ST work/files/ST.txt fullhadronic cat2 noSys noSync MC"
echo $CMD
ST work/files/ST.txt fullhadronic cat2 noSys noSync MC
res=$?
echo =========================================================
echo Exit code: $res
echo =========================================================
echo Done on `date` \| $res -  ST_fullhadronic_cat2_noSys--SingleTopAnalysis [started: $START_TIME ]>> qexe.log
