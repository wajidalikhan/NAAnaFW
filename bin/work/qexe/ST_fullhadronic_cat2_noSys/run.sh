#!/bin/bash
START_TIME=`date`
cd /afs/cern.ch/work/w/wajid/NapoliFW/CMSSW_8_0_20/src/Analysis/NAAnaFW/bin/work/qexe/ST_fullhadronic_cat2_noSys
source environment.sh
cd /afs/cern.ch/work/w/wajid/NapoliFW/CMSSW_8_0_20/src/Analysis/NAAnaFW/bin
CMD="SingleTopAnalysis ST work/files/ST.txt fullhadronic cat2 noSys noSync MC"
echo $CMD
SingleTopAnalysis ST work/files/ST.txt fullhadronic cat2 noSys noSync MC
res=$?
echo =========================================================
echo Exit code: $res
echo =========================================================
echo Done on `date` \| $res -  ST_fullhadronic_cat2_noSys [started: $START_TIME ]>> qexe.log
