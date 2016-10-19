#!/bin/sh
grep -E 'root' < work/files/ST.txt | cut -d ':' -f2 | awk -v prefix="root://cms-xrd-global.cern.ch/" '{print prefix $0}' >  work/files/ST_1.txt
mv work/files/ST_1.txt work/files/ST.txt

grep -E 'root' < work/files/STbar.txt | cut -d ':' -f2 | awk -v prefix="root://cms-xrd-global.cern.ch/" '{print prefix $0}' >  work/files/STbar_1.txt
mv work/files/STbar_1.txt work/files/STbar.txt

#grep -E 'root' < work/files/TT.txt | cut -d ':' -f2 | awk -v prefix="root://cms-xrd-global.cern.ch/" '{print prefix $0}' >  work/files/TT_1.txt
#mv work/files/TT_1.txt work/files/TT.txt
