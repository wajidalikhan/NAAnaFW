#!/bin/sh
grep -E 'root' < work/files/TT.txt | cut -d ':' -f2 | awk -v prefix="root://cms-xrd-global.cern.ch/" '{print prefix $0}' >  work/files/TT_1.txt
mv work/files/TT_1.txt work/files/TT.txt
