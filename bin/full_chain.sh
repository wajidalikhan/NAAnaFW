#!/bin/sh
echo "Renaming the files ..."
python script_rename.py -s ~oiorio/public/xWajid/files/trees/nov09/ -d files/renamed

echo "Renaming the xrd ..."
python script_replacexrd.py -f ./files/renamed/ -o files/final -x xrootd.ba.infn.it -P ST,VV,VJ,TT,SingleMuon

echo "Submitting jobs to cluster ... Remember to merge the output with the script: merge_res.py -l ./res/ -P ST,TT,VJ,VV --rm True"
python new_singletop.py --t3batch -f files/final/ -P ST,TT,VJ,VV -S 10

#Running over the data
#python new_singletop.py --t3batch -f files/final/ -P SingleMuon -d DATA -S 10

