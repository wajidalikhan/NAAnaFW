python new_singletop.py -c muon -s noSys -C CSVT -m local -l /afs/cern.ch/work/o/oiorio/public/xWajid/synch/trees/mc --sync sync -P _synchMC
python new_singletop.py -c muon -s noSys -C CSVT -m local -l /afs/cern.ch/work/o/oiorio/public/xWajid/synch/trees/data --sync sync -P _synchDATA -d DATA

python new_singletop.py -c muonantiiso -m local -l ../test/mcqcd/ --sync sync -P _synchQCD
python new_singletop.py -c electron -m local -l ../test/signal4/ --sync sync -P _synchEl
python new_singletop.py -c muon -m local -l ../test/signal4/ --sync sync -P _synchMu
