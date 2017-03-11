#python makePlots.py -l 1.0 -c muonic -r "../../bin/"
#python makePlots.py -c muonic -n 1 -r "../../bin/"
#python makePlots.py -l 1.0 -c muonic -r "../../bin/"
#python makePlots.py -c muon -n 1 -r "/afs/cern.ch/work/o/oiorio/public/xNAAnaFW2016"
python makePlots.py -c muon -n 1 -r "./"
python makePlots.py -c muon -l 35.82 -r "../../bin" --focusOn h_2j0t_jetpt40_leading
python makePlots.py -c muon -l 35.82 -r ../../bin/ --getSF combine --fileSF mlfitqcd.root --focusOn h_2j1t_mtwcut_etajprime
#Public folder: Wajid: /afs/cern.ch/work/w/wajid/public/ForSingleTop/res 
#Public folder: Orso: /afs/cern.ch/work/o/oiorio/public/xNAAnaFW2016/res

