#python makePlots.py -l 1.0 -c muonic -r "../../bin/"
#python makePlots.py -c muonic -n 1 -r "../../bin/"
#python makePlots.py -l 1.0 -c muonic -r "../../bin/"
#python makePlots.py -c muon -n 1 -r "/afs/cern.ch/work/o/oiorio/public/xNAAnaFW2016"
python makePlots.py -c muon -n 1 -r "./"
python makePlots.py -c muon -l 27.1 -r "../../bin" --focusOn h_2j0t_jetpt40_leading
python makePlots.py -c muon -l 12.89 -r ../../bin/ --getSF combine --fileSF mlfitqcd.root --focusOn h_2j1t_mtwcut_etajprime

python makePlots.py -c muon -l 35.89 -r ./prob --getSF combine --fileSF mlfiteta.root --focusOn h_2j1t_mtwcut_le_sr_etajprime
python makePlots.py -c muon -l 35.89 -r muonch/ --getSF combine --fileSF mlfitqcd.root --focusOn h_2j1t_mtw_le50
