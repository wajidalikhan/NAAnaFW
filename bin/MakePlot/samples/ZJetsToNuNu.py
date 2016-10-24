from utils import *


ZToNuNu = sample()
ZToNuNu.files = outlist (d,"ZJets")
ZToNuNu.jpref = jetLabel 
ZToNuNu.jp = jetLabel 
ZToNuNu.sf = 0.0747229
ZToNuNu.skimEff = 1.
ZToNuNu.sigma = 409.4874 
ZToNuNu.color = ROOT.kYellow -9
ZToNuNu.style = 1
ZToNuNu.fill = 1001
ZToNuNu.leglabel = "ZToNuNu"
ZToNuNu.label = "ZNuNu"


ZToNuNu_HT100to200 = sample()
ZToNuNu_HT100to200.files = outlist (d,"ZJets_HT100to200")
ZToNuNu_HT100to200.jpref = jetLabel 
ZToNuNu_HT100to200.jp = jetLabel 
ZToNuNu_HT100to200.sf = 0.0747229
ZToNuNu_HT100to200.skimEff = 1.
ZToNuNu_HT100to200.sigma = 280.35
#ZToNuNu_HT100to200.sigma = 409.4874
ZToNuNu_HT100to200.color = ROOT.kYellow -9
ZToNuNu_HT100to200.style = 1
ZToNuNu_HT100to200.fill = 1001
ZToNuNu_HT100to200.leglabel = "ZToNuNu"
ZToNuNu_HT100to200.label = "ZJets_HT100to200"

ZToNuNu_HT200to400 = sample()
ZToNuNu_HT200to400.files = outlist (d,"ZJets_HT200to400")
ZToNuNu_HT200to400.jpref = jetLabel 
ZToNuNu_HT200to400.jp = jetLabel 
ZToNuNu_HT200to400.sf = 0.0221710
ZToNuNu_HT200to400.skimEff = 1.
#ZToNuNu_HT200to400.sigma = 110.7792 
ZToNuNu_HT200to400.sigma = 77.67  
ZToNuNu_HT200to400.color = ROOT.kYellow -9
ZToNuNu_HT200to400.style = 1
ZToNuNu_HT200to400.fill = 1001
ZToNuNu_HT200to400.leglabel = "ZToNuNu"
ZToNuNu_HT200to400.label = "ZJets_HT200to400"

ZToNuNu_HT400to600 = sample()
ZToNuNu_HT400to600.files = outlist (d,"ZJets_HT400to600")
ZToNuNu_HT400to600.jpref = jetLabel 
ZToNuNu_HT400to600.jp = jetLabel 
ZToNuNu_HT400to600.sf = 0.0027042
ZToNuNu_HT400to600.skimEff = 1.
#ZToNuNu_HT400to600.sigma = 13.17701
ZToNuNu_HT400to600.sigma = 10.73
ZToNuNu_HT400to600.color = ROOT.kYellow -9
ZToNuNu_HT400to600.style = 1
ZToNuNu_HT400to600.fill = 1001
ZToNuNu_HT400to600.leglabel = "ZToNuNu"
ZToNuNu_HT400to600.label = "ZJets_HT400to600"

ZToNuNu_HT600toInf = sample()
ZToNuNu_HT600toInf.files = outlist (d,"ZJets_HT600toInf")
ZToNuNu_HT600toInf.jpref = jetLabel 
ZToNuNu_HT600toInf.jp = jetLabel 
ZToNuNu_HT600toInf.sf = 0.0009820
ZToNuNu_HT600toInf.skimEff = 1.
#ZToNuNu_HT600toInf.sigma = 4.520187
ZToNuNu_HT600toInf.sigma = 4.116 
ZToNuNu_HT600toInf.color = ROOT.kYellow -9
ZToNuNu_HT600toInf.style = 1
ZToNuNu_HT600toInf.fill = 1001
ZToNuNu_HT600toInf.leglabel = "ZToNuNu"
ZToNuNu_HT600toInf.label = "ZJets_HT600toInf"


ZToNuNu = sample()
ZToNuNu.color = ROOT.kGreen + 2
ZToNuNu.style = 1
ZToNuNu.fill = 1001
ZToNuNu.leglabel = "Z#rightarrow#nu#nu"
ZToNuNu.label = "ZToNuNu"
#ZToNuNu.components = [ ZToNuNu]
ZToNuNu.components = [ ZToNuNu_HT100to200,  ZToNuNu_HT200to400,  ZToNuNu_HT400to600, ZToNuNu_HT600toInf]

