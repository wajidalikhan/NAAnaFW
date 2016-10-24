from utils import *

DYToLL = sample()
DYToLL.files = outlist (d,"DY")
DYToLL.jpref = jetLabel 
DYToLL.jp = jetLabel 
DYToLL.skimEff = 1.
DYToLL.sigma = 213.5357
DYToLL.color = ROOT.kYellow -9
DYToLL.style = 1
DYToLL.fill = 1001
DYToLL.leglabel = "DY"
DYToLL.label = "DY"


DYToLL_HT100to200 = sample()
DYToLL_HT100to200.files = outlist (d,"DYJets_HT100to200")
DYToLL_HT100to200.jpref = jetLabel 
DYToLL_HT100to200.jp = jetLabel 
DYToLL_HT100to200.skimEff = 1.
DYToLL_HT100to200.sigma = 147.40
DYToLL_HT100to200.color = ROOT.kYellow -9
DYToLL_HT100to200.style = 1
DYToLL_HT100to200.fill = 1001
DYToLL_HT100to200.leglabel = "DY"
DYToLL_HT100to200.label = "DYJets_HT100to200"


DYToLL_HT200to400 = sample()
DYToLL_HT200to400.files = outlist (d,"DYJets_HT200to400")
DYToLL_HT200to400.jpref = jetLabel 
DYToLL_HT200to400.jp = jetLabel 
DYToLL_HT200to400.skimEff = 1.
DYToLL_HT200to400.sigma = 40.99
DYToLL_HT200to400.color = ROOT.kYellow -9
DYToLL_HT200to400.style = 1
DYToLL_HT200to400.fill = 1001
DYToLL_HT200to400.leglabel = "DY"
DYToLL_HT200to400.label = "DYJets_HT200to400"

DYToLL_HT400to600 = sample()
DYToLL_HT400to600.files = outlist (d,"DYJets_HT400to600")
DYToLL_HT400to600.jpref = jetLabel 
DYToLL_HT400to600.jp = jetLabel 
DYToLL_HT400to600.skimEff = 1.
DYToLL_HT400to600.sigma = 5.678
DYToLL_HT400to600.color = ROOT.kYellow -9
DYToLL_HT400to600.style = 1
DYToLL_HT400to600.fill = 1001
DYToLL_HT400to600.leglabel = "DY"
DYToLL_HT400to600.label = "DYJets_HT400to600"

DYToLL_HT600toInf = sample()
DYToLL_HT600toInf.files = outlist (d,"DYJets_HT600toInf")
DYToLL_HT600toInf.jpref = jetLabel 
DYToLL_HT600toInf.jp = jetLabel 
DYToLL_HT600toInf.skimEff = 1.
DYToLL_HT600toInf.sigma = 2.198
DYToLL_HT600toInf.color = ROOT.kYellow -9
DYToLL_HT600toInf.style = 1
DYToLL_HT600toInf.fill = 1001
DYToLL_HT600toInf.leglabel = "DY"
DYToLL_HT600toInf.label = "DYJets_HT600toInf"

DY = sample()
DY.color = 8
DY.style = 1
DY.fill = 1001
DY.leglabel = "Z#rightarrowll"
DY.label = "DY"
#DY.components = [ DYToLL]
DY.components = [ DYToLL_HT100to200,  DYToLL_HT200to400,  DYToLL_HT400to600, DYToLL_HT600toInf]
