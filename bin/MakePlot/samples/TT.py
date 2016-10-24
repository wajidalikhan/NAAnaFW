from utils import *

TT = sample()
TT.files = outlist (d,"TT")
TT.files = TT.files[:1]
TT.jpref = jetLabel 
TT.jp = jetLabel 
TT.skimEff = 1.
TT.sigma = 831.76
#TT.color = 9
#TT.color = ROOT.kMagenta + 2
#TT.color = ROOT.kMagenta - 9
TT.color = ROOT.kAzure - 4
TT.style = 1
TT.fill = 1001
TT.leglabel = "t#bar{t}"
TT.label = "TT"
