from utils import *

TTWJetsToLNu = sample()
TTWJetsToLNu.files = outlist (d,"TTWJetsToLNu")
TTWJetsToLNu.jpref = jetLabel 
TTWJetsToLNu.jp = jetLabel 
TTWJetsToLNu.skimEff = 1.
TTWJetsToLNu.sigma = 0.2043
TTWJetsToLNu.color = ROOT.kYellow -9
TTWJetsToLNu.style = 1
TTWJetsToLNu.fill = 1001
TTWJetsToLNu.leglabel = "TTWJetsToLNu"
TTWJetsToLNu.label = "TTWJetsToLNu"

TTWJetsToQQ = sample()
TTWJetsToQQ.files = outlist (d,"TTWJetsToQQ")
TTWJetsToQQ.jpref = jetLabel 
TTWJetsToQQ.jp = jetLabel 
TTWJetsToQQ.skimEff = 1.
TTWJetsToQQ.sigma = 0.4062
TTWJetsToQQ.color = ROOT.kYellow -9
TTWJetsToQQ.style = 1
TTWJetsToQQ.fill = 1001
TTWJetsToQQ.leglabel = "TTWJetsToQQ"
TTWJetsToQQ.label = "TTWJetsToQQ"


TTZToLLNuNu = sample()
TTZToLLNuNu.files = outlist (d,"TTZToLLNuNu")
TTZToLLNuNu.jpref = jetLabel 
TTZToLLNuNu.jp = jetLabel 
TTZToLLNuNu.skimEff = 1.
TTZToLLNuNu.sigma = 0.2529
TTZToLLNuNu.color = ROOT.kYellow -9
TTZToLLNuNu.style = 1
TTZToLLNuNu.fill = 1001
TTZToLLNuNu.leglabel = "TTZToLLNuNu"
TTZToLLNuNu.label = "TTZToLLNuNu"

TTZToQQ = sample()
TTZToQQ.files = outlist (d,"TTZToQQ")
TTZToQQ.jpref = jetLabel 
TTZToQQ.jp = jetLabel 
TTZToQQ.skimEff = 1.
TTZToQQ.sigma = 0.5297
TTZToQQ.color = ROOT.kYellow -9
TTZToQQ.style = 1
TTZToQQ.fill = 1001
TTZToQQ.leglabel = "TTZToQQ"
TTZToQQ.label = "TTZToQQ"

otherBkg = sample()
#otherBkg.color = ROOT.kMagenta + 2
otherBkg.color = ROOT.kBlue + 1
#otherBkg.color = 6
otherBkg.style = 1
otherBkg.fill = 1001
otherBkg.leglabel = "ttV"#splitline{other}{Bkg}"
otherBkg.label = "otherBkg"
otherBkg.components = [ TTWJetsToLNu,TTWJetsToQQ,  TTZToLLNuNu, TTZToQQ]
