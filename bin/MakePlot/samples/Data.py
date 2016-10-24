######################################
#
# Deborah Pinna, August 2015
#
######################################


from utils import *

MET_05Oct = sample()
MET_05Oct.files = outlist (d,"MET_05Oct")
MET_05Oct.skimEff = 1
MET_05Oct.sigma = 1
MET_05Oct.color = ROOT.kBlack
MET_05Oct.jpref = jetLabel 
MET_05Oct.jp = jetLabel
MET_05Oct.style = 1
MET_05Oct.fill = 1001
MET_05Oct.leglabel = "MET_05Oct"
MET_05Oct.label = "MET_05Oct"

MET_Promptv4 = sample()
MET_Promptv4.files = outlist (d,"MET_Promptv4")
MET_Promptv4.skimEff = 1
MET_Promptv4.sigma = 1
MET_Promptv4.color = ROOT.kBlack
MET_Promptv4.jpref = jetLabel
MET_Promptv4.jp = jetLabel
MET_Promptv4.style = 1
MET_Promptv4.fill = 1001
MET_Promptv4.leglabel = "MET_Promptv4"
MET_Promptv4.label = "MET_Promptv4"


SingleEl_Promptv4 = sample()
SingleEl_Promptv4.files = outlist (d,"SingleEl_Promptv4")
SingleEl_Promptv4.skimEff = 1
SingleEl_Promptv4.sigma = 1
SingleEl_Promptv4.color = ROOT.kBlack
SingleEl_Promptv4.jpref = jetLabel 
SingleEl_Promptv4.jp = jetLabel
SingleEl_Promptv4.style = 1
SingleEl_Promptv4.fill = 1001
SingleEl_Promptv4.leglabel = "SingleEl_Promptv4"
SingleEl_Promptv4.label = "SingleEl_Promptv4"


SingleEl_05Oct = sample()
SingleEl_05Oct.files = outlist (d,"SingleEl_05Oct")
SingleEl_05Oct.skimEff = 1
SingleEl_05Oct.sigma = 1
SingleEl_05Oct.color = ROOT.kBlack
SingleEl_05Oct.jpref = jetLabel 
SingleEl_05Oct.jp = jetLabel
SingleEl_05Oct.style = 1
SingleEl_05Oct.fill = 1001
SingleEl_05Oct.leglabel = "SingleEl_05Oct"
SingleEl_05Oct.label = "SingleEl_05Oct"

SingleMu_05Oct = sample()
SingleMu_05Oct.files = outlist (d,"SingleMu_05Oct")
SingleMu_05Oct.skimEff = 1.0
SingleMu_05Oct.sigma = 1
SingleMu_05Oct.jpref = jetLabel 
SingleMu_05Oct.jp = jetLabel
SingleMu_05Oct.color = ROOT.kBlack
SingleMu_05Oct.style = 1
SingleMu_05Oct.fill = 1001
SingleMu_05Oct.leglabel = "SingleMu_05Oct"
SingleMu_05Oct.label = "SingleMu_05Oct"

SingleMu_Promptv4 = sample()
SingleMu_Promptv4.files = outlist (d,"SingleMu_Promptv4")
SingleMu_Promptv4.skimEff = 1.0
SingleMu_Promptv4.sigma = 1
SingleMu_Promptv4.jpref = jetLabel 
SingleMu_Promptv4.jp = jetLabel
SingleMu_Promptv4.color = ROOT.kBlack
SingleMu_Promptv4.style = 1
SingleMu_Promptv4.fill = 1001
SingleMu_Promptv4.leglabel = "SingleMu_Promptv4"
SingleMu_Promptv4.label = "SingleMu_Promptv4"

SingleMu_Run2016B = sample()
SingleMu_Run2016B.files = outlist (d,"SingleMu_Run2016B")
SingleMu_Run2016B.skimEff = 1.0
SingleMu_Run2016B.sigma = 1
SingleMu_Run2016B.jpref = jetLabel 
SingleMu_Run2016B.jp = jetLabel
SingleMu_Run2016B.color = ROOT.kBlack
SingleMu_Run2016B.style = 1
SingleMu_Run2016B.fill = 1001
SingleMu_Run2016B.leglabel = "SingleMu_Run2016B"
SingleMu_Run2016B.label = "SingleMu_Run2016B"


Data = sample()
Data.color = ROOT.kBlack
Data.style = 1
Data.fill = 1001
Data.leglabel = "Data"
Data.label = "Data"
Data.components = [SingleMu_Run2016B]
#Data.components = [SingleMu_05Oct, SingleEl_05Oct, MET_05Oct, MET_Promptv4, SingleMu_Promptv4, SingleEl_Promptv4]
#Data.components = [SingleMu_05Oct, SingleEl_05Oct, MET_05Oct, MET_Promptv4]
#Data.components = [MET]



