######################################
#
# Deborah Pinna, August 2015
#
######################################


from utils import *

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



