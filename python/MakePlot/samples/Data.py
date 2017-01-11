######################################
#
# Deborah Pinna, August 2015
#
######################################


from utils import *

#samples.append("SingleMuon_Run2016G-PromptReco-v1")
#samples.append("SingleMuon_Run2016F-PromptReco-v1")
#samples.append("SingleMuon_Run2016E-PromptReco-v2")
#samples.append("SingleMuon_Run2016D-PromptReco-v2")
#samples.append("SingleMuon_Run2016C-PromptReco-v2")
#samples.append("SingleMuon_Run2016B-PromptReco-v2")


SingleMu_Run2016G = sample()
SingleMu_Run2016G.files = outlist (d,"SingleMuon_Run2016G-PromptReco-v1")
SingleMu_Run2016G.skimEff = 1.0
SingleMu_Run2016G.sigma = -1
SingleMu_Run2016G.jpref = jetLabel 
SingleMu_Run2016G.jp = jetLabel
SingleMu_Run2016G.color = ROOT.kBlack
SingleMu_Run2016G.style = 1
SingleMu_Run2016G.fill = 1001
SingleMu_Run2016G.leglabel = "SingleMuon_Run2016G-PromptReco-v1"
SingleMu_Run2016G.label = "SingleMuon_Run2016G-PromptReco-v1"

SingleMu_Run2016F = sample()
SingleMu_Run2016F.files = outlist (d,"SingleMuon_Run2016F-PromptReco-v1")
SingleMu_Run2016F.skimEff = 1.0
SingleMu_Run2016F.sigma = -1
SingleMu_Run2016F.jpref = jetLabel 
SingleMu_Run2016F.jp = jetLabel
SingleMu_Run2016F.color = ROOT.kBlack
SingleMu_Run2016F.style = 1
SingleMu_Run2016F.fill = 1001
SingleMu_Run2016F.leglabel = "SingleMuon_Run2016F-PromptReco-v1"
SingleMu_Run2016F.label = "SingleMuon_Run2016F-PromptReco-v1"

SingleMu_Run2016E = sample()
SingleMu_Run2016E.files = outlist (d,"SingleMuon_Run2016E-PromptReco-v2")
SingleMu_Run2016E.skimEff = 1.0
SingleMu_Run2016E.sigma = -1
SingleMu_Run2016E.jpref = jetLabel 
SingleMu_Run2016E.jp = jetLabel
SingleMu_Run2016E.color = ROOT.kBlack
SingleMu_Run2016E.style = 1
SingleMu_Run2016E.fill = 1001
SingleMu_Run2016E.leglabel = "SingleMuon_Run2016E-PromptReco-v2"
SingleMu_Run2016E.label = "SingleMuon_Run2016E-PromptReco-v2"

SingleMu_Run2016D = sample()
SingleMu_Run2016D.files = outlist (d,"SingleMuon_Run2016D-PromptReco-v2")
SingleMu_Run2016D.skimEff = 1.0
SingleMu_Run2016D.sigma = -1
SingleMu_Run2016D.jpref = jetLabel 
SingleMu_Run2016D.jp = jetLabel
SingleMu_Run2016D.color = ROOT.kBlack
SingleMu_Run2016D.style = 1
SingleMu_Run2016D.fill = 1001
SingleMu_Run2016D.leglabel = "SingleMuon_Run2016D-PromptReco-v2"
SingleMu_Run2016D.label = "SingleMuon_Run2016D-PromptReco-v2"

SingleMu_Run2016C = sample()
SingleMu_Run2016C.files = outlist (d,"SingleMuon_Run2016C-PromptReco-v2")
SingleMu_Run2016C.skimEff = 1.0
SingleMu_Run2016C.sigma = -1
SingleMu_Run2016C.jpref = jetLabel 
SingleMu_Run2016C.jp = jetLabel
SingleMu_Run2016C.color = ROOT.kBlack
SingleMu_Run2016C.style = 1
SingleMu_Run2016C.fill = 1001
SingleMu_Run2016C.leglabel = "SingleMuon_Run2016C-PromptReco-v2"
SingleMu_Run2016C.label = "SingleMuon_Run2016C-PromptReco-v2"

SingleMu_Run2016B = sample()
SingleMu_Run2016B.files = outlist (d,"SingleMuon_Run2016B-PromptReco-v2")
SingleMu_Run2016B.skimEff = 1.0
SingleMu_Run2016B.sigma = -1
SingleMu_Run2016B.jpref = jetLabel 
SingleMu_Run2016B.jp = jetLabel
SingleMu_Run2016B.color = ROOT.kBlack
SingleMu_Run2016B.style = 1
SingleMu_Run2016B.fill = 1001
SingleMu_Run2016B.leglabel = "SingleMuon_Run2016B-PromptReco-v2"
SingleMu_Run2016B.label = "SingleMuon_Run2016B-PromptReco-v2"


Data = sample()
Data.color = ROOT.kBlack
Data.style = 1
Data.fill = 1001
Data.leglabel = "Data"
Data.label = "Data"
#Data.components = [    SingleMu_Run2016B,
#                       SingleMu_Run2016C,
#                       SingleMu_Run2016D,
#                       SingleMu_Run2016E,
#                       SingleMu_Run2016F,
#                       SingleMu_Run2016G]
Data.components = [    SingleMu_Run2016B,# 5.879141 
                       SingleMu_Run2016C,# 2.645968 
                       SingleMu_Run2016D# 4.353449 
#                       SingleMu_Run2016E,# 4.049732
#                       SingleMu_Run2016F,# 3.147823 
#                       SingleMu_Run2016G,# 7.115969 
]

#Data.components = [SingleMu_05Oct, SingleEl_05Oct, MET_05Oct, MET_Promptv4, SingleMu_Promptv4, SingleEl_Promptv4]
#Data.components = [SingleMu_05Oct, SingleEl_05Oct, MET_05Oct, MET_Promptv4]
#Data.components = [MET]



