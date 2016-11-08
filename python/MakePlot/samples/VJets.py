
######################################
#
# Annapaola de Cosa, January 2015
#
######################################


from utils import *


WJetsToLNu = sample()
WJetsToLNu.files = outlist (d,"WJetsToLNu")
WJetsToLNu.skimEff = 1.
#WJetsToLNu.sigma = 23.1363
WJetsToLNu.sigma =  (11811.4+ 8677.3)*3.
WJetsToLNu.jpref = jetLabel 
WJetsToLNu.jp = jetLabel
WJetsToLNu.color = ROOT.kGreen+2
WJetsToLNu.style = 1
WJetsToLNu.fill = 1001
WJetsToLNu.leglabel = "WJets"
WJetsToLNu.label = "WJetsToLNu"

DYJetsToLL = sample()
DYJetsToLL.files = outlist (d,"DYJetsToLL")
DYJetsToLL.skimEff = 1.
#DYJetsToLL.sigma = 23.1363
DYJetsToLL.sigma =  1921.8*3.
DYJetsToLL.jpref = jetLabel 
DYJetsToLL.jp = jetLabel
DYJetsToLL.color = ROOT.kGreen+2
DYJetsToLL.style = 1
DYJetsToLL.fill = 1001
DYJetsToLL.leglabel = "DYJets"
DYJetsToLL.label = "DYJetsToLL"

WJets = sample()
#WJets.color = 6
#WJets.color = ROOT.kTeal -6
WJets.color = ROOT.kGreen - 2
WJets.style = 1
WJets.fill = 1001
WJets.leglabel = "W + Jets"
WJets.label = "WJets"
WJets.components = [WJetsToLNu]
#WJets.components = [WJets_HT100_200, WJets_HT200_400, WJets_HT400_600, WJets_HT600_800, WJets_HT800_1200, WJets_HT1200_2500, WJets_HT2500_Inf]
# WJets.components = [WJets_HT100_200, WJets_HT200_400,  WJets_HT800_1200, WJets_HT1200_2500 ]

DYJets = sample()
#DYJets.color = 6
DYJets.color = ROOT.kTeal -6
#DYJets.color = ROOT.kGreen - 2
DYJets.style = 1
DYJets.fill = 1001
DYJets.leglabel = "Z + Jets"
DYJets.label = "DYJets"
DYJets.components = [DYJetsToLL]




