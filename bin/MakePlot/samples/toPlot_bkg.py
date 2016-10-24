######################################
#
# Annapaola de Cosa, January 2015
#
######################################

import ROOT
import collections
import os, commands
from SingleTopWJets import *
from ttDM import *
from TT import *
from DYJetsToLL import *
from ZJetsToNuNu import *
from QCD import *
from Data import *
from otherBkg import *

samples = collections.OrderedDict()

#samples["MET_Prompt"] = MET_Prompt
#samples["MET_ReMiniAOD"] = MET_ReMiniAOD17Jul
#samples["SingleMu_MiniAOD17Jul"] = SingleMu_MiniAOD17Jul
#samples["SingleMu_Prompt"] = SingleMu_Prompt
#samples["SingleEl_MiniAOD17Jul"] = SingleEl_MiniAOD17Jul
#samples["SingleEl_Prompt"] = SingleEl_Prompt
#samples["Data"] = Data
#samples["QCD"] = QCD
#samples["DY"] = DY
#samples["ZToNuNu"] = ZToNuNu
#samples["SingleTop"] = SingleTop
#samples["WJets"] = WJets
#samples["TT"] = TT
#samples["otherBkg"] = otherBkg

samples["DMtt_sc_Mchi1Mphi10"] = DMtt_sc_Mchi1Mphi10
samples["DMtt_sc_Mchi1Mphi20"] = DMtt_sc_Mchi1Mphi20
