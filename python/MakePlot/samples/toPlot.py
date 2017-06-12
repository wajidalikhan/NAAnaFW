######################################
#
# Annapaola de Cosa, January 2015
#
######################################

import ROOT
import collections
import os, commands
from SingleTop import *
from VJets import *
from TT import *
#from QCD import *
from DDQCD import *
from VV import *
from WZVV import *
from TT_tW import *
from Data import *

samples = collections.OrderedDict()
samples["Data"] = Data

#----------- Samples Added together 
samples["DDQCD"] = DDQCD
samples["WZJetsVV"] = WZJetsVV
samples["TT_tWST"] = TT_tWST
samples["SingleTop_tchannel"] = ST_tch
#-----------

#-----------Samples in original form
#samples["DDQCD"] = DDQCD
#samples["VV"] = VV
#samples["DYJets"] = DYJets
#samples["WJets"] = WJets
#samples["TT"] = TT
#samples["SingleTop_schannel"] = ST_sch
#samples["SingleTop_tW"] = ST_tW
#samples["SingleTop_tchannel"] = ST_tch
#samples["SingleTop_tchannel_sd"] = ST_tch_sd





