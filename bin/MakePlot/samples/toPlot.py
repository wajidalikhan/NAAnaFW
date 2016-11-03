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
from VV import *
from Data import *

samples = collections.OrderedDict()

samples["Data"] = Data
samples["SingleTop_tchannel"] = ST_tch
samples["SingleTop_schannel"] = ST_sch
#samples["SingleTop_tW"] = ST_tW
samples["WJets"] = WJets
samples["DYJets"] = DYJets
samples["VV"] = VV
samples["TT"] = TT





