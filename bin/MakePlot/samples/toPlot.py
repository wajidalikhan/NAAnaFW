######################################
#
# Annapaola de Cosa, January 2015
#
######################################

import ROOT
import collections
import os, commands
from SingleTopWJets import *
from TT import *
from Data import *

samples = collections.OrderedDict()

samples["Data"] = Data
samples["SingleTop_TChannel"] = SingleTop
#samples["SingleTop_tWChannel"] = SingleTop
samples["TT"] = TT





