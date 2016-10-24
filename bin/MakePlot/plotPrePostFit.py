#!/bin/env python


########################################
###  Author: Annapaola de Cosa             
###  Code to produce pre e post fit plots  
###  from combine output (mlfit.root)
########################################

from plots.plotUtils import *



bkgs_fh = ["TT", "QCD", "SingleTop", "DY","WJets", "ZToNuNu", "otherBkg"]
bkgs_sl = ["TT", "SingleTop", "DY","WJets", "ZToNuNu", "otherBkg"]
bkgs_fhtt = ["TT", "SingleTop", "DY","WJets", "otherBkg"]
bkgs_sltt = ["TT", "SingleTop", "DY","WJets", "otherBkg"]
bkgs_slz = ["TT", "SingleTop", "DY","WJets", "ZToNuNu", "otherBkg"]

regions = {
    "fh":("metFinal","fullhadronic"),
    "fh_tt":("metFinal_SR_1lep","fullhadronic_CR_TT"),
    "fh_VJets":("metFinal_CR5","fullhadronic_CR_VJets"),
    "fh_WJets":("metFinal_CR6nw","fullhadronic_CR_WJets"),
    "sl":("metFinal", "semileptonic"),
    "sl_TT":("metFinal_2lep","semileptonic_CR_TT"),
    "sl_WJets":("metFinal_met_0btag","semileptonic_CR_WJets")
    }
    
#plot("prefit", regions["fh"], bkgs_fh)
plot("fit_b", regions["fh"], bkgs_fh)

#plot("prefit", regions["fh_tt"], bkgs_fhtt)
plot("fit_b", regions["fh_tt"], bkgs_fhtt)

#plot("prefit", regions["fh_VJets"], bkgs_fh)
plot("fit_b", regions["fh_VJets"], bkgs_fh)

#plot("prefit", regions["fh_WJets"], bkgs_fh)
plot("fit_b", regions["fh_WJets"], bkgs_fh)

#plot("prefit", regions["fh_ZJets"], bkgs_fhtt)
#plot("fit_b", regions["fh_ZJets"], bkgs_fhtt)

#plot("prefit", regions["sl"], bkgs_sl)
plot("fit_b", regions["sl"], bkgs_sl)

#plot("prefit", regions["sl_TT"], bkgs_sltt)
plot("fit_b", regions["sl_TT"], bkgs_sltt)

#plot("prefit", regions["sl_WJets"], bkgs_slz)
plot("fit_b", regions["sl_WJets"], bkgs_slz)




    
    




