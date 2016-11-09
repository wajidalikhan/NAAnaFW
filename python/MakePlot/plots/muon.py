met_range_had = (200,480, True)
met_rebin_had = 2
common_rebin = 5 
from common import met_range

settings = {
    #'cutFlow'                  : (None                           , 20, None, (-0.5,7.5) ),
    #    'h_nPV' : ('number of primary vertices' , )
    #    'h_nMu': ('Number of tight Muons', 13, 1, (-0.5,12.5) ),
    #'h_nJets': ('Number of tight jets', 13, 1, (-0.5,12.5) ),
    #'h_nbJets': ('Number of tight b-jets', 11, 1, (-0.5,10.5) ),
    'h_2j1t_topMass': ('Top mass in the 2-jets 1-tag', 200, 5, (100,500) ),

    "h_2j0t_jetpt40_leading": ("2j0t leading jet Pt ",250, common_rebin, (0,500)),
    "h_2j0t_jetpt40_subleading": ("2j0t Sub. leading jet Pt ",250, common_rebin, (0,500)),
    "h_2j0t_mtw": (     "2j0t mtw ",250,common_rebin,(0,500)),
  
    "h_2j1t_bjetpt": (  "2j1t b jet pt",250,common_rebin,(0,500)),
    "h_2j1t_jetpt": (   "2j1t light jet pt" ,250,common_rebin,(0,500)),
    "h_2j1t_ljeteta": (  "2j1t light jet eta",200,5,(-4.7,4.7)),
    #"h_2j1t_nMu": (     "2j1t Number of tight Muons",13,1,(-0.5,12.5)),
    "h_2j1t_MuPt": (    "2j1t muon pt ",100,common_rebin,(0,500)),
    "h_2j1t_MuEta": (   "2j1t muon eta ",100,common_rebin,(-2.1,2.1)),
    "h_2j1t_MuPhi": (   "2j1t muon phi ",100,common_rebin,( -3.2, 3.2)),
    "h_2j1t_MuE": (     "2j1t muon e ",250,common_rebin,(0,500)),
    #"h_2j1t_MuCharge": ("2j1t muon charge ",2,1,(-1,1)),
    "h_2j1t_mtwcut_mtw": (     "2j1t mtw ",250,common_rebin,(0,500)),
    "h_2j1t_mtw": (     "2j1t mtw ",250,common_rebin,(0,500)),
    "h_2j1t_mtwcut_topMass": ( "2j1t top mass",200,5,(100,500)),
    "h_2j1t_topPt": (   "2j1t top pt",250,common_rebin,(0,500)),
    
    "h_3j1t_bjetpt": (  "3j1t b jet pt ",250,common_rebin,(0,500)),
    "h_3j1t_MuPt": (    "muon pt ",250,common_rebin,(0,500)),
    "h_3j1t_MuEta": (   "muon eta ",100,common_rebin,(-2.1,2.1)),
    "h_3j1t_MuPhi": (   "muon phi ",100, common_rebin,(-3.2, 3.2)),
    "h_3j1t_MuE": (     "muon e ",100,common_rebin,(0,500)),
    "h_3j1t_mtw": (     "3j1t mtw ",100,common_rebin,(0,500)),
    
    "h_3j2t_bjetpt": (     "3j2t b jet pt ",250,common_rebin,(0,500)),
    "h_3j2t_2ndbjetpt": (  "3j2t sub lead. b jet pt",250,common_rebin,(0,500)),
    "h_3j2t_MuPt": (       "3j2t muon pt ",250,common_rebin,(0,500)),
    "h_3j2t_MuEta": (      "3j2t muon eta ",100,common_rebin,(-2.1,2.1)),
    "h_3j2t_MuPhi": (      "3j2t muon phi ",100,common_rebin,(-3.2, 3.2)),
    "h_3j2t_MuE": (        "3j2t muon e ",250,common_rebin,(0,500)),
    #"h_3j2t_MuCharge": (   "3j2t muon charge ",2,1,(-1,1)),
    "h_3j2t_mtw": (        "3j2t mtw ",250,common_rebin,(0,500)),
}

store = [
    #"h_nJets", #SR shape
    "h_2j1t_topMass",
    #"h_nbJets", #SR shape
    "h_2j0t_jetpt40_leading", 
    "h_2j0t_jetpt40_subleading", 
    "h_2j0t_mtw", 
    "h_2j1t_bjetpt", 
    "h_2j1t_jetpt", 
    "h_2j1t_jeteta", 
    #"h_2j1t_nMu", 
    "h_2j1t_MuPt", 
    "h_2j1t_MuEta", 
    "h_2j1t_MuPhi", 
    "h_2j1t_MuE", 
    #"h_2j1t_MuCharge", 
    "h_2j1t_mtwcut_mtw", 
    "h_2j1t_mtw", 
    "h_2j1t_topMass", 
    "h_2j1t_mtwcut_topMass",
    "h_2j1t_topPt", 
    "h_3j1t_bjetpt", 
    "h_3j1t_MuPt", 
    "h_3j1t_MuEta", 
    "h_3j1t_MuPhi", 
    "h_3j1t_MuE",
    "h_3j1t_mtw", 
    "h_3j2t_bjetpt", 
    "h_3j2t_2ndbjetpt", 
    "h_3j2t_MuPt", 
    "h_3j2t_MuEta", 
    "h_3j2t_MuPhi", 
    "h_3j2t_MuE", 
    #"h_3j2t_MuCharge", 
    "h_3j2t_mtw", 
    ]
