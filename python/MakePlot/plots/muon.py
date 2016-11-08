met_range_had = (200,480, True)
met_rebin_had = 2 
from common import met_range

settings = {
    #'cutFlow'                  : (None                           , 20, None, (-0.5,7.5) ),

    #    'h_nPV' : ('number of primary vertices' , )
    #    'h_nMu': ('Number of tight Muons', 13, 1, (-0.5,12.5) ),
    'h_nJets': ('Number of tight jets', 13, 1, (-0.5,12.5) ),
    'h_nbJets': ('Number of tight b-jets', 11, 1, (-0.5,10.5) ),
    'h_2j1t_topMass': ('Top mass in the 2-jets 1-tag', 200, 5, (100,500) ),


    "h_2j0t_jetpt40_leading": ("2j0t leading jet Pt ",250, 1, (0,500)),
    "h_2j0t_jetpt40_subleading": ("2j0t Sub. leading jet Pt ",250, 1, (0,500)),
    "h_2j0t_mtw": (     "2j0t mtw ",250,1,(0,500)),
  
  
    "h_2j1t_bjetpt": (  "2j1t b jet pt",250,1,(0,500)),
    "h_2j1t_jetpt": (   "2j1t light jet pt" ,250,1,(0,500)),
    "h_2j1t_ljeteta": (  "2j1t light jet eta",200,5,(-4.7,4.7)),
    "h_2j1t_nMu": (     "2j1t Number of tight Muons",13,1,(-0.5,12.5)),
    "h_2j1t_MuPt": (    "2j1t muon pt ",100,1,(250,500)),
    "h_2j1t_MuEta": (   "2j1t muon eta ",100,1,(-2.1,2.1)),
    "h_2j1t_MuPhi": (   "2j1t muon phi ",100,1,( -3.2, 3.2)),
    "h_2j1t_MuE": (     "2j1t muon e ",250,1,(0,500)),
    "h_2j1t_MuCharge": ("2j1t muon charge ",2,1,(-1,1)),
    "h_2j1t_mtwcut_mtw": (     "2j1t mtw ",250,5,(0,500)),
    "h_2j1t_mtw": (     "2j1t mtw ",250,1,(0,500)),
    "h_2j1t_mtwcut_topMass": ( "2j1t top mass",200,5,(100,500)),
    "h_2j1t_topPt": (   "2j1t top pt",250,1,(0,500)),
    
    
    "h_3j1t_bjetpt": (  "3j1t b jet pt ",250,1,(0,500)),
    "h_3j1t_MuPt": (    "muon pt ",250,1,(0,500)),
    "h_3j1t_MuEta": (   "muon eta ",100,1,(-2.1,2.1)),
    "h_3j1t_MuPhi": (   "muon phi ",100, 1,(-3.2, 3.2)),
    "h_3j1t_MuE": (     "muon e ",100,1,(0,500)),
    "h_3j1t_mtw": (     "3j1t mtw ",100,1,(0,500)),
    
    
    "h_3j2t_bjetpt": (     "3j2t b jet pt ",250,1,(0,500)),
    "h_3j2t_2ndbjetpt": (  "3j2t sub lead. b jet pt",250,1,(0,500)),
    "h_3j2t_MuPt": (       "3j2t muon pt ",250,1,(0,500)),
    "h_3j2t_MuEta": (      "3j2t muon eta ",100,1,(-2.1,2.1)),
    "h_3j2t_MuPhi": (      "3j2t muon phi ",100, 1,(-3.2, 3.2)),
    "h_3j2t_MuE": (        "3j2t muon e ",250,1,(0,500)),
    "h_3j2t_MuCharge": (   "3j2t muon charge ",2,1,(-1,1)),
    "h_3j2t_mtw": (        "3j2t mtw ",250,1,(0,500)),
    

#    'dphi_CR5'     : ('#Delta #phi (j_{1,2},E^{miss}_{T})', 20, 2, None ),
#    'dphi_6j_CR5'  : ('#Delta #phi (j_{1,...,6},E^{miss}_{T})', 20, 2, None ),
#
#    'h_nJets_CR5'  : ('Number of tight jets'         , 20, 1, (4,12)),
#    'h_nbJets_CR5' : ('Number of tight b-jets'       , 20, 1, (2,6) ),
#    'h_jet1Pt_CR5' : ('Leading jet Pt distribution' , 20, 5, (30,490) ),
#    'h_jet2Pt_CR5' : ('Sub-leading jet Pt distribution'  , 20, 5, (30,490)),
#    'h_jet3Pt_CR5' : ('Third Jet Pt distribution'    , 20, 5, (30,230) ),
#
#    #CR6
#    #'metFinal_CR6'             : ('E^{miss}_{T} [GeV]'                          , 20, met_rebin_had, met_range_had ),
#    'muonmetFinal_CR6'             : ('E^{miss}_{T} [GeV]'                          , 20, met_rebin_had, met_range_had ),
#    'electronmetFinal_CR6'             : ('E^{miss}_{T} [GeV]'                          , 20, met_rebin_had, met_range_had ),
#
#    'dphi_CR6'         : ('#Delta #phi (j_{1,2},E^{miss}_{T})', 20, 2, None ),
#    'dphi_6j_CR6'         : ('#Delta #phi (j_{1,...,6},E^{miss}_{T})', 20, 2, None ),
#    'muondphi_CR6'         : ('#Delta #phi (j_{1,2},E^{miss}_{T})', 20, 2, None ),
#    'muondphi_6j_CR6'      : ('#Delta #phi (j_{1,...,6},E^{miss}_{T})', 20, 2, None ),
#    'lectrondphi_CR6'      : ('#Delta #phi (j_{1,2},E^{miss}_{T})', 20, 2, None ),
#    'lectrondphi_6j_CR6'   : ('#Delta #phi (j_{1,...,6},E^{miss}_{T})', 20, 2, None ),
#
#    'h_nJets_CR6'        : ('Number of tight jets'         , 20, 1, (4,12)),
#    'h_nbJets_CR6'       : ('Number of tight b-jets'       , 20, 1, (2,6) ),
#    'h_muonnJets_CR6'      : ('Number of tight jets'         , 20, 1, (4,12)),
#    'h_muonnbJets_CR6'     : ('Number of tight b-jets'       , 20, 1, (2,6) ),
#    'h_electronnJets_CR6'  : ('Number of tight jets'         , 20, 1, (4,12)),
#    'h_electronnbJets_CR6' : ('Number of tight b-jets'       , 20, 1, (2,6) ),
#
#    'h_jet1Pt_CR6'       : ('Leading jet Pt distribution' , 20, 20, (30,490) ),
#    'h_jet2Pt_CR6'       : ('Sub-leading jet Pt distribution'  , 20, 20, (30,490)),
#    'h_jet3Pt_CR6'       : ('Third Jet Pt distribution'    , 20, 20, (30,230) ),
#    'h_muonjet1Pt_CR6'     : ('Leading jet Pt distribution' , 20, 5, (30,490) ),
#    'h_muonjet2Pt_CR6'     : ('Sub-leading jet Pt distribution'  , 20, 5, (30,490)),
#    'h_muonjet3Pt_CR6'     : ('Third Jet Pt distribution'    , 20, 5, (30,230) ),
#    'h_electronjet1Pt_CR6' : ('Leading jet Pt distribution' , 20, 5, (30,490) ),
#    'h_electronjet2Pt_CR6' : ('Sub-leading jet Pt distribution'  , 20, 5, (30,490)),
#    'h_electronjet3Pt_CR6' : ('Third Jet Pt distribution'    , 20, 5, (30,230) ),
#
#    #CR6nw
#    'metFinal_CR6nw'             : ('E^{miss}_{T} [GeV]'                          , 1000, met_rebin_had, met_range_had ),
#    'muonmetFinal_CR6nw'             : ('E^{miss}_{T} [GeV]'                          , 20, met_rebin_had, met_range_had ),
#    'electronmetFinal_CR6nw'             : ('E^{miss}_{T} [GeV]'                          , 20, met_rebin_had, met_range_had ),
#
#    'dphi_CR6nw'         : ('#Delta #phi (j_{1,2},E^{miss}_{T})', 20, 2, None ),
#    'dphi_6j_CR6nw'         : ('#Delta #phi (j_{1,...,6},E^{miss}_{T})', 20, 2, None ),
#    'muondphi_CR6nw'         : ('#Delta #phi (j_{1,2},E^{miss}_{T})', 20, 2, None ),
#    'muondphi_6j_CR6nw'      : ('#Delta #phi (j_{1,...,6},E^{miss}_{T})', 20, 2, None ),
#    'lectrondphi_CR6nw'      : ('#Delta #phi (j_{1,2},E^{miss}_{T})', 20, 2, None ),
#    'lectrondphi_6j_CR6nw'   : ('#Delta #phi (j_{1,...,6},E^{miss}_{T})', 20, 2, None ),
#
#    'h_nJets_CR6nw'        : ('Number of tight jets'         , 20, 1, (4,12)),
#    'h_nbJets_CR6nw'       : ('Number of tight b-jets'       , 20, 1, (2,6) ),
#    'h_muonnJets_CR6nw'      : ('Number of tight jets'         , 20, 1, (4,12)),
#    'h_muonnbJets_CR6nw'     : ('Number of tight b-jets'       , 20, 1, (2,6) ),
#    'h_electronnJets_CR6nw'  : ('Number of tight jets'         , 20, 1, (4,12)),
#    'h_electronnbJets_CR6nw' : ('Number of tight b-jets'       , 20, 1, (2,6) ),
#
#    'h_jet1Pt_CR6nw'       : ('Leading jet Pt distribution' , 20, 20, (30,490) ),
#    'h_jet2Pt_CR6nw'       : ('Sub-leading jet Pt distribution'  , 20, 20, (30,490)),
#    'h_jet3Pt_CR6nw'       : ('Third Jet Pt distribution'    , 20, 20, (30,230) ),
#    'h_muonjet1Pt_CR6nw'     : ('Leading jet Pt distribution' , 20, 5, (30,490) ),
#    'h_muonjet2Pt_CR6nw'     : ('Sub-leading jet Pt distribution'  , 20, 5, (30,490)),
#    'h_muonjet3Pt_CR6nw'     : ('Third Jet Pt distribution'    , 20, 5, (30,230) ),
#    'h_electronjet1Pt_CR6nw' : ('Leading jet Pt distribution' , 20, 5, (30,490) ),
#    'h_electronjet2Pt_CR6nw' : ('Sub-leading jet Pt distribution'  , 20, 5, (30,490)),
#    'h_electronjet3Pt_CR6nw' : ('Third Jet Pt distribution'    , 20, 5, (30,230) ),
#
#    #CR7
#    #'metFinal_CR7'             : ('E^{miss}_{T} [GeV]'                          , 20, met_rebin_had, met_range_had ),
#    'muonmetFinal_CR7'             : ('E^{miss}_{T} [GeV]'                          , 20, met_rebin_had, met_range_had ),
#    'electronmetFinal_CR7'             : ('E^{miss}_{T} [GeV]'                          , 20, met_rebin_had, met_range_had ),
#    #'mixmetFinal_CR7'             : ('E^{miss}_{T} [GeV]'                          , 20, met_rebin_had, met_range_had ),
#
#    'dphi_CR7'             : ('#Delta #phi (j_{1,2},E^{miss}_{T})', 20, 2, None ),
#    'dphi_6j_CR7'          : ('#Delta #phi (j_{1,...,6},E^{miss}_{T})', 20, 2, None ),
#    'muondphi_CR7'         : ('#Delta #phi (j_{1,2},E^{miss}_{T})', 20, 2, None ),
#    'muondphi_6j_CR7'      : ('#Delta #phi (j_{1,...,6},E^{miss}_{T})', 20, 2, None ),
#    'electrondphi_CR7'     : ('#Delta #phi (j_{1,2},E^{miss}_{T})', 20, 2, None ),
#    'electrondphi_6j_CR7'  : ('#Delta #phi (j_{1,...,6},E^{miss}_{T})', 20, 2, None ),
#    #'mixdphi_CR7'         : ('#Delta #phi (j_{1,2},E^{miss}_{T})', 20, 1, None ),
#    #'mixdphi_6j_CR7'         : ('#Delta #phi (j_{1,...,6},E^{miss}_{T})', 20, 1, None ),
#    
#    'h_mu1Pt_CR7'     : ('Leading muon Pt distribution' , 20, 40, (30,230) ),
#    'h_mu2Pt_CR7'     : ('Sub-leading muon Pt distribution' , 20, 40, (30,230) ),
#    'h_el1Pt_CR7'     : ('Leading electron Pt distribution' , 20, 40, (30,230) ),
#    'h_el2Pt_CR7'     : ('Sub-leading electron Pt distribution' , 20, 40, (30,230) ),
#    
#    'h_nJets_CR7'          : ('Number of tight jets'         , 20, 1, (4,12)),
#    'h_nbJets_CR7'         : ('Number of tight b-jets'       , 20, 1, (2,6) ),
#    'h_muonnJets_CR7'      : ('Number of tight jets'         , 20, 1, (4,12)),
#    'h_muonnbJets_CR7'     : ('Number of tight b-jets'       , 20, 1, (2,6) ),
#    'h_electronnJets_CR7'  : ('Number of tight jets'         , 20, 1, (4,12)),
#    'h_electronnbJets_CR7' : ('Number of tight b-jets'       , 20, 1, (2,6) ),
#    'h_mixnJets_CR7'        : ('Number of tight jets'         , 20, 1, (4,12)),
#    'h_mixnbJets_CR7'       : ('Number of tight b-jets'       , 20, 1, (2,6) ),
#    
#    'h_jet1Pt_CR7'         : ('Leading jet Pt distribution' , 20, 50, (30,490) ),
#    'h_jet2Pt_CR7'         : ('Sub-leading jet Pt distribution'  , 20, 50, (30,490)),
#    'h_jet3Pt_CR7'         : ('Third Jet Pt distribution'    , 20, 50, (30,230) ),
#    'h_muonjet1Pt_CR7'     : ('Leading jet Pt distribution' , 20, 50, (30,490) ),
#    'h_muonjet2Pt_CR7'     : ('Sub-leading jet Pt distribution'  , 20, 50, (30,490)),
#    'h_muonjet3Pt_CR7'     : ('Third Jet Pt distribution'    , 20, 50, (30,230) ),
#    'h_electronjet1Pt_CR7' : ('Leading jet Pt distribution' , 20, 50, (30,490) ),
#    'h_electronjet2Pt_CR7' : ('Sub-leading jet Pt distribution'  , 20, 50, (30,490)),
#    'h_electronjet3Pt_CR7' : ('Third Jet Pt distribution'    , 20, 50, (30,230) ),
#    'h_mixjet1Pt_CR7'       : ('Leading jet Pt distribution' , 20, 20, (30,490) ),
#    'h_mixjet2Pt_CR7'       : ('Sub-leading jet Pt distribution'  , 20, 20, (30,490)),
#    'h_mixjet3Pt_CR7'       : ('Third Jet Pt distribution'    , 20, 20, (30,230) ),
#
#    #CR7nw
#    'metFinal_CR7nw'             : ('E^{miss}_{T} [GeV]'                          , 1000, met_rebin_had, met_range_had ),
#    'muonmetFinal_CR7nw'             : ('E^{miss}_{T} [GeV]'                          , 20, met_rebin_had, met_range_had ),
#    'electronmetFinal_CR7nw'             : ('E^{miss}_{T} [GeV]'                          , 20, met_rebin_had, met_range_had ),
#    'mixmetFinal_CR7nw'             : ('E^{miss}_{T} [GeV]'                          , 20, met_rebin_had, met_range_had ),
#
#    'dphi_CR7nw'             : ('#Delta #phi (j_{1,2},E^{miss}_{T})', 20, 2, None ),
#    'dphi_6j_CR7nw'          : ('#Delta #phi (j_{1,...,6},E^{miss}_{T})', 20, 2, None ),
#    'muondphi_CR7nw'         : ('#Delta #phi (j_{1,2},E^{miss}_{T})', 20, 2, None ),
#    'muondphi_6j_CR7nw'      : ('#Delta #phi (j_{1,...,6},E^{miss}_{T})', 20, 2, None ),
#    'electrondphi_CR7nw'     : ('#Delta #phi (j_{1,2},E^{miss}_{T})', 20, 2, None ),
#    'electrondphi_6j_CR7nw'  : ('#Delta #phi (j_{1,...,6},E^{miss}_{T})', 20, 2, None ),
#    #'mixdphi_CR7nw'         : ('#Delta #phi (j_{1,2},E^{miss}_{T})', 20, 1, None ),
#    #'mixdphi_6j_CR7nw'         : ('#Delta #phi (j_{1,...,6},E^{miss}_{T})', 20, 1, None ),
#    
#    'h_mu1Pt_CR7nw'     : ('Leading muon Pt distribution' , 20, 40, (30,230) ),
#    'h_mu2Pt_CR7nw'     : ('Sub-leading muon Pt distribution' , 20, 40, (30,230) ),
#    'h_el1Pt_CR7nw'     : ('Leading electron Pt distribution' , 20, 40, (30,230) ),
#    'h_el2Pt_CR7nw'     : ('Sub-leading electron Pt distribution' , 20, 40, (30,230) ),
#    
#    'h_nJets_CR7nw'          : ('Number of tight jets'         , 20, 1, (4,12)),
#    'h_nbJets_CR7nw'         : ('Number of tight b-jets'       , 20, 1, (2,6) ),
#    'h_muonnJets_CR7nw'      : ('Number of tight jets'         , 20, 1, (4,12)),
#    'h_muonnbJets_CR7nw'     : ('Number of tight b-jets'       , 20, 1, (2,6) ),
#    'h_electronnJets_CR7nw'  : ('Number of tight jets'         , 20, 1, (4,12)),
#    'h_electronnbJets_CR7nw' : ('Number of tight b-jets'       , 20, 1, (2,6) ),
#    'h_mixnJets_CR7nw'        : ('Number of tight jets'         , 20, 1, (4,12)),
#    'h_mixnbJets_CR7nw'       : ('Number of tight b-jets'       , 20, 1, (2,6) ),
#    
#    'h_jet1Pt_CR7nw'         : ('Leading jet Pt distribution' , 20, 50, (30,490) ),
#    'h_jet2Pt_CR7nw'         : ('Sub-leading jet Pt distribution'  , 20, 50, (30,490)),
#    'h_jet3Pt_CR7nw'         : ('Third Jet Pt distribution'    , 20, 50, (30,230) ),
#    'h_muonjet1Pt_CR7nw'     : ('Leading jet Pt distribution' , 20, 50, (30,490) ),
#    'h_muonjet2Pt_CR7nw'     : ('Sub-leading jet Pt distribution'  , 20, 50, (30,490)),
#    'h_muonjet3Pt_CR7nw'     : ('Third Jet Pt distribution'    , 20, 50, (30,230) ),
#    'h_electronjet1Pt_CR7nw' : ('Leading jet Pt distribution' , 20, 50, (30,490) ),
#    'h_electronjet2Pt_CR7nw' : ('Sub-leading jet Pt distribution'  , 20, 50, (30,490)),
#    'h_electronjet3Pt_CR7nw' : ('Third Jet Pt distribution'    , 20, 50, (30,230) ),
#    'h_mixjet1Pt_CR7nw'       : ('Leading jet Pt distribution' , 20, 20, (30,490) ),
#    'h_mixjet2Pt_CR7nw'       : ('Sub-leading jet Pt distribution'  , 20, 20, (30,490)),
#    'h_mixjet3Pt_CR7nw'       : ('Third Jet Pt distribution'    , 20, 20, (30,230) ),
#
}


store = [
    "h_nJets", #SR shape
    "h_2j1t_topMass",
    "h_nbJets", #SR shape

    "h_2j0t_jetpt40_leading", 
    "h_2j0t_jetpt40_subleading", 
    "h_2j0t_mtw", 
  
  
    "h_2j1t_bjetpt", 
    "h_2j1t_jetpt", 
    "h_2j1t_jeteta", 
    "h_2j1t_nMu", 
    "h_2j1t_MuPt", 
    "h_2j1t_MuEta", 
    "h_2j1t_MuPhi", 
    "h_2j1t_MuE", 
    "h_2j1t_MuCharge", 
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
    "h_3j2t_MuCharge", 
    "h_3j2t_mtw", 
    ]
