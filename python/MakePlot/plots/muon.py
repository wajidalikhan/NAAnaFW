met_range_had = (200,480, True)
met_rebin_had = 2
common_rebin = 5 
rebin_sr = 10
from common import met_range

settings = {
#  TH1F *h_cutFlow = new TH1F("h_cutFlow","cutflow",10,-0.5,9.5);
    'h_cutFlow'    : ("cutFlow",10, 1, (-0.5,9.5) ),
    #(None                           , 20, None, (-0.5,7.5) ),
    #'h_nPV' : ('number of primary vertices' , )
    #    'h_nMu': ('Number of tight Muons', 13, 1, (-0.5,12.5) ),
    #'h_nJets': ('Number of tight jets', 13, 1, (-0.5,12.5) ),
    #'h_nbJets': ('Number of tight b-jets', 11, 1, (-0.5,10.5) ),
#    'h_2j1t_topMass': ('Top mass in the 2-jets 1-tag', 200, 5, (100,500) ),
    #2j0t
    "h_2j0t_jetpt40_leading": ("2j0t leading jet P_T ",125, common_rebin, (0,500)),
    "h_2j0t_jetpt40_subleading": ("2j0t Sub. leading jet Pt ",125, common_rebin, (0,250)),
    "h_2j0t_mtw": (     "2j0t mtw ",250,common_rebin*2,(0,250)),
    "h_2j0t_MuPt": (    "2j0t muon pt ",100,2,(0,500)),
    "h_2j0t_mtwcut_MuPt": (    "2j0t muon pt with mtwcut",100,2,(0,500)),
#    "h_2j0t_mtwcut_jeteta40_leading":( "2j0t leading jet eta", 200,common_rebin,(-4.7,4.7)),
    "h_2j0t_mtwcut_incljeteta":( "2j0t leading inclusive jet eta", 200,common_rebin,(-4.7,4.7)),


    "h_2j0t_mtwcut_jetabseta40_leading":( "2j0t leading jet eta", 100,common_rebin,(0.0,4.7)),
    "h_2j0t_mtwcut_incljetabseta":( "2j0t leading inclusive jet eta", 100,common_rebin,(0.0,4.7)),
    #"h_2j0t_muIso":("2j0t Muon Isolation",40,common_rebin,(0.0,1.0)),
    "h_2j0t_dR_lepjetpt40_1st":("2j0t #DeltaR(Lep,Jet) ",50,1,(0.0,5)),
    "h_2j0t_dPhi_lepjetpt40_1st":("2j0t #Delta#phi(Lep,Jet) ",50,1,(0.0,3.2)),
    "h_2j0t_dEta_lepjetpt40_1st":("2j0t #Delta#eta(Lep,Jet) ",50,2,(-3.2,3.2)),

    #2j1t  
    "h_2j1t_bjetpt": (  "2j1t b jet pt",250,common_rebin*2,(0,500)),
    "h_2j1t_ljetpt": (   "2j1t light jet pt" ,250,common_rebin*2,(0,500)),
    "h_2j1t_ljeteta": (  "2j1t light jet eta",100,5,(-4.7,4.7)),
    #"h_2j1t_nMu": (     "2j1t Number of tight Muons",13,1,(-0.5,12.5)),
    "h_2j1t_MuPt": (    "2j1t muon pt ",100,2,(26,100)),
    "h_2j1t_mtwcut_MuPt": (    "2j1t muon pt ",250,2,(0,500)),
    "h_2j1t_MuEta": (   "2j1t muon eta ",100,common_rebin,(-2.1,2.1)),
    "h_2j1t_MuPhi": (   "2j1t muon phi ",100,common_rebin,( -3.2, 3.2)),
    "h_2j1t_MuE": (     "2j1t muon e ",250,common_rebin,(0,500)),
    #"h_2j1t_MuCharge": ("2j1t muon charge ",2,1,(-1,1)),
    "h_2j1t_mtw": (     "2j1t mtw ",250,common_rebin*2,(0,250)),
    "h_2j1t_topPt": (   "2j1t top pt",250,common_rebin,(0,500)),
#    "h_2j1t_etajprime": ( "2j1t eta j' ",100,common_rebin,(0,4.7)),

    "h_2j1t_mtwcut_mtw": (     "2j1t mtw ",250,common_rebin,(0,500)),
    "h_2j1t_mtwcut_etajprime": ( "2j1t eta j' ",100,common_rebin,(0,4.7)),

    "h_2j1t_mtwcut_nextrajets": ( "2j1t nextrajets ",5,1,(-0.5,4.5)),
    "h_2j1t_mtwcut_topMass": ( "2j1t top mass",200,5,(100,500)),
    "h_2j1t_mtwcut_leadingextrajetpt": (     "2j1t leading extra jet pt ",100,1,(0,200)),
    "h_2j1t_mtwcut_leadingextrajetcsv": (     "2j1t leading extra jet csv discriminator ",100,common_rebin,(0,1.0)),
    "h_2j1t_mtwcut_topMassExtra": ( "2j1t top mass from the extra leading jet as b",200,5,(100,500)),

    "h_2j1t_mtwcut_sr_nextrajets": ( "2j1t nextrajets ",5,1,(-0.5,4.5)),
    "h_2j1t_mtwcut_sr_topMass": ( "2j1t top mass",200,rebin_sr,(100,500)),
    "h_2j1t_mtwcut_sr_leadingextrajetpt": (     "2j1t leading extra jet pt ",100,1,(0,200)),
    "h_2j1t_mtwcut_sr_leadingextrajetcsv": (     "2j1t leading extra jet csv discriminator ",100,rebin_sr,(0,1.0)),
    "h_2j1t_mtwcut_sr_topMassExtra": ( "2j1t top mass from the extra leading jet as b",200,rebin_sr,(100,500)),
    #"h_2j1t_mtwcut_sr_MuPt": (    "2j1t muon pt in sr with mtwcut",250,2,(0,500)),
    
    #"h_2j1t_muIso":("2j1t Muon Isolation",40,1,(0.0,1.0)),
    "h_2j1t_dR_lepjetpt40_1st":("2j1t #DeltaR(Lep,Jet) ",50,1,(0.0,5)),
    "h_2j1t_dPhi_lepjetpt40_1st":("2j1t #Delta#phi(Lep,Jet) ",50,1,(0.0,3.2)),
    "h_2j1t_dEta_lepjetpt40_1st":("2j1t #Delta#eta(Lep,Jet) ",50,2,(-3.2,3.2)),

    #3j1t
    "h_3j1t_bjetpt": (  "3j1t b jet pt ",250,common_rebin,(0,500)),
    "h_3j1t_MuPt": (    "muon pt ",250,2,(0,500)),
    "h_3j1t_mtwcut_MuPt": (    "3j1t muon pt ",250,2,(0,500)),
    "h_3j1t_MuEta": (   "muon eta ",100,common_rebin,(-2.1,2.1)),
    "h_3j1t_MuPhi": (   "muon phi ",100, common_rebin,(-3.2, 3.2)),
    "h_3j1t_MuE": (     "muon e ",100,common_rebin,(0,500)),
    "h_3j1t_mtw": (     "3j1t mtw ",250,common_rebin*2,(0,250)),
    #"h_3j1t_muIso":("3j2t Muon Isolation",40,1,(0.0,1.0)),
    "h_3j1t_dR_lepjetpt40_1st":("3j1t #DeltaR(Lep,Jet) ",50,1,(0.0,5)),
    "h_3j1t_dPhi_lepjetpt40_1st":("3j1t #Delta#phi(Lep,Jet) ",50,1,(0.0,3.2)),
    "h_3j1t_dEta_lepjetpt40_1st":("3j1t #Delta#eta(Lep,Jet) ",50,2,(-3.2,3.2)),
    
    
    #3j2t
    "h_3j2t_bjetpt": (     "3j2t b jet pt ",250,2,(0,500)),
    "h_3j2t_2ndbjetpt": (  "3j2t sub lead. b jet pt",250,2,(0,500)),
    "h_3j2t_mtwcut_MuPt": (    "3j1t muon pt ",250,2,(0,500)),
    "h_3j2t_MuPt": (       "3j2t muon pt ",250,2,(0,500)),
    "h_3j2t_MuEta": (      "3j2t muon eta ",100,2,(-2.1,2.1)),
    "h_3j2t_MuPhi": (      "3j2t muon phi ",100,2,(-3.2, 3.2)),
    "h_3j2t_MuE": (        "3j2t muon e ",250,2,(0,500)),
    #"h_3j2t_MuCharge": (   "3j2t muon charge ",2,1,(-1,1)),
    
    #"h_3j2t_mtw": (        "3j2t mtw ",250,common_rebin*2,(0,250)),
    "h_3j2t_ljetpt": (  "3j2t light jet pt ",100,2,(0,500)),
    "h_3j2t_ljetpt":(    "3j2t light jet pt ",100,2,(0,500)),
    "h_3j2t_ljeteta":(    "3j2t light jet eta ",100,1,(0,4.7)),
    "h_3j2t_etajprime":(      "3j2t eta j' ",100,2,(0.,4.7)),
    "h_3j2t_nextrajets":(   "3j2t number of extra jets ",5,1,(-0.5,4.5)),
    #"h_3j2t_muIso":("3j2t Muon Isolation",40,1,(0.0,1.0)),
    "h_3j2t_dR_lepjetpt40_1st":("3j2t #DeltaR(Lep,Jet) ",50,1,(0.0,5)),
    "h_3j2t_dPhi_lepjetpt40_1st":("3j2t #Delta#phi(Lep,Jet) ",50,1,(0.0,3.2)),
    "h_3j2t_dEta_lepjetpt40_1st":("3j2t #Delta#eta(Lep,Jet) ",50,2,(-3.2,3.2)),
}

store = [
    #"h_nJets", #SR shape
    "h_2j1t_topMass",
    #"h_nbJets", #SR shape
    "h_2j0t_jetpt40_leading", 
    "h_2j0t_jetpt40_subleading", 
    "h_2j0t_mtw", 
    "h_2j0t_mtwcut_MuPt", 
    "h_2j0t_MuPt", 
    #"h_2j0t_muIso",
    "h_2j0t_dR_lepjetpt40_1st",
    "h_2j0t_dPhi_lepjetpt40_1st",
    "h_2j0t_dEta_lepjetpt40_1st",


    "h_2j0t_mtwcut_jeteta40_leading",
    "h_2j0t_mtwcut_incljeteta",


    "h_2j0t_mtwcut_jetabseta40_leading",
    "h_2j0t_mtwcut_incljetabseta",
    "h_2j1t_bjetpt", 
    "h_2j1t_ljetpt", 
    "h_2j1t_ljeteta", 
    #"h_2j1t_nMu", 
    "h_2j1t_mtwcut_MuPt", 
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
    #"h_2j1t_muIso",
    "h_2j1t_dR_lepjetpt40_1st",
    "h_2j1t_dPhi_lepjetpt40_1st",
    "h_2j1t_dEta_lepjetpt40_1st",
    
    "h_3j1t_bjetpt", 
    "h_3j1t_MuPt", 
    "h_3j1t_mtwcut_MuPt", 
    "h_3j1t_MuEta", 
    "h_3j1t_MuPhi", 
    "h_3j1t_MuE",
    "h_3j1t_mtw", 
    #"h_3j1t_muIso",
    "h_3j1t_dR_lepjetpt40_1st",
    "h_3j1t_dPhi_lepjetpt40_1st",
    "h_3j1t_dEta_lepjetpt40_1st",
    
    "h_3j2t_bjetpt", 
    "h_3j2t_2ndbjetpt", 
    "h_3j2t_MuPt", 
    "h_3j2t_mtwcut_MuPt", 
    "h_3j2t_MuEta", 
    "h_3j2t_MuPhi", 
    "h_3j2t_MuE", 
    #"h_3j2t_MuCharge", 
    #"h_3j2t_mtw", 
#    "h_2j1t_etajprime",
    "h_3j2t_ljetpt",
    "h_3j2t_ljetpt",
    "h_3j2t_ljeteta",
    "h_3j2t_etajprime",
    "h_3j2t_nextrajets",

    #"h_3j2t_muIso",
    "h_3j2t_dR_lepjetpt40_1st",
    "h_3j2t_dPhi_lepjetpt40_1st",
    "h_3j2t_dEta_lepjetpt40_1st",
    
    "h_2j1t_mtwcut_mtw",
    "h_2j1t_mtwcut_topMass",
    "h_2j1t_mtwcut_etajprime",

    "h_2j1t_mtwcut_nextrajets",
    "h_2j1t_mtwcut_topMass",
    "h_2j1t_mtwcut_leadingextrajetpt",
    "h_2j1t_mtwcut_leadingextrajetcsv",
    "h_2j1t_mtwcut_topMassExtra",

    "h_2j1t_mtwcut_sr_nextrajets",
    "h_2j1t_mtwcut_sr_topMass",
    "h_2j1t_mtwcut_sr_leadingextrajetpt",
    "h_2j1t_mtwcut_sr_leadingextrajetcsv",
    "h_2j1t_mtwcut_sr_topMassExtra",
    #"h_2j1t_mtwcut_sr_MuPt",
    ]
