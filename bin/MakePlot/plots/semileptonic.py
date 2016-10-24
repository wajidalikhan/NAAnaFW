 
from common import met_range
met_rebin = 4
# title,scale,rebin, usrrng
settings = {
    #'cutFlow'                : (None                           , 20, None, (-0.5,8.5) ),

    #N-1 plots
    'met'                    : ('E^{miss}_{T} [GeV]', 1, 1, met_range),  
    'muonmet'               : ('E^{miss}_{T} [GeV]'                          , 20, 1, met_range ), 
    'electronmet'            : ('E^{miss}_{T} [GeV]'                          , 20, 1, met_range ),

    'mt2w'                   : ('M^{W}_{T2} [GeV]'                         , 20, 10, (50,450) ),
    'mt'                     : ('M_{T} [GeV]'                        , 100, 5, (0,290) ),
    #'muonmt2w'        : ('M^{W}_{T2} [GeV]'                         , 20, 10, (50,450) ),
    #'muonmt'          : ('M_{T} [GeV]'                        , 20, 5, (0,290) ),
    #'electronmt2w'        : ('M^{W}_{T2} [GeV]'                         , 20, 5, (50,450) ),
    #'electronmt'          : ('M_{T} [GeV]'                        , 20, 2, (0,200) ),

    'muondphi'            : ('#Delta #phi (j_{1,2},E^{miss}_{T})', 20, 2, None ),
    'muondphi_6j'         : ('#Delta #phi (j_{1,...,6},E^{miss}_{T})', 20, 2, None ),
    #'electrondphi'            : ('#Delta #phi (j_{1,2},E^{miss}_{T})', 20, 1, None ),
    #'electrondphi_6j'         : ('#Delta #phi (j_{1,...,6},E^{miss}_{T})', 20, 1, None ),

    #SR shape
    'metFinal'               : ('E^{miss}_{T} [GeV]'                 , 20, met_rebin, met_range),
    #'metFinal_noMT2W'               : ('E^{miss}_{T} [GeV]'                 , 20, met_rebin, met_range),
    'muonmetFinal'          : ('E^{miss}_{T} [GeV]'                          , 20, met_rebin, met_range ),
    'electronmetFinal'       : ('E^{miss}_{T} [GeV]'                          , 20, met_rebin, met_range ),

    #pre-selection
    'met_preS'               : ('E^{miss}_{T} [GeV]'                 , 100, met_rebin, met_range),
    'muonmet_preS'          : ('E^{miss}_{T} [GeV]'                          , 100, met_rebin, met_range ),
    'electronmet_preS'       : ('E^{miss}_{T} [GeV]'                          , 100, met_rebin, met_range ),

    'mt2w_preS'              : ('M^{W}_{T2} [GeV]'                         , 100, 2, (50,450) ),
    'muonmt2w_preS'        : ('M^{W}_{T2} [GeV]'                         , 100, 4, (50,450)),
    'electronmt2w_preS'        : ('M^{W}_{T2} [GeV]'                         , 100, 4, (50,450) ),

    'mt_preaS'                : ('M_{T} [GeV]'                        , 100, 20, (0,290) ),
    'mt_EC_preS'                     : ('M_{T} [GeV] EC'                        , 100, 20, (0,290) ),
    'mt_BC_preS'                     : ('M_{T} [GeV] BC'                        , 100, 20, (0,290) ),
    'mt_EC_hadpreS'                     : ('M_{T} [GeV] EC'                        , 100, 20, (0,290) ),
    'mt_BC_hadpreS'                     : ('M_{T} [GeV] BC'                        , 100, 20, (0,290) ),
    'muonmt_preS'          : ('M_{T} [GeV]'                        , 100, 4, (0,290)),
    'electronmt_preS'          : ('M_{T} [GeV]'                        , 100, 4, (0,300) ),

    'h_muonnJets_preS'       : ('Number of tight jets'         , 100, 1, (3,11) ),
    'h_muonnbJets_preS'      : ('Number of tight b-jets'       , 100, 1, (1,6) ),
    'h_electronnJets_preS'       : ('Number of tight jets'       , 100, 1, (3,11) ),
    'h_electronnbJets_preS'      : ('Number of tight b-jets'     , 100, 1, (1,6) ),

    'muondphi_preS'       : ('#Delta #phi (j_{1,2},E^{miss}_{T})', 100, 2, None ),    
    'muondphi_6j_preS'    : ('#Delta #phi (j_{1,...,6},E^{miss}_{T})', 100, 2, None ),
    'electrondphi_preS'      : ('#Delta #phi (j_{1,2},E^{miss}_{T})', 100, 2, None ),
    'electrondphi_6j_preS'   : ('#Delta #phi (j_{1,...,6},E^{miss}_{T})', 100, 2, None ),

    'h_lep1Pt_preS'          : ('p_{T} leading lep'           , 100, 40, (30,190) ),
    'h_lep1Eta_preS'          : ('#eta leading lep'           , 100, 2, (-2.4,2.4) ),
    #'h_lep1Phi_preS'          : ('#phi leading lep'           , 100, 2, (-3.2,3.2) ),
    #'h_lep1E_preS'          : ('E leading lep'           , 100, 30, (30,390) ),

    'h_mu1Pt_preS'          : ('p_{T} leading muon'           , 100, 40, (30,190) ),
    'h_mu1Eta_preS'          : ('#eta leading muon'           , 100, 2, (-2.4,2.4) ),
    #'h_mu1Phi_preS'          : ('#phi leading muon'           , 100, 2, (-3.2,3.2) ),
    #'h_mu1E_preS'          : ('E leading muon'           , 100, 30, (30,390) ),

    'h_el1Pt_preS'          : ('p_{T} leading electron'           , 100, 40, (30,190) ),
    'h_el1Eta_preS'          : ('#eta leading electron'           , 100, 2, (-2.4,2.4) ),
    #'h_el1Phi_preS'          : ('#phi leading electron'           , 100, 2, (-3.2,3.2) ),
    #'h_el1E_preS'          : ('E leading electron'           , 100, 30, (30,350) ),

    'h_muonjet1Pt_preS'     : ('Leading jet Pt distribution'  , 100, 20, (30,430) ),
    'h_muonjet2Pt_preS'     : ('Sub-leading Jet Pt distribution' , 100, 20, (30,430)),
    'h_muonjet3Pt_preS'     : ('Third Jet Pt distribution'    , 100, 20, (30,260)  ),
    'h_electronjet1Pt_preS'     : ('Leading jet Pt distribution'  , 100, 20, (30,430) ),
    'h_electronjet2Pt_preS'     : ('Sub-leading Jet Pt distribution' , 100, 20, (30,430)),
    'h_electronjet3Pt_preS'     : ('Third Jet Pt distribution'    , 100, 20, (30,260)  ),

    #CR1
    'metFinal_2lep'          : ('E^{miss}_{T} [GeV]'                 , 200, met_rebin, met_range ),# Temp bug in cpp to fix
    'muonmetFinal_2lep'      : ('E^{miss}_{T} [GeV]'                          , 20, met_rebin, met_range ),
    'electronmetFinal_2lep'  : ('E^{miss}_{T} [GeV]'                          , 20, met_rebin, met_range ),
    'mixmetFinal_2lep'       : ('E^{miss}_{T} [GeV]'                          , 20, met_rebin, met_range ),

    'h_dphi_CRtt2_lep'       : ('#Delta #phi (j_{1,2},E^{miss}_{T})', 20, 2, None ),
    'muondphi_CR1'         : ('#Delta #phi (j_{1,2},E^{miss}_{T})', 20, 2, None ),
    'muondphi_6j_CR1'      : ('#Delta #phi (j_{1,...,6},E^{miss}_{T})', 20, 2, None ),
    'muonmt2w_CR1'         : ('M^{W}_{T2} [GeV]'                         , 20, 5, (50,450) ),
    'muonmt_CR1'           : ('M_{T} [GeV]'                        , 20, 10, (0,500) ),
    'electrondphi_CR1'     : ('#Delta #phi (j_{1,2},E^{miss}_{T})', 20, 2, None ),
    'electrondphi_6j_CR1'  : ('#Delta #phi (j_{1,...,6},E^{miss}_{T})', 20, 2, None ),
    'electronmt2w_CR1'     : ('M^{W}_{T2} [GeV]'                         , 20, 5, (50,450) ),
    'electronmt_CR1'       : ('M_{T} [GeV]'                        , 20, 10, (0,500) ),
    'mixdphi_CR1'          : ('#Delta #phi (j_{1,2},E^{miss}_{T})', 20, 2, None ),
    'mixdphi_6j_CR1'       : ('#Delta #phi (j_{1,...,6},E^{miss}_{T})', 20, 2, None ),
    'mixmt2w_CR1'          : ('M^{W}_{T2} [GeV]'                         , 20, 10, (50,450) ),
    'mixmt_CR1'            : ('M_{T} [GeV]'                        , 20, 20, (0,500) ),

    'h_lep1Pt_CR1'          : ('p_{T} leading lep'           , 20, 40, (30,190) ),
    'h_lep2Pt_CR1'          : ('p_{T} sub-leading lep'           , 20, 40, (30,190) ),
    'h_mu1Pt_CR1'          : ('p_{T} leading muon'           , 20, 40, (30,190) ),
    'h_mu2Pt_CR1'          : ('p_{T} sub-leading muon'           , 20, 40, (30,190) ),
    'h_el1Pt_CR1'          : ('p_{T} leading electron'           , 20, 40, (30,190) ),
    'h_el2Pt_CR1'          : ('p_{T} sub-leading electron'           , 20, 40, (30,190) ),

    'h_nJets_CRtt2_lep'      : ('Number of tight jets'         , 20, 1, (2,11)),
    'h_nbJets_CRtt2_lep'     : ('Number of tight b-jets'       , 20, 1, (1,6)),
    'h_muonnJets_CR1'      : ('Number of tight jets'         , 20, 1, (3,11) ),
    'h_muonnbJets_CR1'     : ('Number of tight b-jets'       , 20, 1, (1,6) ),
    'h_electronnJets_CR1'  : ('Number of tight jets'         , 20, 1, (3,11) ),
    'h_electronnbJets_CR1' : ('Number of tight b-jets'       , 20, 1, (1,6) ),
    'h_mixnJets_CR1'       : ('Number of tight jets'         , 20, 1, (3,11) ),
    'h_mixnbJets_CR1'      : ('Number of tight b-jets'       , 20, 1, (1,6) ),

    'h_jet1Pt_CRtt2_lep'     : ('Leading jet Pt distribution'  , 20, 20,  (30,430) ),
    'h_jet2Pt_CRtt2_lep'     : ('Sub-leading Jet Pt distribution' , 20, 20,  (30,430)),
    'h_jet3Pt_CRtt2_lep'     : ('Third Jet Pt distribution'    , 20, 20,  (30,260) ),
    'h_muonjet1Pt_CR1'     : ('Leading jet Pt distribution'  , 20, 50, (30,430) ),
    'h_muonjet2Pt_CR1'     : ('Sub-leading Jet Pt distribution' , 20, 50, (30,430)),
    'h_muonjet3Pt_CR1'     : ('Third Jet Pt distribution'    , 20, 50, (30,260)  ),
    'h_electronjet1Pt_CR1' : ('Leading jet Pt distribution'  , 20, 50, (30,430) ),
    'h_electronjet2Pt_CR1' : ('Sub-leading Jet Pt distribution' , 20, 50, (30,430)),
    'h_electronjet3Pt_CR1' : ('Third Jet Pt distribution'    , 20, 50, (30,260)  ),
    'h_mixjet1Pt_CR1'      : ('Leading jet Pt distribution'  , 20, 20,  (30,430) ),
    'h_mixjet2Pt_CR1'      : ('Sub-leading Jet Pt distribution' , 20, 20,  (30,430)),
    'h_mixjet3Pt_CR1'      : ('Third Jet Pt distribution'    , 20, 20,  (30,260)  ),

    #Z window
    #'metFinal_2lep_Z_nobtag' : ('E^{miss}_{T} [GeV]'                          , 20, met_rebin, met_range ),#
    #'muonmetFinal_2lep_Z_nobtag' : ('E^{miss}_{T} [GeV]'                          , 20, met_rebin, met_range ),#
    #'electronmetFinal_2lep_Z_nobtag' : ('E^{miss}_{T} [GeV]'                          , 20, met_rebin, met_range ),#
    #'mixmetFinal_2lep_Z_nobtag' : ('E^{miss}_{T} [GeV]'                          , 20, met_rebin, met_range ),#

    'h_nJets_CRZ_lep'        : ('Number of tight jets'         , 20, 1, (2,11) ),
    'h_nbJets_CRZ_lep'       : ('Number of tight b-jets'       , 20, 1, (1,6) ),
    'h_dphi_CRZ_lep'         : ('#Delta #phi (j_{1,2},E^{miss}_{T})', 20, 2, None ),

    'h_jet1Pt_CRZ_lep'       : ('Leading jet Pt distribution'  , 20, 20,  (30,430) ),
    'h_jet2Pt_CRZ_lep'       : ('Sub-leading Jet Pt distribution' , 20, 20,  (30,330) ),
    'h_jet3Pt_CRZ_lep'       : ('Third Jet Pt distribution'    , 20, 20,  (30,260)  ),


    #CR2
    'metFinal_met_0btag'     : ('E^{miss}_{T} [GeV]'                          , 100, met_rebin, met_range ),#
    'muonmetFinal_met_0btag' : ('E^{miss}_{T} [GeV]'                          , 20, met_rebin, met_range ),
    'electronmetFinal_met_0btag'     : ('E^{miss}_{T} [GeV]'                  , 20, met_rebin, met_range ),

    'h_mt2w_CRwj_lep'      : ('M^{W}_{T2} [GeV]'                         , 20, 10, (50,450) ),
    'h_mt_CRwj_lep'        : ('M_{T} [GeV]'                        , 20, 20, (160,500) ),

    'h_dphi_CRwj_lep'        : ('#Delta #phi (j_{1,2},E^{miss}_{T})', 20, 2, None ),
    'dphi_6j_CRwj_lep'        : ('#Delta #phi (j_{1,...,6},E^{miss}_{T})', 20, 2, None ),
    'muondphi_CR2'         : ('#Delta #phi (j_{1,2},E^{miss}_{T})', 20, 2, None ),
    'muondphi_6j_CR2'      : ('#Delta #phi (j_{1,...,6},E^{miss}_{T})', 20, 2, None ),
    'electrondphi_CR2'     : ('#Delta #phi (j_{1,2},E^{miss}_{T})', 20, 2, None ),
    'electrondphi_6j_CR2'  : ('#Delta #phi (j_{1,...,6},E^{miss}_{T})', 20, 2, None ),

    'h_lepPt_CRwj_lep'          : ('p_{T} lepton'           , 20, 30, (30,190) ),
    'h_mu1Pt_CR2'          : ('p_{T} muon'           , 20, 30, (30,190) ),
    'h_el1Pt_CR2'          : ('p_{T} electron'           , 20, 30, (30,190) ),

    'h_nJets_CRwj_lep'       : ('Number of tight jets'         , 20, 1, (2,11) ),
    'h_nbJets_CRwj_lep'      : ('Number of tight b-jets'       , 20, 1, None ),
    'h_muonnJets_CR2'      : ('Number of tight jets'         , 20, 1, (3,11) ),
    'h_muonnbJets_CR2'     : ('Number of tight b-jets'       , 20, 1, (1,6) ),
    'h_electronnJets_CR2'  : ('Number of tight jets'         , 20, 1, (3,11) ),
    'h_electronnbJets_CR2' : ('Number of tight b-jets'       , 20, 1, (1,6) ),
    
    'h_jet1Pt_CRwj_lep'      : ('Leading jet Pt distribution'  , 20, 20,  (30,430) ),
    'h_jet2Pt_CRwj_lep'      : ('Sub-leading Jet Pt distribution' , 20, 20,  (30,390) ),
    'h_jet3Pt_CRwj_lep'      : ('Third Jet Pt distribution'    , 20, 20, (30,260)  ),
    'h_muonjet1Pt_CR2'     : ('Leading jet Pt distribution'  , 20, 50, (30,430) ),
    'h_muonjet2Pt_CR2'     : ('Sub-leading Jet Pt distribution' , 20, 50, (30,430)),
    'h_muonjet3Pt_CR2'     : ('Third Jet Pt distribution'    , 20, 50, (30,260)  ),
    'h_electronjet1Pt_CR2' : ('Leading jet Pt distribution'  , 20, 50, (30,430) ),
    'h_electronjet2Pt_CR2' : ('Sub-leading Jet Pt distribution' , 20, 50, (30,430)),
    'h_electronjet3Pt_CR2' : ('Third Jet Pt distribution'    , 20, 50, (30,260)  ),
    }


#plot to store in histo folder
store = [
    'metFinal',
    'metFinal_2lep', #CR1
    'metFinal_met_0btag' #CR2
    #'metFinal_noMT2W', #SR shape no mt2w cut
    #    'h_topSLCosBM',
    #    'h_topSLCosLM',
    #    'h_topSLCosTM',
    #    'h_topSLMass',
    #    'h_topSLPt',
    #    'h_toSLCosLBM',
    #    'metFinal_2lep',
    #    'metFinal_2lep_Full',
    #    'metFinal_2lep_noZ',
    #    'metFinal_2lep_Z_Nobtag',
    #    'metFinal_cosTM00',
    #    'metFinal_cosTM02',
    #    'metFinal_cosTM04',
    #    'metFinal_cosTMp02',
    ##    'metFinal_cosTMp04',
    #    'metFinal_sl3Jets',
    #    'metFinal_sl4Jets',
    #    'metFinal_tag_CR1',
    #    'metFinal_tag_CR2',
    #    'metFinal_tag_CR3',
    #    'metFinal_tag_CR4',
    #    'metFinal_untag_CR1',
    #    'metFinal_untag_CR2',
    #    'metFinal_untag_CR3',
    #    'metFinal_untag_CR4'
 ]
