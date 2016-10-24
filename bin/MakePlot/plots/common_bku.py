met_range = (160,480, True)
# met_range = (160,520)

# title,    scale,  rebin, usrrng
settings = {
    'met'                    : ('MET', 10, 2, met_range),
    'dphi'                   : ('#Delta #phi (j_{1,...,6},E^{miss}_{T})', 10, None, None ),
    'dphi_6j'                : ('#Delta #phi (j_{1,...,6},E^{miss}_{T})', 10, None, None ),
    'dphi_6j_preS'           : ('#Delta #phi (j_{1,...,6},E^{miss}_{T})', 10, None, None ),
    'dphi_preS'              : ('#Delta #phi (j_{1,2},E^{miss}_{T})', 10, None, None ),
    'h_bdt_bjcsv'            : (None                  , 10, 2, None ),
    'h_bdt_dphij1b'          : (None                  , 10, 10, None ),
    'h_bdt_dphij2b'          : (None                  , 10, 10, None ),
    'h_bdt_jet1csv'          : (None                  , 10, 2, None ),
    'h_bdt_jet2csv'          : (None                  , 10, 2, None ),
    'h_bdt_prob'             : (None                  , 10, 2, None ),
    'h_bdt_qgl1'             : (None                  , 10, 2, None ),
    'h_bdt_qgl2'             : (None                  , 10, 2, None ),
    'h_jetQGL_preS'          : ('QGL discriminant'    , 10, 2, None ),
    'h_nbJets_preS'          : ('# b jets'            , 10, 1, (2, 6) ),
    'h_nJets_preS'           : ('# jets'              , 10, 1, (3, 11) ),
    'h_topMassPostFit'       : ('Mass_top - postfit'  , 10, 10, None ),
    'h_topMassPostFit_preS'  : ('Mass_top - postfit'  , 10, 10, None ),
    'h_topMassPreFit'        : ('Mass_top - prefit'   , 10, 10, None ),
    'h_topMassPreFit_preS'   : ('Mass_top - prefit'  , 10, 10, None ),
    'met_preS'               : ('MET'                 , 10, 2, met_range),
    'metFinal'               : ('MET'                 , 10, 2, met_range),
    'metFinal_Angular'       : ('MET'                 , 10, 2, met_range),
    'metFinal_Angular_tag'   : ('MET'                 , 10, 2, met_range),
    'metFinal_Angular_untag' : ('MET'                 , 10, 2, met_range),
    'metFinal_tag'           : ('MET'                 , 10, 2, met_range),
    'metFinal_untag'         : ('MET'                 , 10, 2, met_range),
    'mva_first'              : (None                  , 10, 10, None ),
}


store = [
    'met_preS',
    'metFinal',
    'metFinal_tag',
    'metFinal_untag',
    "metFinal_Angular",
    "metFinal_Angular_tag",
    "metFinal_Angular_untag",
]
