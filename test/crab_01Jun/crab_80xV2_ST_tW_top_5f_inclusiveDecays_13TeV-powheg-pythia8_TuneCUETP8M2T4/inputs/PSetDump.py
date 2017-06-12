import FWCore.ParameterSet.Config as cms

process = cms.Process("ttDManalysisTrees")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/o/oiorio/public/xWajid/synch/mc/B2GSynchMC.root')
)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

process.DMTreesDumper = cms.EDAnalyzer("DMAnalysisTreeMaker",
    HBHEFilter = cms.InputTag("HBHENoiseFilterResultProducer","HBHENoiseFilterResultRun2Loose"),
    HBHEIsoFilter = cms.InputTag("HBHENoiseFilterResultProducer","HBHEIsoNoiseFilterResult"),
    JECVersion = cms.string('Spring16_23Sep2016V2'),
    ak8jetSubjetIndex0 = cms.InputTag("jetsAK8CHS","jetAK8CHStopSubjetIndex0"),
    ak8jetSubjetIndex1 = cms.InputTag("jetsAK8CHS","jetAK8CHStopSubjetIndex1"),
    ak8jetSubjetIndex2 = cms.InputTag("jetsAK8CHS","jetAK8CHStopSubjetIndex2"),
    ak8jetSubjetIndex3 = cms.InputTag("jetsAK8CHS","jetAK8CHStopSubjetIndex3"),
    applyRes = cms.untracked.bool(False),
    boostedTopsLabel = cms.string('jetsAK8CHS'),
    boostedTopsSubjetsLabel = cms.string('subjetsAK8CHS'),
    changeJECs = cms.untracked.bool(True),
    channelInfo = cms.PSet(
        PhotonTriggers = cms.vstring(''),
        SingleElTriggers = cms.vstring('HLT_Ele32_eta2p1_WPTight_Gsf_v1', 
            'HLT_Ele32_eta2p1_WPTight_Gsf_v2', 
            'HLT_Ele32_eta2p1_WPTight_Gsf_v3', 
            'HLT_Ele32_eta2p1_WPTight_Gsf_v4', 
            'HLT_Ele32_eta2p1_WPTight_Gsf_v5', 
            'HLT_Ele32_eta2p1_WPTight_Gsf_v6', 
            'HLT_Ele32_eta2p1_WPTight_Gsf_v7', 
            'HLT_Ele32_eta2p1_WPTight_Gsf_v8', 
            'HLT_Ele32_eta2p1_WPTight_Gsf_v9', 
            'HLT_Ele27_eta2p1_WPTight_Gsf_v1', 
            'HLT_Ele27_eta2p1_WPTight_Gsf_v2', 
            'HLT_Ele27_eta2p1_WPTight_Gsf_v3', 
            'HLT_Ele27_eta2p1_WPTight_Gsf_v4', 
            'HLT_Ele27_eta2p1_WPTight_Gsf_v5', 
            'HLT_Ele27_eta2p1_WPTight_Gsf_v6', 
            'HLT_Ele27_eta2p1_WPTight_Gsf_v7', 
            'HLT_Ele27_eta2p1_WPTight_Gsf_v8', 
            'HLT_Ele27_eta2p1_WPTight_Gsf_v9'),
        SingleMuTriggers = cms.vstring('HLT_IsoMu22_v1', 
            'HLT_IsoMu22_v2', 
            'HLT_IsoMu22_v3', 
            'HLT_IsoMu22_v4', 
            'HLT_IsoMu24_v1', 
            'HLT_IsoMu24_v2', 
            'HLT_IsoMu24_v3', 
            'HLT_IsoMu24_v4', 
            'HLT_IsoTkMu24_v1', 
            'HLT_IsoTkMu24_v2', 
            'HLT_IsoTkMu24_v3', 
            'HLT_IsoTkMu24_v4'),
        addLHAPDFWeights = cms.untracked.bool(False),
        channel = cms.string('ttDM'),
        crossSection = cms.double(1),
        doTopDecayReshaping = cms.untracked.bool(False),
        doTopResweighting = cms.untracked.bool(False),
        doTopReweighting = cms.untracked.bool(True),
        doWReweighting = cms.untracked.bool(True),
        getPartonTop = cms.untracked.bool(True),
        getPartonW = cms.untracked.bool(True),
        hadronicTriggers = cms.vstring(''),
        maxWeights = cms.untracked.int32(111),
        metFilters = cms.vstring('Flag_CSCTightHaloFilter', 
            'Flag_goodVertices', 
            'Flag_eeBadScFilter'),
        originalEvents = cms.double(1),
        useLHE = cms.untracked.bool(False),
        useLHEWeights = cms.untracked.bool(False)
    ),
    cutOnTriggers = cms.untracked.bool(False),
    doPreselection = cms.untracked.bool(True),
    doResolvedTopHad = cms.untracked.bool(False),
    doResolvedTopSemiLep = cms.untracked.bool(False),
    eleLabel = cms.string('electrons'),
    eventLabel = cms.string('eventShapePFVars'),
    eventNumber = cms.InputTag("eventInfo","evtInfoEventNumber"),
    genJets = cms.InputTag("genJetsAK4"),
    genParticles = cms.InputTag("filteredPrunedGenParticles"),
    genprod = cms.InputTag("generator"),
    isData = cms.untracked.bool(False),
    jetKeysAK4CHS = cms.InputTag("jetKeysAK4CHS"),
    jetScanCuts = cms.vdouble(30),
    jetsLabel = cms.string('jetsAK4CHS'),
    lhes = cms.InputTag("externalLHEProducer"),
    lumiBlock = cms.InputTag("eventInfo","evtInfoLumiBlock"),
    metBits = cms.InputTag("METUserData","triggerBitTree"),
    metLabel = cms.string('metFull'),
    metNames = cms.InputTag("METUserData","triggerNameTree"),
    muLabel = cms.string('muons'),
    muonKeys = cms.InputTag("muonKeys"),
    partE = cms.InputTag("genPart","genPartE"),
    partEta = cms.InputTag("genPart","genPartEta"),
    partID = cms.InputTag("genPart","genPartID"),
    partMom0ID = cms.InputTag("genPart","genPartMom0ID"),
    partPhi = cms.InputTag("genPart","genPartPhi"),
    partPt = cms.InputTag("genPart","genPartPt"),
    partStatus = cms.InputTag("genPart","genPartStatus"),
    photonLabel = cms.string('photons'),
    physicsObjects = cms.VPSet(cms.PSet(
        categories = cms.vstring('Tight', 
            'TightAntiIso', 
            'Loose'),
        label = cms.string('muons'),
        maxInstances = cms.untracked.int32(10),
        prefix = cms.string('mu'),
        saveBaseVariables = cms.untracked.bool(False),
        saveNoCat = cms.untracked.bool(False),
        scanCuts = cms.vstring('Iso04_0p06_LE', 
            'Iso04_0p15_LE', 
            'Iso04_0p06_GE', 
            'Iso04_0p15_GE'),
        singleD = cms.VInputTag(),
        singleF = cms.VInputTag(),
        singleI = cms.VInputTag(),
        systCats = cms.vstring(''),
        toSave = cms.vstring('muE', 
            'muPt', 
            'muEta', 
            'muPhi', 
            'muIso04', 
            'muCharge', 
            'muIsTightMuon', 
            'muIsLooseMuon', 
            'allExtra'),
        variablesD = cms.VInputTag(),
        variablesF = cms.VInputTag(cms.InputTag("muons","muE"), cms.InputTag("muons","muPt"), cms.InputTag("muons","muEta"), cms.InputTag("muons","muPhi"), cms.InputTag("muons","muCharge"), 
            cms.InputTag("muons","muIsLooseMuon"), cms.InputTag("muons","muIsSoftMuon"), cms.InputTag("muons","muIsTightMuon"), cms.InputTag("muons","muGenMuonCharge"), cms.InputTag("muons","muGenMuonEta"), 
            cms.InputTag("muons","muGenMuonPt"), cms.InputTag("muons","muGenMuonE"), cms.InputTag("muons","muGenMuonPhi"), cms.InputTag("muons","muGlbTrkNormChi2"), cms.InputTag("muons","muInTrkNormChi2"), 
            cms.InputTag("muons","muIsGlobalMuon"), cms.InputTag("muons","muIsPFMuon"), cms.InputTag("muons","muIsTrackerMuon"), cms.InputTag("muons","muIso04"), cms.InputTag("muons","muNumberMatchedStations"), 
            cms.InputTag("muons","muNumberOfPixelLayers"), cms.InputTag("muons","muNumberOfValidTrackerHits"), cms.InputTag("muons","muNumberTrackerLayers"), cms.InputTag("muons","muNumberValidMuonHits"), cms.InputTag("muons","muNumberValidPixelHits"), 
            cms.InputTag("muons","muSumChargedHadronPt"), cms.InputTag("muons","muSumNeutralHadronPt"), cms.InputTag("muons","muSumPUPt"), cms.InputTag("muons","muSumPhotonPt")),
        variablesI = cms.VInputTag()
    ), 
        cms.PSet(
            categories = cms.vstring('Tight', 
                'TightAntiIso', 
                'Veto'),
            label = cms.string('electrons'),
            maxInstances = cms.untracked.int32(10),
            prefix = cms.string('el'),
            saveBaseVariables = cms.untracked.bool(False),
            saveNoCat = cms.untracked.bool(False),
            scanCuts = cms.vstring(),
            singleD = cms.VInputTag(),
            singleF = cms.VInputTag(),
            singleI = cms.VInputTag(),
            systCats = cms.vstring(''),
            toSave = cms.vstring('elE', 
                'elPt', 
                'elEta', 
                'elPhi', 
                'elIso03', 
                'elisTight', 
                'elvidTightnoiso', 
                'elCharge', 
                'elisMedium', 
                'elisLoose', 
                'elisVeto', 
                'elSCEta', 
                'allExtra', 
                'elDz', 
                'elDxy', 
                'elvidVeto', 
                'elvidLoose', 
                'elvidMedium', 
                'elvidTight'),
            variablesD = cms.VInputTag(),
            variablesF = cms.VInputTag(cms.InputTag("electrons","elE"), cms.InputTag("electrons","elPt"), cms.InputTag("electrons","elEta"), cms.InputTag("electrons","elPhi"), cms.InputTag("electrons","elCharge"), 
                cms.InputTag("electrons","elEta"), cms.InputTag("electrons","elSCEta"), cms.InputTag("electrons","elHoE"), cms.InputTag("electrons","elIso03"), cms.InputTag("electrons","eldEtaIn"), 
                cms.InputTag("electrons","eldPhiIn"), cms.InputTag("electrons","elmissHits"), cms.InputTag("electrons","elfull5x5siee"), cms.InputTag("electrons","elooEmooP"), cms.InputTag("electrons","elhasMatchedConVeto"), 
                cms.InputTag("electrons","elvidVeto"), cms.InputTag("electrons","elvidTight"), cms.InputTag("electrons","elvidTightnoiso"), cms.InputTag("electrons","elvidMedium"), cms.InputTag("electrons","elDz"), 
                cms.InputTag("electrons","elDxy"), cms.InputTag("electrons","elvidVeto"), cms.InputTag("electrons","elvidLoose"), cms.InputTag("electrons","elvidMedium"), cms.InputTag("electrons","elvidTight")),
            variablesI = cms.VInputTag()
        ), 
        cms.PSet(
            categories = cms.vstring('Tight'),
            label = cms.string('jetsAK4CHS'),
            maxInstances = cms.untracked.int32(20),
            prefix = cms.string('jetAK4CHS'),
            saveBaseVariables = cms.untracked.bool(False),
            saveNoCat = cms.untracked.bool(False),
            scanCuts = cms.vstring('CorrPt_20'),
            singleD = cms.VInputTag(),
            singleF = cms.VInputTag(),
            singleI = cms.VInputTag(),
            systCats = cms.vstring('JESUp', 
                'JESDown'),
            toSave = cms.vstring('jetAK4CHSE', 
                'jetAK4CHSPt', 
                'jetAK4CHSEta', 
                'jetAK4CHSPhi', 
                'jetAK4CHSGenJetPt', 
                'jetAK4CHSGenJetEta', 
                'jetAK4CHSCMVAv2', 
                'jetAK4CHSCSVv2', 
                'jetAK4CHSHadronFlavour', 
                'jetAK4CHSPartonFlavour', 
                'jetAK4CHSQGL', 
                'jetAK4CHSjecFactor0', 
                'jetAK4CHSpileupJetIdRMS', 
                'jetAK4CHSpileupJetIdbeta', 
                'jetAK4CHSpileupJetIdbetaClassic', 
                'jetAK4CHSpileupJetIdbetaStar', 
                'jetAK4CHSpileupJetIdbetaStarClassic', 
                'allExtra'),
            variablesD = cms.VInputTag(),
            variablesF = cms.VInputTag(cms.InputTag("jetsAK4CHS","jetAK4CHSE"), cms.InputTag("jetsAK4CHS","jetAK4CHSPt"), cms.InputTag("jetsAK4CHS","jetAK4CHSEta"), cms.InputTag("jetsAK4CHS","jetAK4CHSPhi"), cms.InputTag("jetsAK4CHS","jetAK4CHSPartonFlavour"), 
                cms.InputTag("jetsAK4CHS","jetAK4CHSPhi"), cms.InputTag("jetsAK4CHS","jetAK4CHSCSVv2"), cms.InputTag("jetsAK4CHS","jetAK4CHSCMVAv2"), cms.InputTag("jetsAK4CHS","jetAK4CHSCharge"), cms.InputTag("jetsAK4CHS","jetAK4CHSGenJetCharge"), 
                cms.InputTag("jetsAK4CHS","jetAK4CHSGenJetE"), cms.InputTag("jetsAK4CHS","jetAK4CHSGenJetEta"), cms.InputTag("jetsAK4CHS","jetAK4CHSGenJetPhi"), cms.InputTag("jetsAK4CHS","jetAK4CHSGenJetPt"), cms.InputTag("jetsAK4CHS","jetAK4CHSGenPartonCharge"), 
                cms.InputTag("jetsAK4CHS","jetAK4CHSGenPartonE"), cms.InputTag("jetsAK4CHS","jetAK4CHSGenPartonEta"), cms.InputTag("jetsAK4CHS","jetAK4CHSGenPartonPhi"), cms.InputTag("jetsAK4CHS","jetAK4CHSGenPartonPt"), cms.InputTag("jetsAK4CHS","jetAK4CHSHadronFlavour"), 
                cms.InputTag("jetsAK4CHS","jetAK4CHSneutralHadronEnergyFrac"), cms.InputTag("jetsAK4CHS","jetAK4CHSneutralEmEnergyFrac"), cms.InputTag("jetsAK4CHS","jetAK4CHSchargedEmEnergyFrac"), cms.InputTag("jetsAK4CHS","jetAK4CHSchargedHadronEnergyFrac"), cms.InputTag("jetsAK4CHS","jetAK4CHSneutralMultiplicity"), 
                cms.InputTag("jetsAK4CHS","jetAK4CHSchargedMultiplicity"), cms.InputTag("jetsAK4CHS","jetAK4CHSneutralEmEnergy"), cms.InputTag("jetsAK4CHS","jetAK4CHSjecFactor0"), cms.InputTag("jetsAK4CHS","jetAK4CHSjetArea")),
            variablesI = cms.VInputTag()
        ), 
        cms.PSet(
            categories = cms.vstring(),
            label = cms.string('metFull'),
            maxInstances = cms.untracked.int32(1),
            prefix = cms.string('metFull'),
            saveBaseVariables = cms.untracked.bool(True),
            scanCuts = cms.vstring(),
            singleD = cms.VInputTag(),
            singleF = cms.VInputTag(),
            singleI = cms.VInputTag(),
            systCats = cms.vstring(),
            toSave = cms.vstring(),
            variablesD = cms.VInputTag(),
            variablesF = cms.VInputTag(cms.InputTag("metFull","metFulluncorPt"), cms.InputTag("metFull","metFulluncorPhi"), cms.InputTag("metFull","metFullPt"), cms.InputTag("metFull","metFullPhi"), cms.InputTag("metFull","metFullPx"), 
                cms.InputTag("metFull","metFullPy")),
            variablesI = cms.VInputTag()
        )),
    runNumber = cms.InputTag("eventInfo","evtInfoRunNumber"),
    systematics = cms.vstring('noSyst'),
    triggerBits = cms.InputTag("TriggerUserData","triggerBitTree"),
    triggerNames = cms.InputTag("TriggerUserData","triggerNameTree"),
    triggerPrescales = cms.InputTag("TriggerUserData","triggerPrescaleTree"),
    useMETFilters = cms.untracked.bool(False),
    useTriggers = cms.untracked.bool(True),
    vertexChi2 = cms.InputTag("vertexInfo","chi"),
    vertexNdof = cms.InputTag("vertexInfo","ndof"),
    vertexRho = cms.InputTag("vertexInfo","rho"),
    vertexZ = cms.InputTag("vertexInfo","z")
)


process.analysisPath = cms.Path(process.DMTreesDumper)


process.MessageLogger = cms.Service("MessageLogger",
    FrameworkJobReport = cms.untracked.PSet(
        FwkJob = cms.untracked.PSet(
            limit = cms.untracked.int32(10000000),
            optionalPSet = cms.untracked.bool(True)
        ),
        default = cms.untracked.PSet(
            limit = cms.untracked.int32(0)
        ),
        optionalPSet = cms.untracked.bool(True)
    ),
    categories = cms.untracked.vstring('FwkJob', 
        'FwkReport', 
        'FwkSummary', 
        'Root_NoDictionary', 
        'HLTrigReport'),
    cerr = cms.untracked.PSet(
        FwkJob = cms.untracked.PSet(
            limit = cms.untracked.int32(0),
            optionalPSet = cms.untracked.bool(True)
        ),
        FwkReport = cms.untracked.PSet(
            limit = cms.untracked.int32(10000000),
            optionalPSet = cms.untracked.bool(True),
            reportEvery = cms.untracked.int32(1)
        ),
        FwkSummary = cms.untracked.PSet(
            limit = cms.untracked.int32(10000000),
            optionalPSet = cms.untracked.bool(True),
            reportEvery = cms.untracked.int32(1)
        ),
        INFO = cms.untracked.PSet(
            limit = cms.untracked.int32(0)
        ),
        Root_NoDictionary = cms.untracked.PSet(
            limit = cms.untracked.int32(0),
            optionalPSet = cms.untracked.bool(True)
        ),
        default = cms.untracked.PSet(
            limit = cms.untracked.int32(10000000)
        ),
        noTimeStamps = cms.untracked.bool(False),
        optionalPSet = cms.untracked.bool(True),
        threshold = cms.untracked.string('INFO')
    ),
    cerr_stats = cms.untracked.PSet(
        optionalPSet = cms.untracked.bool(True),
        output = cms.untracked.string('cerr'),
        threshold = cms.untracked.string('WARNING')
    ),
    cout = cms.untracked.PSet(
        placeholder = cms.untracked.bool(True)
    ),
    debugModules = cms.untracked.vstring(),
    debugs = cms.untracked.PSet(
        placeholder = cms.untracked.bool(True)
    ),
    default = cms.untracked.PSet(

    ),
    destinations = cms.untracked.vstring('warnings', 
        'errors', 
        'infos', 
        'debugs', 
        'cout', 
        'cerr'),
    errors = cms.untracked.PSet(
        placeholder = cms.untracked.bool(True)
    ),
    fwkJobReports = cms.untracked.vstring('FrameworkJobReport'),
    infos = cms.untracked.PSet(
        Root_NoDictionary = cms.untracked.PSet(
            limit = cms.untracked.int32(0),
            optionalPSet = cms.untracked.bool(True)
        ),
        optionalPSet = cms.untracked.bool(True),
        placeholder = cms.untracked.bool(True)
    ),
    statistics = cms.untracked.vstring('cerr_stats'),
    suppressDebug = cms.untracked.vstring(),
    suppressInfo = cms.untracked.vstring(),
    suppressWarning = cms.untracked.vstring(),
    warnings = cms.untracked.PSet(
        placeholder = cms.untracked.bool(True)
    )
)


process.TFileService = cms.Service("TFileService",
    fileName = cms.string('treesTTJets.root')
)


process.CastorDbProducer = cms.ESProducer("CastorDbProducer")


process.EcalLaserCorrectionService = cms.ESProducer("EcalLaserCorrectionService")


process.SiStripRecHitMatcherESProducer = cms.ESProducer("SiStripRecHitMatcherESProducer",
    ComponentName = cms.string('StandardMatcher'),
    NSigmaInside = cms.double(3.0),
    PreFilter = cms.bool(False)
)


process.StripCPEfromTrackAngleESProducer = cms.ESProducer("StripCPEESProducer",
    ComponentName = cms.string('StripCPEfromTrackAngle'),
    ComponentType = cms.string('StripCPEfromTrackAngle'),
    parameters = cms.PSet(
        mLC_P0 = cms.double(-0.326),
        mLC_P1 = cms.double(0.618),
        mLC_P2 = cms.double(0.3),
        mTEC_P0 = cms.double(-1.885),
        mTEC_P1 = cms.double(0.471),
        mTIB_P0 = cms.double(-0.742),
        mTIB_P1 = cms.double(0.202),
        mTID_P0 = cms.double(-1.427),
        mTID_P1 = cms.double(0.433),
        mTOB_P0 = cms.double(-1.026),
        mTOB_P1 = cms.double(0.253),
        maxChgOneMIP = cms.double(6000.0),
        useLegacyError = cms.bool(False)
    )
)


process.hcal_db_producer = cms.ESProducer("HcalDbProducer",
    dump = cms.untracked.vstring(''),
    file = cms.untracked.string('')
)


process.siPixelQualityESProducer = cms.ESProducer("SiPixelQualityESProducer",
    ListOfRecordToMerge = cms.VPSet(cms.PSet(
        record = cms.string('SiPixelQualityFromDbRcd'),
        tag = cms.string('')
    ), 
        cms.PSet(
            record = cms.string('SiPixelDetVOffRcd'),
            tag = cms.string('')
        ))
)


process.siStripBackPlaneCorrectionDepESProducer = cms.ESProducer("SiStripBackPlaneCorrectionDepESProducer",
    BackPlaneCorrectionDeconvMode = cms.PSet(
        label = cms.untracked.string('deconvolution'),
        record = cms.string('SiStripBackPlaneCorrectionRcd')
    ),
    BackPlaneCorrectionPeakMode = cms.PSet(
        label = cms.untracked.string('peak'),
        record = cms.string('SiStripBackPlaneCorrectionRcd')
    ),
    LatencyRecord = cms.PSet(
        label = cms.untracked.string(''),
        record = cms.string('SiStripLatencyRcd')
    )
)


process.siStripGainESProducer = cms.ESProducer("SiStripGainESProducer",
    APVGain = cms.VPSet(cms.PSet(
        Label = cms.untracked.string(''),
        NormalizationFactor = cms.untracked.double(1.0),
        Record = cms.string('SiStripApvGainRcd')
    ), 
        cms.PSet(
            Label = cms.untracked.string(''),
            NormalizationFactor = cms.untracked.double(1.0),
            Record = cms.string('SiStripApvGain2Rcd')
        )),
    AutomaticNormalization = cms.bool(False),
    appendToDataLabel = cms.string(''),
    printDebug = cms.untracked.bool(False)
)


process.siStripLorentzAngleDepESProducer = cms.ESProducer("SiStripLorentzAngleDepESProducer",
    LatencyRecord = cms.PSet(
        label = cms.untracked.string(''),
        record = cms.string('SiStripLatencyRcd')
    ),
    LorentzAngleDeconvMode = cms.PSet(
        label = cms.untracked.string('deconvolution'),
        record = cms.string('SiStripLorentzAngleRcd')
    ),
    LorentzAnglePeakMode = cms.PSet(
        label = cms.untracked.string('peak'),
        record = cms.string('SiStripLorentzAngleRcd')
    )
)


process.siStripQualityESProducer = cms.ESProducer("SiStripQualityESProducer",
    ListOfRecordToMerge = cms.VPSet(cms.PSet(
        record = cms.string('SiStripDetVOffRcd'),
        tag = cms.string('')
    ), 
        cms.PSet(
            record = cms.string('SiStripDetCablingRcd'),
            tag = cms.string('')
        ), 
        cms.PSet(
            record = cms.string('RunInfoRcd'),
            tag = cms.string('')
        ), 
        cms.PSet(
            record = cms.string('SiStripBadChannelRcd'),
            tag = cms.string('')
        ), 
        cms.PSet(
            record = cms.string('SiStripBadFiberRcd'),
            tag = cms.string('')
        ), 
        cms.PSet(
            record = cms.string('SiStripBadModuleRcd'),
            tag = cms.string('')
        ), 
        cms.PSet(
            record = cms.string('SiStripBadStripRcd'),
            tag = cms.string('')
        )),
    PrintDebugOutput = cms.bool(False),
    ReduceGranularity = cms.bool(False),
    ThresholdForReducedGranularity = cms.double(0.3),
    UseEmptyRunInfo = cms.bool(False),
    appendToDataLabel = cms.string('')
)


process.sistripconn = cms.ESProducer("SiStripConnectivity")


process.stripCPEESProducer = cms.ESProducer("StripCPEESProducer",
    ComponentName = cms.string('stripCPE'),
    ComponentType = cms.string('SimpleStripCPE'),
    parameters = cms.PSet(

    )
)


process.GlobalTag = cms.ESSource("PoolDBESSource",
    DBParameters = cms.PSet(
        authenticationPath = cms.untracked.string(''),
        authenticationSystem = cms.untracked.int32(0),
        messageLevel = cms.untracked.int32(0),
        security = cms.untracked.string('')
    ),
    DumpStat = cms.untracked.bool(False),
    ReconnectEachRun = cms.untracked.bool(False),
    RefreshAlways = cms.untracked.bool(False),
    RefreshEachRun = cms.untracked.bool(False),
    RefreshOpenIOVs = cms.untracked.bool(False),
    connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
    globaltag = cms.string('80X_mcRun2_asymptotic_2016_miniAODv2'),
    pfnPostfix = cms.untracked.string(''),
    pfnPrefix = cms.untracked.string('frontier://FrontierProd/'),
    snapshotTime = cms.string(''),
    toGet = cms.VPSet()
)


process.es_hardcode = cms.ESSource("HcalHardcodeCalibrations",
    GainWidthsForTrigPrims = cms.bool(False),
    HERecalibration = cms.bool(False),
    HEreCalibCutoff = cms.double(20.0),
    HFRecalibration = cms.bool(False),
    iLumi = cms.double(-1.0),
    testHFQIE10 = cms.bool(False),
    toGet = cms.untracked.vstring('GainWidths')
)


process.prefer("es_hardcode")

