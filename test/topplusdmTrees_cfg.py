import commands, os
### *****************************************************************************************
### Usage:
###
### cmsRun topplusdmanaEDMntuples_cfg.py maxEvts=N sample="mySample/sample.root" version="1"7 outputLabel="myoutput"
###
### Default values for the options are set:
### maxEvts     = -1
### sample      = 'file:/scratch/decosa/ttDM/testSample/tlbsm_53x_v3_mc_10_1_qPV.root'
### outputLabel = 'analysisTTDM.root'
### *****************************************************************************************
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as opts



options = opts.VarParsing ('analysis')

options.register('maxEvts',
                 -1,# default value: process all events
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.int,
                 'Number of events to process')

SingleElTriggers = []
SingleMuTriggers = []
hadronTriggers = []


chan = "MET_Prompt"

chan = "TTbarDMJets_scalar_Mchi-50_Mphi-50"
#chan = "TTbarDMJets_scalar_Mchi-1_Mphi-50"
#chan = "TTbarDMJets_scalar_Mchi-10_Mphi-50"
chan = "DY"
chan = "WJ"
filedir= "/tmp/oiorio/"
cmd = "ls "+filedir+"/"+chan+"/"

status,ls_la = commands.getstatusoutput( cmd )
listFiles = ls_la.split(os.linesep)
files = []
#files = ["file:re-MiniAOD_17Jul/"+l for l in listFiles]
#files = ["file:"+filedir+"MET_Prompt/"+l for l in listFiles]
files = ["file:"+filedir+"/"+chan+"/"+l for l in listFiles]
options.register('sample',
                 #                 files,
                 #'root://xrootd.ba.infn.it//store/user/grauco/Ntuples_Fwk76v1_wPileUpJet/QCD_Pt-15to3000_TuneCUETP8M1_Flat_13TeV_pythia8/b2ganafw76x_QCD_Pt-15to3000_TuneCUETP8M1_Flat_13TeV_pythia8_JetMet_18Mar_v1/160318_104000/0000/B2GEDMNtuple_1.root',
                 #/store/user/grauco/Ntuples_Fwk76v1_wPileUpJet/QCD_Pt-15to3000_TuneCUETP8M1_Flat_13TeV_pythia8/b2ganafw76x_QCD_Pt-15to3000_TuneCUETP8M1_Flat_13TeV_pythia8_JetMet_18Mar_v1/160318_104000/0000/B2GEDMNtuple_1.root
                 #'root://xrootd.ba.infn.it//store/user/grauco/Ntuples_Fwk76v1_wPileUpJet/DoubleMuon/b2ganafw76x_DoubleMuon_02Mar_v3/160302_145236/0000/B2GEDMNtuple_101.root',
                 
                 #'root://xrootd.ba.infn.it//store/user/grauco/Ntuples_Fwk76v1/BprimeBToHB_M-1000_TuneCUETP8M1_13TeV-madgraph-pythia8/b2ganafw76x_BprimeBToHB_M-1000_TuneCUETP8M1_13TeV-madgraph-pythia8_23Feb_v1/160223_072932/0000/B2GEDMNtuple_1.root',
                 #'root://xrootd.ba.infn.it//store/user/algomez/QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_B2GAnaFW_v76x_v1p0/160222_145142/0000/B2GEDMNtuple_100.root',
                 #'root://xrootd.unl.edu//store/user/decosa/ttDM/CMSSW_7_4_15/MET/MET_Run2015D_miniAODv2_13Nov/151113_154306/0000/B2GEDMNtuple_1.root',
                 #'root://xrootd.unl.edu//store/user/dpinna/Samples_v2_Fw7_4_15/TTbarDMJets_scalar_Mchi-1_Mphi-10_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/DMtt_Mchi1Mphi10_v211Nov/151111_102628/0000/B2GEDMNtuple_1.root',
                 #'root://xrootd.unl.edu//store/user/oiorio/ttDM/samples/Nov2015/ST_t-channel_top_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1/EDMNTUPLE_12Nov/151113_183052/0000/B2GEDMNtuple_1.root',
#                 ['root://xrootd.unl.edu//store/user/dpinna/Samples_v2_Fw7_4_15/TTbarDMJets_scalar_Mchi-1_Mphi-10_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/DMtt_Mchi1Mphi10_v211Nov/1511#11_102628/0000/B2GEDMNtuple_1.root',
#                 'root://xrootd.unl.edu//store/user/dpinna/Samples_v2_Fw7_4_15/TTbarDMJets_scalar_Mchi-1_Mphi-10_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/DMtt_Mchi1Mphi10_v211Nov/15111#1_102628/0000/B2GEDMNtuple_2.root',
#                 'root://xrootd.unl.edu//store/user/dpinna/Samples_v2_Fw7_4_15/TTbarDMJets_scalar_Mchi-1_Mphi-10_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/DMtt_Mchi1Mphi10_v211Nov/15111#1_102628/0000/B2GEDMNtuple_3.root',],
                 #'root://xrootd.unl.edu//store/user/jkarancs/SusyAnalysis/B2GEdmNtuple/DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/B2GAnaFW_v74x_V8p4_RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v2/151111_094305/0000/B2GEDMNtuple_1.root',
#                 'root://xrootd.ba.infn.it//store/user/jkarancs/SusyAnalysis/B2GEdmNtuple/DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/B2GAnaFW_v74x_V8p4_RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v2/151111_094305/0000/B2GEDMNtuple_26.root',
#/store/user/decosa/ttDM/CMSSW_7_4_15/QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/QCD_HT300to500/151111_112957/0000/B2GEDMNtuple_1.root',
                 #'root://xrootd.ba.infn.it//store/user/decosa/ttDM/CMSSW_7_4_12/MET/MET_Run2015D_v3/150927_052306/0000/B2GEDMNtuple_1.root',

#                 'root://xrootd.ba.infn.it//store/user/decosa/ttDM/CMSSW_7_4_12/MET/MET_Run2015D_v5/151012_123712/0000/B2GEDMNtuple_104.root',
#                 'root://xrootd.ba.infn.it//store/user/dpinna/TTDM_Fwv7.4.x_v6.1_25ns/SingleElectron/SingleElectron_20Oct_v2/151020_082044/0000/B2GEDMNtuple_1.root',
#                 "file:/tmp/oiorio/B2G.root",
#                 "file:/tmp/oiorio/B2GWJetsNLO.root",
#                /user/oiorio/ttDM/samples/Nov2015/ST_t-channel_top_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1/EDMNTUPLE_12Nov/151 'root://xrootd.ba.infn.it//store/user/decosa/ttDM/CMSSW_7_4_12/SingleMuon/SingleMuon_Run2015D_v4/151012_123201/0000/B2GEDMNtuple_1.root',
#                 'root://xrootd.ba.infn.it//store/user/dpinna/TTDM_Fwv7.4.x_v6.1_25ns/TTbarDMJets_pseudoscalar_Mchi-1_Mphi-500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/psDMtt_Mchi1Mphi500_1Oct/150930_221122/0000/B2GEDMNtuple_1.root',
                 
#                 'root://xrootd.ba.infn.it//store/user/grauco/Ntuples_Fwk76v1_wPileUpJet/QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp/b2ganafw76x_QCD_Pt-15to7000_TuneCUETHS1_Flat_13TeV_herwigpp_02Mar_v3/160302_134653/0000/B2GEDMNtuple_1.root',
                 'file:/tmp/oiorio/B2GEDMNtuple.root',
                 #'file:B2GEDMNtuple_QCD.root',
                 #'file:B2GEDMNtuple_DoubleMuon.root',  
#                 'file:B2GEDMNtuple_test.root',
#                'file:/tmp/oiorio/BGEDMNtuple_TT.root',
#                'file:/tmp/oiorio/B2GEDMNtuple_QCD.root',
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.string,
                 'Sample to analyze')

options.register('version',
                 #'53',
                 '71',
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.string,
                 'ntuple version (53 or 71)')

options.register('outputLabel',
                 #                 'treesTTJetsNew_ttDMMCChangeJECs.root',
                 #'treesTTJetsNew_ttDMMCChangeJECs_DirectV4_MCCorr.root',
                 'treesTTJets.root',
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.string,
                 'Output label')

options.register('isData',
                 False,
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.bool,
                 'Is data?')

options.register('applyRes',
                 False,
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.bool,
                 'ApplyResiduals?')

options.register('addPartonInfo',
                 True,
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.bool,
                 'Add parton info??')



options.register('changeJECs',
                 False,
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.bool,
                 'Apply new JECs?')

options.register('LHE',
                 False,
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.bool,
                 'Keep LHEProducts')

options.register('lhesource',
                 'externalLHEProducer',
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.string,
                 'LHEProducts source')

options.register('channel',
                 'ttDM',
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.string,
                 'channel for weight evaluation'
                 )

options.parseArguments()

if(options.isData):options.LHE = False
options.LHE = False 
if(not options.isData): options.applyRes = False

l = ["singleTrigger"+str(s) for s in xrange(15)]
l = l + ["trigger2"+str(s) for s in xrange(15)]

SingleElTriggers = ["HLT_Ele27_eta2p1_WP75_Gsf"]
SingleElTriggers = SingleElTriggers + ["HLT_Ele27_eta2p1_WP75_Gsf_v"+str(s) for s in xrange(15)]

PhotonTriggers = ["HLT_Photon155"]
PhotonTriggers = PhotonTriggers + ["HLT_Photon155_v"+str(s) for s in xrange(15)]

SingleMuTriggers = ["HLT_IsoMu20","HLT_IsoTkMu20"]
SingleMuTriggers = SingleMuTriggers + ["HLT_IsoMu20_v"+str(s) for s in xrange(15)]
SingleMuTriggers = SingleMuTriggers + ["HLT_IsoTkMu20_v"+str(s) for s in xrange(15)]

hadronTriggers = ["HLT_PFMETNoMu120_NoiseCleaned_PFMHTNoMu120_IDTight","HLT_PFMETNoMu120_JetIdCleaned_PFMHTNoMu120_IDTight","HLT_PFMETNoMu90_NoiseCleaned_PFMHTNoMu90_IDTight","HLT_PFMETNoMu90_JetIdCleaned_PFMHTNoMu90_IDTight", "HLT_PFHT475_v2", "HLT_PFHT400_v2", "HLT_PFHT350_v3", "HLT_PFHT200_v2", "HLT_PFHT300_v2", "HLT_PFHT250_v2", "HLT_PFHT800_v2", "HLT_PFHT650_v3", "HLT_PFHT600_v3","HLT_AK8PFHT700_TrimR0p1PT0p03Mass50_v3", "HLT_Mu17_Mu8_v1", "HLT_ZeroBias_v2","HLT_Photon22_v2 ","HLT_Photon30_v3 ","HLT_Photon36_v3 ","HLT_Photon50_v3 ","HLT_Photon75_v3 ","HLT_Photon90_v3 ","HLT_Photon120_v3","HLT_Photon165_HE10_v3"]
#hadronTriggers = ["HLT_PFMETNoMu120_NoiseCleaned_PFMHTNoMu120_IDTight","HLT_PFMETNoMu120_JetIdCleaned_PFMHTNoMu120_IDTight"]

hadronTriggers = hadronTriggers + ["HLT_PFMETNoMu120_NoiseCleaned_PFMHTNoMu120_IDTight_v"+str(s) for s in xrange(15)]
hadronTriggers = hadronTriggers + ["HLT_PFMETNoMu120_JetIdCleaned_PFMHTNoMu120_IDTight_v"+str(s) for s in xrange(15)]
hadronTriggers = hadronTriggers + ["HLT_PFMETNoMu90_NoiseCleaned_PFMHTNoMu90_IDTight_v"+str(s) for s in xrange(15)]
hadronTriggers = hadronTriggers + ["HLT_PFMETNoMu90_JetIdCleaned_PFMHTNoMu90_IDTight_v"+str(s) for s in xrange(15)]

if(options.isData):

    SingleElTriggers = ["HLT_Ele27_eta2p1_WPLoose_Gsf", "HLT_Ele23_WPLoose_Gsf"]
    SingleElTriggers = SingleElTriggers + ["HLT_Ele27_eta2p1_WPLoose_Gsf_v"+str(s) for s in xrange(15)]
    SingleElTriggers = SingleElTriggers + ["HLT_Ele23_WPLoose_Gsf_v"+str(s) for s in xrange(15)]

    SingleMuTriggers = ["HLT_IsoMu20","HLT_IsoTkMu20"]
    SingleMuTriggers = SingleMuTriggers + ["HLT_IsoMu20_v"+str(s) for s in xrange(15)]
    SingleMuTriggers = SingleMuTriggers + ["HLT_IsoTkMu20_v"+str(s) for s in xrange(15)]

    hadronTriggers = ["HLT_PFMETNoMu120_NoiseCleaned_PFMHTNoMu120_IDTight","HLT_PFMETNoMu120_JetIdCleaned_PFMHTNoMu120_IDTight","HLT_PFMETNoMu90_JetIdCleaned_PFMHTNoMu90_IDTight","HLT_PFMETNoMu90_NoiseCleaned_PFMHTNoMu90_IDTight","HLT_PFMETNoMu90_NoiseCleaned_PFMHTNoMu90_NoID", "HLT_PFMETNoMu90_PFMHTNoMu90_IDTight", "HLT_PFHT475_v2", "HLT_PFHT400_v2", "HLT_PFHT350_v3", "HLT_PFHT200_v2", "HLT_PFHT300_v2", "HLT_PFHT250_v2", "HLT_PFHT650_v3", "HLT_PFHT600_v3", "HLT_PFHT800_v2", "HLT_AK8PFHT700_TrimR0p1PT0p03Mass50_v3", "HLT_Mu17_Mu8_v1", "HLT_ZeroBias_v2", "HLT_Photon22_v2 ","HLT_Photon30_v3 ","HLT_Photon36_v3 ","HLT_Photon50_v3 ","HLT_Photon75_v3 ","HLT_Photon90_v3 ","HLT_Photon120_v3", "HLT_Photon165_HE10_v3"]

    hadronTriggers = hadronTriggers+ ["HLT_PFMETNoMu120_NoiseCleaned_PFMHTNoMu120_IDTight_v"+str(s) for s in xrange(15)]
    hadronTriggers = hadronTriggers+ ["HLT_PFMETNoMu120_JetIdCleaned_PFMHTNoMu120_IDTight_v"+str(s) for s in xrange(15)]
    hadronTriggers = hadronTriggers+ ["HLT_PFMETNoMu90_JetIdCleaned_PFMHTNoMu90_IDTight_v"+str(s) for s in xrange(15)]
    hadronTriggers = hadronTriggers+ ["HLT_PFMETNoMu90_NoiseCleaned_PFMHTNoMu90_IDTight_v"+str(s) for s in xrange(15)]
    hadronTriggers = hadronTriggers+ ["HLT_PFMETNoMu90_NoiseCleaned_PFMHTNoMu90_NoID_v"+str(s) for s in xrange(15)]
    hadronTriggers = hadronTriggers+ ["HLT_PFMETNoMu90_PFMHTNoMu90_IDTight_v"+str(s) for s in xrange(15)]



process = cms.Process("ttDManalysisTrees")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.MessageLogger.categories.append('HLTrigReport')
### Output Report
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
### Number of maximum events to process
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvts) )
### Source file
process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring(
        options.sample
        )
)

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

from Configuration.AlCa.GlobalTag import GlobalTag as customiseGlobalTag
#process.GlobalTag = customiseGlobalTag(process.GlobalTag, globaltag = 'auto:run2_mc_50nsGRun')
#process.GlobalTag = customiseGlobalTag(process.GlobalTag, globaltag = '76X_mcRun2_asymptotic_v12')
process.GlobalTag = customiseGlobalTag(process.GlobalTag, globaltag = '80X_mcRun2_asymptotic_2016_miniAODv2')
#process.GlobalTag.connect   = 'frontier://FrontierProd/CMS_COND_31X_GLOBALTAG'
process.GlobalTag.pfnPrefix = cms.untracked.string('frontier://FrontierProd/')

### Rootplizer

process.TFileService = cms.Service("TFileService", fileName = cms.string(options.outputLabel))

process.load("Analysis.NAAnaFW.topplusdmedmRootTreeMaker_cff")
#process.load("ttDM.treeDumper.topplusdmedmRootTreeMaker_with_cat_cff")
#process.load("B2GAnaFW.B2GAnaFW.topplusdmedmRootTreeMaker_cff")

process.DMTreesDumper.channelInfo.SingleElTriggers=cms.vstring(SingleElTriggers)
process.DMTreesDumper.channelInfo.SingleMuTriggers=cms.vstring(SingleMuTriggers)
process.DMTreesDumper.channelInfo.hadronicTriggers=cms.vstring(hadronTriggers)



if options.addPartonInfo:
    if options.isData: #G
        process.DMTreesDumper.channelInfo.getPartonW=cms.untracked.bool(True)
        process.DMTreesDumper.channelInfo.getPartonTop=cms.untracked.bool(True)
        process.DMTreesDumper.channelInfo.doWReweighting=cms.untracked.bool(True)
        process.DMTreesDumper.channelInfo.doTopReweighting=cms.untracked.bool(True)

    if not options.isData: #G                                                                                                                            
        process.DMTreesDumper.channelInfo.getPartonW=cms.untracked.bool(True)
        process.DMTreesDumper.channelInfo.getPartonTop=cms.untracked.bool(True)
        process.DMTreesDumper.channelInfo.doWReweighting=cms.untracked.bool(True)
        process.DMTreesDumper.channelInfo.doTopReweighting=cms.untracked.bool(True)
#process.DMTreesDumper.lhes =cms.InputTag("externalLHEProducer")
process.DMTreesDumper.lhes =cms.InputTag(options.lhesource)
#process.DMTreesDumper.lhes =cms.InputTag(options.lhesource)
process.DMTreesDumper.changeJECs = cms.untracked.bool(options.changeJECs)
process.DMTreesDumper.isData = cms.untracked.bool(options.isData)#This adds the L2L3Residuals
process.DMTreesDumper.applyRes = cms.untracked.bool(options.applyRes)#This adds the L2L3Residuals

#G
process.DMTreesDumper.channelInfo.useLHE = cms.untracked.bool(False)
#G
process.DMTreesDumper.channelInfo.useLHEWeights = cms.untracked.bool(False)

if options.channel == "ttbar":
    process.DMTreesDumper.getPartonTop  = cms.untracked.bool(True)

if options.channel == "wzjets":
    print "channel is " + options.channel 
    process.DMTreesDumper.channelInfo.getPartonW  = cms.untracked.bool(True)
    process.DMTreesDumper.channelInfo.getParticleWZ  = cms.untracked.bool(True)

if options.isData:
    process.DMTreesDumper.channelInfo.useLHE = cms.untracked.bool(False)
    process.DMTreesDumper.channelInfo.useLHEWeights = cms.untracked.bool(False)

if not options.isData:
    process.DMTreesDumper.channelInfo.useLHE = cms.untracked.bool(False)
    process.DMTreesDumper.channelInfo.useLHEWeights = cms.untracked.bool(False)

if(options.isData): del process.DMTreesDumper.physicsObjects[-1]
process.analysisPath = cms.Path(
    process.DMTreesDumper
    )


#if(options.isData):
for p in process.DMTreesDumper.physicsObjects:
    if(p.prefix == cms.string("el")):
        p.variablesF.append(cms.InputTag("electrons","elvidVeto"))
        p.variablesF.append(cms.InputTag("electrons","elvidLoose"))
        p.variablesF.append(cms.InputTag("electrons","elvidMedium"))
        p.variablesF.append(cms.InputTag("electrons","elvidTight"))
            #print "yes"
        p.toSave.append("elvidVeto")
        p.toSave.append( "elvidLoose")
        p.toSave.append("elvidMedium")
        p.toSave.append("elvidTight")
            
