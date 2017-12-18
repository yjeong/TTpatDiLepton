import FWCore.ParameterSet.Config as cms
process = cms.Process("TtbarDiLeptonAnalyzer")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )
process.options.allowUnscheduled = cms.untracked.bool(True)

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(
	'root://cms-xrdr.sdfarm.kr:1094///xrd/store/group/CAT/v8-0-8/ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/171126_165300/0000/catTuple_1.root',#SingleTbar_tW 
	#'root://cms-xrdr.sdfarm.kr:1094///xrd/store/group/CAT/v8-0-8/ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/171126_165041/0000/catTuple_1.root',#SingleTop_tW
 )
)
#from CATTools.Validation.commonTestInput_cff import commonTestCATTuples
#process.source.fileNames = commonTestCATTuples["bkg"]

process.load("CATTools.CatAnalyzer.ttll.ttbarDileptonKinSolutionAlgos_cff")
process.load("CATTools.CatAnalyzer.filters_cff")
process.load("CATTools.CatAnalyzer.topPtWeightProducer_cfi")
process.load("CATTools.CatAnalyzer.flatGenWeights_cfi")
#process.load("CATTools.CatProducer.mcTruthTop.partonTop_cfi")#except when running SingleTop, and SingleTbar
from CATTools.CatAnalyzer.leptonSF_cff import *

process.ttbarDileptonKinAlgoPSetDESYSmeared.inputTemplatePath = cms.string("CATTools/CatAnalyzer/data/desyKinRecoInput.root")
process.ttbarDileptonKinAlgoPSetDESYSmeared.maxLBMass = cms.double(180)
process.ttbarDileptonKinAlgoPSetDESYSmeared.mTopInput = cms.double(172.5)
process.ttbarDileptonKinAlgoPSetDESYSmearedPseudoTop = process.ttbarDileptonKinAlgoPSetDESYSmeared.clone()
process.ttbarDileptonKinAlgoPSetDESYSmearedPseudoTop.inputTemplatePath = cms.string("CATTools/CatAnalyzer/data/KoreaKinRecoInput_pseudo.root")
process.ttbarDileptonKinAlgoPSetDESYSmearedPseudoTop.maxLBMass = cms.double(360)
process.ttbarDileptonKinAlgoPSetDESYSmearedPseudoTop.mTopInput = cms.double(172.5)

process.cattree = cms.EDAnalyzer("TtbarDiLeptonAnalyzer",
    recoFilters = cms.InputTag("filterRECO"),
    nGoodVertex = cms.InputTag("catVertex","nGoodPV"),
    lumiSelection = cms.InputTag("lumiMask"),
    genweight = cms.InputTag("flatGenWeights"),
    pdfweights = cms.InputTag("flatGenWeights","pdf"),
    scaleupweights = cms.InputTag("flatGenWeights","scaleup"),
    scaledownweights = cms.InputTag("flatGenWeights","scaledown"),
    topPtWeight = cms.InputTag("topPtWeight"),
    puweight = cms.InputTag("pileupWeight"),
    puweight_up = cms.InputTag("pileupWeight","up"),
    puweight_dn = cms.InputTag("pileupWeight","dn"),
    trigMUEL = cms.InputTag("filterTrigMUEL"),
    trigMUMU = cms.InputTag("filterTrigMUMU"),
    trigELEL = cms.InputTag("filterTrigELEL"),

    jets = cms.InputTag("catJets"),
    mets = cms.InputTag("catMETs"),
    vertices = cms.InputTag("catVertex"),
    muon = cms.PSet(
        src = cms.InputTag("catMuons"),
        effSF = muonSFTight,
    ),
    electron = cms.PSet(
        src = cms.InputTag("catElectrons"),
        effSF = electronSFCutBasedIDMediumWP,#electronSFWP90,
    ),
    mcLabel = cms.InputTag("prunedGenParticles"),

    partonTop_channel = cms.InputTag("partonTop","channel"),
    partonTop_modes = cms.InputTag("partonTop", "modes"),
    partonTop_genParticles = cms.InputTag("partonTop"),

    particleLevel = cms.InputTag("particleLevel"),

## Dstar begin
    d0s    = cms.InputTag("catDstars","D0Cand"),
    dstars = cms.InputTag("catDstars","DstarCand"),
    Jpsis  = cms.InputTag("catDstars","JpsiCand"),
    matchingDeltaR = cms.double(0.15),
## Dstar end

    #solver = process.ttbarDileptonKinAlgoPSetCMSKin,
    solver = process.ttbarDileptonKinAlgoPSetDESYSmeared,
    solverPseudoTop = process.ttbarDileptonKinAlgoPSetDESYSmearedPseudoTop,
    #solver = process.ttbarDileptonKinAlgoPSetDESYMassLoop,
)
#process.cattree.solver.tMassStep = 1
if cms.string('DESYSmeared') == process.cattree.solver.algo:
    process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
        cattree = cms.PSet(
            initialSeed = cms.untracked.uint32(123456),
            engineName = cms.untracked.string('TRandom3')
        )
    )

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("cattree.root"))

process.p = cms.Path(process.cattree)
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
