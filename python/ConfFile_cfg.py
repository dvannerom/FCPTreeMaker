import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

## Load the standard set of configuration modules
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

### MC global tag
#process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_2016_TrancheIV_v6'
### Data global tag
#process.GlobalTag.globaltag = '80X_dataRun2_2016SeptRepro_v6'
process.GlobalTag.globaltag = '80X_dataRun2_2016LegacyRepro_v4'

## Message Logger settings
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        ###### MC ######
        ### Signal ###
        #'file:ZprimeTo2Tauprime_MZ91p2_MTau1_QTau3Over3_pythia8_step2.root'
        #'file:ZprimeTo2Tauprime_MZ91p2_MTau1_QTau2Over3_pythia8_step2.root'
        #'file:ZprimeTo2Tauprime_MZ91p2_MTau1_QTau1Over3_pythia8_step2.root'
        #'file:ZprimeTo2Tauprime_MZ91p2_MTau1_QTau3Over3_ISR100_pythia8_step2.root'
        #'file:ZprimeTo2Tauprime_MZ91p2_MTau1_QTau2Over3_ISR100_pythia8_step2.root'
        #'file:ZprimeTo2Tauprime_MZ91p2_MTau1_QTau1Over3_ISR100_pythia8_step2.root'
        #'file:ZprimeTo2Tauprime_MZ91p2_MTau100_QTau3Over3_pythia8_step2.root'
        #'file:ZprimeTo2Tauprime_MZ91p2_MTau100_QTau2Over3_pythia8_step2.root'
        #'file:ZprimeTo2Tauprime_MZ91p2_MTau100_QTau1Over3_pythia8_step2.root'
        #'file:ZprimeTo2Tauprime_MZ91p2_MTau1_QTau3Over3_pythia8_step2_test.root'
        ### DYJetsToLL ###
        #'/store/mc/RunIISummer16DR80Premix/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/008FB50D-4DD2-E611-A238-0CC47A4C8E8A.root',
        #'/store/mc/RunIISummer16DR80Premix/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/00B114DD-60D2-E611-867A-0CC47A7C3612.root'
        #'/store/mc/RunIISummer16DR80Premix/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/110000/001AC973-60E2-E611-B768-001E67586A2F.root',
        #'/store/mc/RunIISummer16DR80Premix/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/110000/00376B09-8EE2-E611-A6F8-003048F5ADF6.root',
        #'/store/mc/RunIISummer16DR80Premix/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/110000/005AE20E-25E3-E611-90E1-FA163ECDDD3D.root'
        ### WJetsToLNu ###
        #'/store/mc/RunIISummer16DR80Premix/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/023963AD-CBBE-E611-AC12-D4AE526A048B.root',
        #'/store/mc/RunIISummer16DR80Premix/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/04444FF8-C3BE-E611-828A-0CC47A7D9966.root',
        #'/store/mc/RunIISummer16DR80Premix/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/049966D9-E9BE-E611-A935-70106F49CBD8.root',
        #'/store/mc/RunIISummer16DR80Premix/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/04CD3B59-D4BE-E611-A74B-00266CF3E0A4.root',
        #'/store/mc/RunIISummer16DR80Premix/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/061311C1-DBBE-E611-A958-002590AC4C49.root',
        #'/store/mc/RunIISummer16DR80Premix/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/066A6CA6-9EBE-E611-AD42-0CC47A7EED28.root',
        #'/store/mc/RunIISummer16DR80Premix/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/06754142-C3BE-E611-85E9-E41D2D08DFF0.root',
        #'/store/mc/RunIISummer16DR80Premix/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/0817EA07-43BE-E611-A42D-0CC47A4C8E1C.root',
        #'/store/mc/RunIISummer16DR80Premix/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/08A5BC66-C4BE-E611-821A-0025904C7DF0.root',
        #'/store/mc/RunIISummer16DR80Premix/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/08B98FD8-43BE-E611-A30C-0025905A6122.root'

        ###### Data ######
        '/store/data/Run2016H/SingleMuon/AOD/07Aug17-v1/910001/FAD5C0BF-368F-E711-858D-F04DA275BFA7.root '
    )
)

## Create output file
## Setup the service to make a ROOT TTree
process.TFileService = cms.Service("TFileService",
    fileName = cms.string("tuple.root"))
    #fileName = cms.string("MTau1_QTau3Over3_noSelection.root"))
    #fileName = cms.string("MTau1_QTau2Over3_noSelection.root"))
    #fileName = cms.string("MTau1_QTau1Over3_noSelection.root"))
    #fileName = cms.string("MTau1_QTau3Over3_ISR100_noSelection.root"))
    #fileName = cms.string("MTau1_QTau2Over3_ISR100_noSelection.root"))
    #fileName = cms.string("MTau1_QTau1Over3_ISR100_noSelection.root"))
    #fileName = cms.string("MTau100_QTau3Over3_noSelection.root"))
    #fileName = cms.string("MTau100_QTau2Over3_noSelection.root"))
    #fileName = cms.string("MTau100_QTau1Over3_noSelection.root"))
    #fileName = cms.string("DYJetsToLL_M-10To50.root"))
    #fileName = cms.string("DYJetsToLL_M-50.root"))
    #fileName = cms.string("WJetsToLNu.root"))

process.demo = cms.EDAnalyzer('FCPTreeMaker',
    isMC = cms.bool(False),
    triggerResults = cms.InputTag("TriggerResults", "", "HLT"),
    metFilters = cms.InputTag("TriggerResults", "", "RECO"),
    jets = cms.InputTag("ak4PFJetsCHS", "", "RECO"),
    muons = cms.InputTag("muons", "", "RECO"),
    tracks = cms.InputTag("generalTracks", "", "RECO"),
    primaryVertices = cms.InputTag("offlinePrimaryVertices", "", "RECO"),
    dedx = cms.InputTag("dedxHitInfo", "", "RECO"),
    met = cms.InputTag("pfMet", "", "RECO"),
    genparticles = cms.InputTag("genParticles", "", "HLT"),
)


process.p = cms.Path(process.demo)
