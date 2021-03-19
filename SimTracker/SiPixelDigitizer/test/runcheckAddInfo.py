#
import FWCore.ParameterSet.Config as cms

process = cms.Process("digiTest")

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.GeometrySimDB_cff')
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
# process.load("Configuration.StandardSequences.Services_cff")

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
#    input = cms.untracked.int32(1)
)

process.MessageLogger = cms.Service("MessageLogger",
    debugModules = cms.untracked.vstring('checkAddInfo'),
    destinations = cms.untracked.vstring('cout'),
#    destinations = cms.untracked.vstring("log","cout"),
    cout = cms.untracked.PSet(
        threshold = cms.untracked.string('ERROR')
    )
#    log = cms.untracked.PSet(
#        threshold = cms.untracked.string('DEBUG')
#    )
)

process.source = cms.Source("PoolSource",
    fileNames =  cms.untracked.vstring(
#    'file:/afs/cern.ch/work/d/dkotlins/public/MC/mu/pt100_71_pre7/digis/digis2_postls171.root'
    'file:/opt/sbg/cms/ui3_data1/ccollard/Pixels/test_CaroPreMix_PURun2_orig.root'
    )
)


process.TFileService = cms.Service("TFileService",
    fileName = cms.string('histo.root')
)

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2021_realistic', '')

  
process.analysis = cms.EDAnalyzer("checkAddInfo",
    Verbosity = cms.untracked.bool(True),
    src = cms.InputTag("simSiPixelDigis"),
#    src = cms.InputTag("mix"),
#    src = cms.InputTag("siPixelDigis"),
)

process.p = cms.Path(process.analysis)

