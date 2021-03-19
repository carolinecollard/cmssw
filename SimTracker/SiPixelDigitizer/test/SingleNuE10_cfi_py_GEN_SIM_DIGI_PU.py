# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: SingleNuE10_cfi.py --conditions auto:phase1_2021_realistic --pileup_input das:/RelValMinBias_14TeV/CMSSW_10_6_1-106X_mcRun3_2021_realistic_v1_rsb-v1/GEN-SIM -n 10 --era Run3 --eventcontent PREMIX --procModifiers premix_stage1 --relval 9000,50 -s GEN,SIM,DIGI:pdigi_valid --beamspot Run3RoundOptics25ns13TeVLowSigmaZ --datatier PREMIX --pileup Run3_Flat55To75_PoissonOOTPU --geometry DB:Extended
import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Run3_cff import Run3
from Configuration.ProcessModifiers.premix_stage1_cff import premix_stage1

process = cms.Process('DIGI',Run3,premix_stage1)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
#process.load('SimGeneral.MixingModule.mix_Run3_Flat55To75_PoissonOOTPU_cfi')
# 1 seul PU par event
#process.load('SimGeneral.MixingModule.mix_FIX_average_cfi')
# run2 scenario
process.load('SimGeneral.MixingModule.mix_2018_25ns_UltraLegacy_PoissonOOTPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.GeometrySimDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedRun3RoundOptics25ns13TeVLowSigmaZ_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10),
    output = cms.optional.untracked.allowed(cms.int32,cms.PSet)
)

# Input source
process.source = cms.Source("EmptySource")

process.options = cms.untracked.PSet(
    FailPath = cms.untracked.vstring(),
    IgnoreCompletely = cms.untracked.vstring(),
    Rethrow = cms.untracked.vstring(),
    SkipEvent = cms.untracked.vstring(),
    allowUnscheduled = cms.obsolete.untracked.bool,
    canDeleteEarly = cms.untracked.vstring(),
    emptyRunLumiMode = cms.obsolete.untracked.string,
    eventSetup = cms.untracked.PSet(
        forceNumberOfConcurrentIOVs = cms.untracked.PSet(

        ),
        numberOfConcurrentIOVs = cms.untracked.uint32(1)
    ),
    fileMode = cms.untracked.string('FULLMERGE'),
    forceEventSetupCacheClearOnNewRun = cms.untracked.bool(False),
    makeTriggerResults = cms.obsolete.untracked.bool,
    numberOfConcurrentLuminosityBlocks = cms.untracked.uint32(1),
    numberOfConcurrentRuns = cms.untracked.uint32(1),
    numberOfStreams = cms.untracked.uint32(0),
    numberOfThreads = cms.untracked.uint32(1),
    printDependencies = cms.untracked.bool(False),
    sizeOfStackForThreadsInKB = cms.optional.untracked.uint32,
    throwIfIllegalParameter = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(False)
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('SingleNuE10_cfi.py nevts:10'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.PREMIXoutput = cms.OutputModule("PoolOutputModule",
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    ),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('PREMIX'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('/opt/sbg/cms/ui3_data1/ccollard/Pixels/test_CaroPreMix.root'),
    outputCommands = process.PREMIXEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition

# Other statements
process.mix.input.fileNames = cms.untracked.vstring(['/store/relval/CMSSW_10_6_1/RelValMinBias_14TeV/GEN-SIM/106X_mcRun3_2021_realistic_v1_rsb-v1/10000/0A6E8B7E-490A-774F-B080-EA948DE7A10B.root', '/store/relval/CMSSW_10_6_1/RelValMinBias_14TeV/GEN-SIM/106X_mcRun3_2021_realistic_v1_rsb-v1/10000/2F4B471B-CB4E-1840-A774-4AE7E4220D9D.root', '/store/relval/CMSSW_10_6_1/RelValMinBias_14TeV/GEN-SIM/106X_mcRun3_2021_realistic_v1_rsb-v1/10000/40A5955D-1071-F64B-984D-F7EAA72A6C77.root', '/store/relval/CMSSW_10_6_1/RelValMinBias_14TeV/GEN-SIM/106X_mcRun3_2021_realistic_v1_rsb-v1/10000/70672451-102A-8243-ACF0-C4F545C049AB.root', '/store/relval/CMSSW_10_6_1/RelValMinBias_14TeV/GEN-SIM/106X_mcRun3_2021_realistic_v1_rsb-v1/10000/97553220-ED91-4E48-B59F-ED44053F8621.root', '/store/relval/CMSSW_10_6_1/RelValMinBias_14TeV/GEN-SIM/106X_mcRun3_2021_realistic_v1_rsb-v1/10000/A810739C-94E2-E44C-8E78-EDFF5775947D.root'])
process.XMLFromDBSource.label = cms.string("Extended")
process.genstepfilter.triggerConditions=cms.vstring("generation_step")
process.mix.digitizers = cms.PSet(process.theDigitizersValid)
process.mix.digitizers.pixel.UseReweighting = cms.bool(True)
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2021_realistic', '')

process.generator = cms.EDProducer("FlatRandomEGunProducer",
    AddAntiParticle = cms.bool(False),
    PGunParameters = cms.PSet(
        MaxE = cms.double(10.01),
        MaxEta = cms.double(2.5),
        MaxPhi = cms.double(3.14159265359),
        MinE = cms.double(9.99),
        MinEta = cms.double(-2.5),
        MinPhi = cms.double(-3.14159265359),
        PartID = cms.vint32(12)
    ),
    Verbosity = cms.untracked.int32(0),
    firstRun = cms.untracked.uint32(1),
    psethack = cms.string('single Nu E 10')
)


# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.simulation_step = cms.Path(process.psim)
process.digitisation_step = cms.Path(process.pdigi_valid)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.PREMIXoutput_step = cms.EndPath(process.PREMIXoutput)

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.simulation_step,process.digitisation_step,process.endjob_step,process.PREMIXoutput_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)
# filter all path with the production filter sequence
for path in process.paths:
	getattr(process,path).insert(0, process.generator)

#process.load("FWCore.MessageService.MessageLogger_cfi")
#process.MessageLogger.cerr.threshold = "DEBUG"
#process.MessageLogger = cms.Service(
#    "MessageLogger",
#    destinations = cms.untracked.vstring(
#        'detailedInfo',
#         'critical'
#         ),
#    detailedInfo = cms.untracked.PSet(
#        threshold = cms.untracked.string('DEBUG')
#         ),
#    debugModules = cms.untracked.vstring( 
#        '*'
##        'PreMixingSiPixelWorker',
##        'SiPixelDigitizer' 
#        )
#    )

# Customisation from command line

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
