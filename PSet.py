import FWCore.ParameterSet.Config as cms

process = cms.Process("TESTPROD")

process.maxEvents = cms.untracked.PSet(
     input = cms.untracked.int32(100)
)

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 50000

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
                            )
)
process.source.fileNames.append(
    'file:///home/mbluj/work/data/NanoAOD/80X_with941/VBFHToTauTau_M-125_80X/v3/test80X_NANO_1.root'
)
#process.p = cms.Path()

#process.schedule = cms.Schedule(process.p)
