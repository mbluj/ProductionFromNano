#!/usr/bin/env python

import os

from ROOT import gSystem, TChain, TSystem, TFile
from PSet import process

#doSvFit = True
doSvFit = False
applyRecoil=True
#applyRecoil=False
if doSvFit :
    print "Run with SVFit computation"
if applyRecoil :
    print "Apply MET recoil corrections"

#Some system have problem runnig compilation (missing glibc-static library?).
#First we try to compile, and only ther we start time consuming cmssw
status = gSystem.CompileMacro('HTTEvent.cxx')
status *= gSystem.CompileMacro('NanoEventsSkeleton.C')
gSystem.Load('$CMSSW_BASE/lib/$SCRAM_ARCH/libTauAnalysisClassicSVfit.so')
gSystem.Load('$CMSSW_BASE/lib/$SCRAM_ARCH/libTauAnalysisSVfitTF.so')
gSystem.Load('$CMSSW_BASE/lib/$SCRAM_ARCH/libHTT-utilitiesRecoilCorrections.so')
status *= gSystem.CompileMacro('HTauTauTreeFromNanoBase.C')
status *= gSystem.CompileMacro('HMuTauhTreeFromNano.C')
status *= gSystem.CompileMacro('HTauhTauhTreeFromNano.C')

print "Compilation status: ",status
if status==0:
    exit(-1)

from ROOT import HMuTauhTreeFromNano, HTauhTauhTreeFromNano
fileNames = [
    "0E6F4B78-CC12-E811-B37D-FA163EA12C78.root",
    "50BE09DD-CC12-E811-869D-F04DA27542B9.root",
    "844BE355-CD12-E811-8871-FA163ED9B872.root",
    "DEBF5F61-CC12-E811-B47A-0CC47AA9943A.root",
    "5A038C2A-CC12-E811-B729-7845C4FC3B8D.root",
]
for name in fileNames:
    aFile = "file:///home/mbluj/work/data/NanoAOD/80X_with944/VBFHToTauTau_M125_13TeV_powheg_pythia8/RunIISummer16NanoAOD_PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2-v1/"+name
    print "Adding file: ",aFile
    print "Making the MuTau tree"
    aROOTFile = TFile.Open(aFile)
    aTree = aROOTFile.Get("Events")
    print "TTree entries: ",aTree.GetEntries()
    HMuTauhTreeFromNano(aTree,doSvFit,applyRecoil).Loop()
    print "Making the TauTau tree"
    aROOTFile = TFile.Open(aFile)
    aTree = aROOTFile.Get("Events")
    HTauhTauhTreeFromNano(aTree,doSvFit,applyRecoil).Loop()

#Produce framework report required by CRAB
print "Generate framework report for CRAB"
#Empty list of input files to avoid CMSSW exception due to incorrect input
process.source.fileNames = []
#Produce new configuration file with an updated source
outFile = open("PSetTmp.py","w")
outFile.write(process.dumpPython())
outFile.close()
command = "cmsRun -j FrameworkJobReport.xml -p PSetTmp.py"
os.system(command)
os.system("rm PSetTmp.py")

exit(0)
