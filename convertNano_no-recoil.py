#!/usr/bin/env python

import os

from ROOT import gSystem, TChain, TSystem, TFile
from PSet import process

#doSvFit = True
doSvFit = False
#applyRecoil=True
applyRecoil=False
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

from ROOT import vector
vlumis = vector('string')()
for lumi in process.source.lumisToProcess:
    vlumis.push_back(lumi)

for aFile in process.source.fileNames:
    aFile = aFile.replace("/store","root://cms-xrd-global.cern.ch///store")
    print "Adding file: ",aFile
    print "Making the MuTau tree"
    aROOTFile = TFile.Open(aFile)
    aTree = aROOTFile.Get("Events")
    print "TTree entries: ",aTree.GetEntries()
    HMuTauhTreeFromNano(aTree,doSvFit,applyRecoil,vlumis).Loop()
    print "Making the TauTau tree"
    aROOTFile = TFile.Open(aFile)
    aTree = aROOTFile.Get("Events")
    HTauhTauhTreeFromNano(aTree,doSvFit,applyRecoil,vlumis).Loop()

#Merge files.
command = "hadd -f HTTMT_HTauTauAnalysis.root HTTMT_*.root"
os.system(command)
command = "hadd -f HTTTT_HTauTauAnalysis.root HTTTT_*.root"
os.system(command)

print "Done!", "Processed ",len(process.source.fileNames), "files"

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
