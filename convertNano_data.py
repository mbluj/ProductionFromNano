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
    "0409055F-7A0C-E811-B0C3-0025901D08F0.root",
]

#lumisToProcess = process.source.lumisToProcess
import FWCore.ParameterSet.Config as cms
import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.ParameterSet.Types as CfgTypes
lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange())
# get JSON file correctly parced
jsonFile = "https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt"
removeJson=False
if jsonFile.find("http")>-1:
    command = "wget "+jsonFile
    os.system(command)
    jsonFile=jsonFile.split("/")[-1]
    removeJson=True
myList = LumiList.LumiList(filename = jsonFile).getCMSSWString().split(',')
lumisToProcess.extend(myList)
if removeJson:
    command = "rm "+jsonFile
    os.system(command)
from ROOT import vector
vlumis = vector('string')()
for lumi in lumisToProcess:
    vlumis.push_back(lumi)

for name in fileNames:
    aFile = "file:///home/mbluj/work/data/NanoAOD/80X_with944/SingleMuon/Run2016E/05Feb2018-v1/"+name
    print "Adding file: ",aFile
    print "Making the MuTau tree"
    aROOTFile = TFile.Open(aFile)
    aTree = aROOTFile.Get("Events")
    print "TTree entries: ",aTree.GetEntries()
    HMuTauhTreeFromNano(aTree,doSvFit,applyRecoil,vlumis).Loop()
    print "Making the TauTau tree"
    aROOTFile = TFile.Open(aFile)
    aTree = aROOTFile.Get("Events")
    #HTauhTauhTreeFromNano(aTree,doSvFit,applyRecoil,vlumis).Loop()

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
