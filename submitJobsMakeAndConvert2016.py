#!/usr/bin/env python

import os, re
import commands
import math
import urllib

from crab3 import *
from mergeROOTFiles import *


submitJobs = True
mergeJobs = not submitJobs
#########################################
#########################################
def prepareCrabCfg(dataset,
                   inputDBS,
                   crabCfgName,
                   eventsPerJob,
                   jsonFile,
                   storage_element,
                   publish_data_suffix):

    workdir = publish_data_suffix
    shortName = dataset.split("/")[1]
    if dataset.split("/")[2].find("Run201")!=-1:
        shortName += "_"+dataset.split("/")[2]

    shortName = shortName.replace("-","_")
    shortName = shortName.split("_")[0]+shortName.split("_")[1]+shortName.split("_")[2]

    if dataset.find("PromptReco-v")!=-1:
        shortName+= "_v"+dataset[dataset.find("PromptReco-v")+12:dataset.find("PromptReco-v")+13]

    if dataset.find("23Sep2016-v")!=-1:
        shortName+= "_v"+dataset[dataset.find("23Sep2016-v")+11:dataset.find("23Sep2016-v")+12]

    if dataset.find("03Feb2017")!=-1:
        patternEnd = dataset.find("/MINIAOD")
        shortName+= dataset[dataset.find("03Feb2017")+9:patternEnd]

    if dataset.find("ext")!=-1:
        shortName+= "_"+dataset[dataset.find("ext"):dataset.find("ext")+4]

    if dataset.find("part")!=-1:
        shortName+= "_"+dataset[dataset.find("part"):dataset.find("part")+6]

    if dataset.find("t-channel")!=-1:
        shortName+= "_"+dataset[dataset.find("channel")+7:dataset.find("channel")+15]

    shortName = shortName.rstrip("-")
    shortName+="_"+publish_data_suffix

    ##Modify CRAB3 configuration
    config.JobType.psetName = 'PSet.py'
    isWZH = False
    if dataset.split("/")[2].find("JetsToLL")!=-1 or dataset.split("/")[2].find("JetsToLNu")!=-1 or dataset.split("/")[2].find("HToTauTau")!=-1:
        isWZH = True
    if isWZH:
        config.JobType.scriptExe = 'convertNano_recoil.py' #RECOIL
    else:
        config.JobType.scriptExe = 'convertNano_no-recoil.py' #No-RECOIL

    config.JobType.disableAutomaticOutputCollection = True

    config.JobType.outputFiles = ['HTTMT_HTauTauAnalysis.root', 'HTTTT_HTauTauAnalysis.root']
    config.JobType.inputFiles = [
        'NanoEventsSkeleton.C', 'NanoEventsSkeleton.h',
        'HTauTauTreeFromNanoBase.C', 'HTauTauTreeFromNanoBase.h',
        'HMuTauhTreeFromNano.C', 'HMuTauhTreeFromNano.h',
        'HTauhTauhTreeFromNano.C', 'HTauhTauhTreeFromNano.h',
        'HTTEvent.cxx', 'HTTEvent.h',
        'AnalysisEnums.h', 'PropertyEnum.h', 'TriggerEnum.h', 'SelectionBitsEnum.h', 'JecUncEnum.h',
        'zpt_weights_summer2016.root', 'zpt_weights_2016_BtoH.root',
        'Summer16_23Sep2016V4_MC_UncertaintySources_AK4PFchs.txt'
    ]

    config.Site.storageSite = storage_element
    config.General.requestName = shortName

    config.Data.inputDataset = dataset
    config.Data.outLFNDirBase = '/store/user/mbluj/WAWNTupleFromNanoAOD/'+publish_data_suffix+"/"
    config.Data.outputDatasetTag = shortName
    config.Data.inputDBS = inputDBS

    #DYJets
    if dataset.split("/")[2].find("Jets")!=-1:
        eventsPerJob = 40000
    #DY and W 3,4 Jets
    if dataset.split("/")[2].find("3Jets")!=-1 or dataset.split("/")[2].find("4Jets")!=-1:
        eventsPerJob = 1000

    #config.Data.splitting = 'EventAwareLumiBased'
    #config.Data.unitsPerJob = eventsPerJob
    config.Data.splitting = 'FileBased'
    config.Data.unitsPerJob = 1 #number of files per jobs
    config.Data.totalUnits = -1

    config.Data.lumiMask=""
    if dataset.split("/")[2].find("Run201")!=-1:
        command = "wget "+jsonFile
        os.system(command)
        config.Data.lumiMask=jsonFile.split("/")[-1]
        config.JobType.psetName = 'analyzerData.py'
    out = open('crabTmp.py','w')
    out.write(config.pythonise_())
    out.close()
    #FIXMEos.system("crab submit -c crabTmp.py")
    os.system("rm -f "+jsonFile.split("/")[-1])
#########################################
#########################################
eventsPerJob = 100000 #Wjets and DYJets hardoced in code above
#eventsPerJob = 200000#4Mu analysis

#FIXME: datasets refer to MiniAOD for now, kept as a placeholder!!!
from datasetsMoriond17 import datasets


##TEST
datasets = [
    "/VBFHToTauTau_M125_13TeV_powheg_pythia8/bluj-RunIISummer16NanoAOD940-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v3-8c9412e29164719d2de5b19c3daa1df7/USER",
]

###############
jsonFileReReco = "https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt"
########################################################
if submitJobs:
    for dataset in datasets:
        jsonFile2016 = jsonFileReReco

        prepareCrabCfg(crabCfgName="crab3.py",
                       dataset=dataset,
                       #inputDBS = 'global',
                       inputDBS = 'phys03',
                       eventsPerJob=eventsPerJob,
                       jsonFile=jsonFile2016,
                       storage_element="T2_PL_Swierk",
                       publish_data_suffix = "test_v1")                  
########################################################
########################################################
## Merge output ROOT files.
########################################################
if mergeJobs:
    for dataset in datasets:
        mergeDataset(dataset=dataset, publish_data_suffix = "test_v1",
                                      outputDir="/mnt/home/mbluj/work/data/WAWNTuples/FromNano/2016/NTUPLES_23_01_2018/")

#for a in v1/*v7_SM*; do crab resubmit -d $a; doneQ
#for a in v1/*Run2016*v7_SM*; do crab report -d $a; done

#mergeJSON.py crab_SingleMuonRun2016*23*/results/processedLumis.json crab_SingleMuonRun2016H*/results/processedLumis.json > processedLumis_SingleMuon.json
#mergeJSON.py crab_TauRun2016*23*/results/processedLumis.json crab_TauRun2016H*/results/processedLumis.json > processedLumis_Tau.json
#for a in *json; do echo $a >>  lumi.out; ~/.local/bin/brilcalc lumi --normtag /afs/cern.ch/user/l/lumipro/public/normtag_file/normtag_DATACERT.json -i $a; done >>  lumi.out
#grep -A 5 'Summary\|Run2016' lumi.out

'''

'''
