from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'Dummy'
config.General.workArea = 'v1'
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'PSet.py'
config.JobType.sendExternalFolder = True

config.section_("Data")
config.Data.inputDataset = 'Dummy'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1 #number of files per jobs
config.Data.totalUnits =  -1 #number of event
config.Data.outLFNDirBase = '/store/user/mbluj/WAWNTupleFromNanoAOD/'
config.Data.publication = False
config.Data.outputDatasetTag = 'Dummy'



config.section_("Site")
config.Site.storageSite = 'T2_PL_Swierk'
config.Site.blacklist = ['T2_KR_*','T2_CN_*','T2_BR_*','T2_US_Florida','T2_US_UCSD']
