from CRABClient.UserUtilities import config, getUsernameFromSiteDB
#submit with 'python crab.py'
#Don't write to my directory (schoef), though

config = config()
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'ConfFile_cfg.py'
config.JobType.outputFiles = ['tuple.root']
config.Data.inputDBS = 'global'
config.Data.lumiMask =  'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'

config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 10
#config.Data.totalUnits = 300
config.Data.publication = False


#config.Data.ignoreLocality = True                                                                                                                                                                     
#config.Site.whitelist = ['T1_US_FNAL_Disk'] 

#config.Data.outLFNDirBase = '' 
#config.Data.publishDataName = ''

#config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())  #'/store/group/phys_jetmet/beranek/private4TSkim_24Jul2015_update/'
config.Site.storageSite = 'T2_BE_IIHE'
#config.Site.blacklist = ['T2_DE_DESY','T2_PL_Swierk']

datasets=[
'/SingleMuon/Run2016B-07Aug17_ver1-v1/AOD',
'/SingleMuon/Run2016B-07Aug17_ver2-v1/AOD',
'/SingleMuon/Run2016C-07Aug17-v1/AOD',
'/SingleMuon/Run2016D-07Aug17-v1/AOD',
'/SingleMuon/Run2016E-07Aug17-v1/AOD',
'/SingleMuon/Run2016F-07Aug17-v1/AOD',
'/SingleMuon/Run2016G-07Aug17-v1/AOD',
'/SingleMuon/Run2016H-07Aug17-v1/AOD'
]

if __name__ == '__main__':
    from CRABAPI.RawCommand import crabCommand
    for dataset in datasets:
        config.Data.inputDataset = dataset
        config.General.requestName = dataset.rstrip('/').lstrip('/').replace('/','_')
#        print config.General.requestName
        crabCommand('submit', config = config)
