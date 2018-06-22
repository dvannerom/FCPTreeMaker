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
#config.Data.lumiMask =  'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions16/13TeV/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON.txt'

config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = 30000
config.Data.publication = False
#config.Data.totalUnits = 20000000

#config.Data.ignoreLocality = True                                                                                                                                                                     
#config.Site.whitelist = ['T1_US_FNAL_Disk'] 

#config.Data.outLFNDirBase = '' 
#config.Data.publishDataName = ''

#config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())  #'/store/group/phys_jetmet/beranek/private4TSkim_24Jul2015_update/'
config.Site.storageSite = 'T2_BE_IIHE'

datasets=[
'/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16DR80Premix-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/AODSIM',
#'/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16DR80Premix-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/AODSIM',
#'/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16DR80Premix-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/AODSIM'
]

if __name__ == '__main__':
    from CRABAPI.RawCommand import crabCommand
    for dataset in datasets:
        config.Data.inputDataset = dataset
        config.General.requestName = 'DYJetsToLL_M-10To50-amcatnloFXFX-pythia8_RunIISummer16-PUMoriond17_v6-v1_AODSIM_v2'
#        config.General.requestName = 'DYJetsToLL_M-50-amcatnloFXFX-pythia8_RunIISummer16-PUMoriond17_v6_ext2-v1_AOD_v2'
#        config.General.requestName = 'WJetsToLNu_amcatnloFXFX-pythia8_RunIISummer16-PUMoriond17_v6-v1_AODSIM_v2'
#        print config.General.requestName
        crabCommand('submit', config = config)
