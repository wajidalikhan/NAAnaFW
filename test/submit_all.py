#!/usr/bin/env python
#usage: python submit_all.py -c topplusdmTrees_cfg.py -f samples.txt -d crab_workdir
"""
This is a small script that submits a config over many datasets
"""
import os
from optparse import OptionParser

def getOptions() :
    """
    Parse and return the arguments provided by the user.
    """
    usage = ('usage: python submit_all.py -c CONFIG -d DIR ')

    parser = OptionParser(usage=usage)    
    parser.add_option("-c", "--config", dest="config",
        help=("The crab script you want to submit "),
        metavar="CONFIG")
    parser.add_option("-d", "--dir", dest="dir",
        help=("The crab directory you want to use "),
        metavar="DIR")
    parser.add_option("-f", "--datasets", dest="datasets",
        help=("File listing datasets to run over"),
        metavar="FILE")
    (options, args) = parser.parse_args()


    if options.config == None or options.dir == None:
        parser.error(usage)
    
    return options
    

def main():

    options = getOptions()

    from WMCore.Configuration import Configuration
    config = Configuration()

    from CRABAPI.RawCommand import crabCommand
    from httplib import HTTPException

    config.section_("General")
    config.General.workArea = options.dir
   
    config.section_("JobType")
    config.JobType.pluginName = 'Analysis'
    config.JobType.psetName = options.config
    config.JobType.allowUndistributedCMSSW = True
    config.JobType.pyCfgParams = ["isData=False", "changeJECs=True"]
    config.JobType.inputFiles = ["Spring16_25nsV6_MC_L2L3Residual_AK4PFchs.txt",#Old MC
                                 "Spring16_25nsV6_MC_L3Absolute_AK4PFchs.txt",
                                 "Spring16_25nsV6_MC_L1FastJet_AK4PFchs.txt",
                                 "Spring16_25nsV6_MC_L2Relative_AK4PFchs.txt",
                                 "Spring16_25nsV10_MC_L2L3Residual_AK4PFchs.txt",#New MC
                                 "Spring16_25nsV10_MC_L3Absolute_AK4PFchs.txt",
                                 "Spring16_25nsV10_MC_L1FastJet_AK4PFchs.txt",
                                 "Spring16_25nsV10_MC_L2Relative_AK4PFchs.txt",
                                 "Spring16_25nsV6_DATA_UncertaintySources_AK4PFchs.txt",#Old Data
                                 "Spring16_25nsV6_DATA_L1FastJet_AK4PFchs.txt",
                                 "Spring16_25nsV6_DATA_L1RC_AK4PFchs.txt",
                                 "Spring16_25nsV6_DATA_L2Relative_AK4PFchs.txt",
                                 "Spring16_25nsV6_DATA_L3Absolute_AK4PFchs.txt",
                                 "Spring16_25nsV6_DATA_L2L3Residual_AK4PFchs.txt",
                                 "Spring16_25nsV10BCD_DATA_UncertaintySources_AK4PFchs.txt",#RunBCD
                                 "Spring16_25nsV10BCD_DATA_L1FastJet_AK4PFchs.txt",
                                 "Spring16_25nsV10BCD_DATA_L1RC_AK4PFchs.txt",
                                 "Spring16_25nsV10BCD_DATA_L2Relative_AK4PFchs.txt",
                                 "Spring16_25nsV10BCD_DATA_L3Absolute_AK4PFchs.txt",
                                 "Spring16_25nsV10BCD_DATA_L2L3Residual_AK4PFchs.txt",
                                 "Spring16_25nsV10E_DATA_UncertaintySources_AK4PFchs.txt",#RunE
                                 "Spring16_25nsV10E_DATA_L1FastJet_AK4PFchs.txt",
                                 "Spring16_25nsV10E_DATA_L1RC_AK4PFchs.txt",
                                 "Spring16_25nsV10E_DATA_L2Relative_AK4PFchs.txt",
                                 "Spring16_25nsV10E_DATA_L3Absolute_AK4PFchs.txt",
                                 "Spring16_25nsV10E_DATA_L2L3Residual_AK4PFchs.txt",
                                 "Spring16_25nsV10F_DATA_UncertaintySources_AK4PFchs.txt",#RunF
                                 "Spring16_25nsV10F_DATA_L1FastJet_AK4PFchs.txt",
                                 "Spring16_25nsV10F_DATA_L1RC_AK4PFchs.txt",
                                 "Spring16_25nsV10F_DATA_L2Relative_AK4PFchs.txt",
                                 "Spring16_25nsV10F_DATA_L3Absolute_AK4PFchs.txt",
                                 "Spring16_25nsV10F_DATA_L2L3Residual_AK4PFchs.txt",
                                 "Spring16_25nsV10p2_DATA_UncertaintySources_AK4PFchs.txt",#RunGH
                                 "Spring16_25nsV10p2_DATA_L1FastJet_AK4PFchs.txt",
                                 "Spring16_25nsV10p2_DATA_L1RC_AK4PFchs.txt",
                                 "Spring16_25nsV10p2_DATA_L2Relative_AK4PFchs.txt",
                                 "Spring16_25nsV10p2_DATA_L3Absolute_AK4PFchs.txt",
                                 "Spring16_25nsV10p2_DATA_L2L3Residual_AK4PFchs.txt",
                                 ]


    config.section_("Data")
    config.Data.inputDataset = None
    config.Data.ignoreLocality = True
    #    config.Data.outLFNDirBase = '/store/user/oiorio/ttDM/samples/2016/Oct/'
    config.Data.outLFNDirBase = '/store/user/oiorio/ttDM/trees/2016/Dec/12Dec/'
    config.Data.inputDBS = 'phys03'
    config.Data.splitting = 'FileBased'
    #    config.Data.totalUnits = -1
    config.Data.unitsPerJob = 5
    config.Data.publication = False

    config.section_("Site")
    config.Site.storageSite = 'T2_CH_CSCS'
    #config.Site.whitelist = ['T2_CH_CERN','T2_IT_*','T2_DE_*','T2_CH_*']


    print 'Using config ' + options.config
    print 'Writing to directory ' + options.dir

    
    def submit(config):
        try:
            crabCommand('submit', config = config)
        except HTTPException, hte:
            print 'Cannot execute command'
            print hte.headers
            print hte
            print HTTPException

    #############################################################################################
    ## From now on that's what users should modify: this is the a-la-CRAB2 configuration part. ##
    #############################################################################################

    datasetsFile = open( options.datasets )
    jobsLines = datasetsFile.readlines()
    jobs = []
    for ijob in jobsLines :
        s = ijob.rstrip()
        jobs.append( s )
        print '  --> added ' + s
        
        #if s.find('ST_t-') > 0 :
        #   print ' single top\n   => add lhesource=source'
        config.JobType.pyCfgParams = ["isData=False", "changeJECs=True"]
        #        config.JobType.pyCfgParams = ["isData=False", "changeJECs=False","lhesource=source"]
        print s.split('/')[1]

    for ijob, job in enumerate(jobs) :

        ptbin = job.split('/')[1]
        cond = job.split('/')[2]#B2GAnaFW_80x_V2p0

        print "jobsplit11 ", job.split('/')[1]
#        print "cfgparams befor", config.JobType.pyCfgParams
        if ("QCD_Pt" in job.split('/')[1]) or ("ST_tW" in job.split('/')[1]):
            config.JobType.pyCfgParams= ["isData=False", "changeJECs=True","channel=qcd"]
        else:
            config.JobType.pyCfgParams = ["isData=False", "changeJECs=True"]
        #        ptbin
        config.General.requestName = '80xV2_' + ptbin #+ cond[10:] 
        config.Data.inputDataset = job
        config.Data.outputDatasetTag = '80xV2_' + ptbin #+ cond[10:]
        print 'Submitting ' + config.General.requestName + ', dataset = ' + job
        print 'Configuration :'
        print config
        try :
            from multiprocessing import Process
            p = Process(target=submit, args=(config,))
            p.start()
            p.join()
            #submit(config)
        except :
            print 'Not submitted.'
        





if __name__ == '__main__':
    main()            
