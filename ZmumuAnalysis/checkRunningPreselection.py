###################################################################
#author L.Perrini
###################################################################

import os
import sys
import glob
import re, datetime
import codecs, jinja2, stat, ntpath, sys, logging

def modification_date(filename):
    t = os.path.getmtime(filename)
    return str(datetime.datetime.fromtimestamp(t))
###################################################################
#### Parameters to be changed for each production

files =["DYJetsLowMass","DY1Jets","DY2Jets","DY3Jets","DY4Jets","DYJets",
 "EWKWMinus","EWKWPlus","EWKZToLL","EWKZToNuNu",
 "GluGluHToTauTau_M125","ST_tW_antitop","ST_tW_top","ST_t_antitop","ST_t_top",
 "SingleMuon_Run2015C","SingleMuon_Run2015D","TT","VBFHToTauTau_M125",  
 "W1Jets", "W2Jets", "W3Jets", "W4Jets", "WJets",
 "WW", "WZJets", "WZTo1L1Nu2Q", "WZTo1L3Nu", "WZTo2L2Q", "ZZTo2L2Q", "ZZTo4L"
]

###################################################################
#### Automated script starting

# dataset block definition
#NumberFilesBegin = "==="
StorageDir = "/hdfs/local/lucia/PreselectedNtuples/"
LogsDir = "launchJobs/cfgs/"

nfiles = len(files)

#for ff in files:
#   for item in os.listdir(StorageDir+ff+"/"):
#      statinfo = os.stat(StorageDir+ff+"/"+item)
#      fileSize = statinfo.st_size
#      fileNumber = item.split("_")[1].split(".root")
#      d = modification_date(StorageDir+ff+"/"+item)
#      if(d.find("2016-09-01") < 0):
#         print "sbatch --partition=short --output=/home/lucia/ZtautauXsection/CMSSW_7_6_3_patch2/src/ZmumuAnalysis/launchJobs/logs/"+ff+"/log_"+fileNumber[0]+"-%a.out /home/lucia/ZtautauXsection/CMSSW_7_6_3_patch2/src/ZmumuAnalysis/launchJobs/cfgs/"+ff+"/cfg_"+fileNumber[0]+".sh"

for ff in files:
   filesPreselected = int(len(glob.glob(StorageDir+ff+"/*.root") ))
   logFileToPreselect = int(len(glob.glob(LogsDir+ff+"/*.sh") ))
   if (filesPreselected<logFileToPreselect):
      for item in os.listdir(LogsDir+ff+"/"):
         fileEnd    = item.split("_")
         fileNumber = item.split("_")[1].split(".sh")
         if not glob.glob(StorageDir+ff+"/HTauTauAnalysis_"+fileNumber[0]+".root"):
            print "sbatch --partition=short --output=/home/lucia/ZtautauXsection/CMSSW_7_6_3_patch2/src/ZmumuAnalysis/launchJobs/logs/"+ff+"/log_"+fileNumber[0]+"-%a.out /home/lucia/ZtautauXsection/CMSSW_7_6_3_patch2/src/ZmumuAnalysis/launchJobs/cfgs/"+ff+"/cfg_"+fileNumber[0]+".sh"
 
