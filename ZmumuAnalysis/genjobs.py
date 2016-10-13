import codecs, jinja2, os, stat, ntpath, sys, logging
from collections import OrderedDict as OD

########################## DEFINE SAMPLES ##############################

samples = OD()

### SM Higgs 125
samples["VBFHToTauTau_M125"] = OD([
  ("path", ["/hdfs/local/lucia/FinalNtuples/VBFHToTauTau_M125"]),
  ("nof_files", [151]), # expected
  ("is_mc", True),
  ("use_it", True)
])
samples["GluGluHToTauTau_M125"] = OD([
  ("path", ["/hdfs/local/lucia/FinalNtuples/GluGluHToTauTau_M125"]),
  ("nof_files", [151]), # expected
  ("is_mc", True),
  ("use_it", True)
])

### Diboson
samples["ZZTo2L2Q"] = OD([
  ("path", ["/hdfs/local/lucia/FinalNtuples/ZZTo2L2Q"]),
  ("nof_files", [3104]), # expected
  ("is_mc", True),
  ("use_it", True)
])
samples["ZZTo4L"] = OD([
  ("path", ["/hdfs/local/andres/Ntuples_Lucia/ZZTo4L"]),
  ("nof_files", [2141]), # expected
  ("is_mc", True),
  ("use_it", True)
])
samples["WW"] = OD([
  ("path", ["/hdfs/local/andres/Ntuples_Lucia/WW"]),
  ("nof_files", [1036]), # expected
  ("is_mc", True),
  ("use_it", True)
])
samples["WZTo2L2Q"] = OD([
  ("path", ["/hdfs/local/lucia/FinalNtuples/WZTo2L2Q"]),
  ("nof_files", [5144]), # expected 5228
  ("is_mc", True),
  ("use_it", True)
])
samples["WZJets"] = OD([
  ("path", ["/hdfs/local/lucia/FinalNtuples/WZJets"]),
  ("nof_files", [2536]), # expected
  ("is_mc", True),
  ("use_it", True)
])
samples["WZTo1L3Nu"] = OD([
  ("path", ["/hdfs/local/lucia/FinalNtuples/WZTo1L3Nu"]),
  ("nof_files", [341]), # expected 343
  ("is_mc", True),
  ("use_it", True)
])
samples["WZTo1L1Nu2Q"] = OD([
  ("path", ["/hdfs/local/andres/Ntuples_Lucia/WZTo1L1Nu2Q"]),
  ("nof_files", [4010]), # expected
  ("is_mc", True),
  ("use_it", True)
])
### Single-top
samples["ST_t_antitop"] = OD([
  ("path", ["/hdfs/local/lucia/FinalNtuples/ST_t_antitop"]),
  ("nof_files", [325]), # expected
  ("is_mc", True),
  ("use_it", True)
])
samples["ST_t_top"] = OD([
  ("path", ["/hdfs/local/lucia/FinalNtuples/ST_t_top"]),
  ("nof_files", [660]), # expected
  ("is_mc", True),
  ("use_it", True)
])
samples["ST_tW_antitop"] = OD([
  ("path", ["/hdfs/local/lucia/FinalNtuples/ST_tW_antitop"]),
  ("nof_files", [200]), # expected
  ("is_mc", True),
  ("use_it", True)
])
samples["ST_tW_top"] = OD([
  ("path", ["/hdfs/local/lucia/FinalNtuples/ST_tW_top"]),
  ("nof_files", [200]), # expected
  ("is_mc", True),
  ("use_it", True)
])

### TT
samples["TT"] = OD([
  ("path", ["/hdfs/local/lucia/FinalNtuples/TT_ext3",
            "/hdfs/local/lucia/FinalNtuples/TT_ext4"]),
  ("nof_files", [4864,
                 9269]),
  ("is_mc", True),
  ("use_it", True)
])
 
### W+Jets
samples["WJets"] = OD([
  ("path", ["/hdfs/local/karl/FinalNtuples/WJets"]),
  ("nof_files", [9515]), # expected
  ("is_mc", True),
  ("use_it", True)
])

samples["W1Jets"] = OD([
  ("path", ["/hdfs/local/andres/Ntuples_Lucia/W1Jets"]),
  ("nof_files", [9304]), # expected 9316
  ("is_mc", True),
  ("use_it", True)
])

samples["W2Jets"] = OD([
  ("path", ["/hdfs/local/karl/FinalNtuples/W2Jets"]),
  ("nof_files", [6156]), # expected
  ("is_mc", True),
  ("use_it", True)
])

samples["W3Jets"] = OD([
  ("path", ["/hdfs/local/karl/FinalNtuples/W3Jets"]),
  ("nof_files", [3786]), # expected
  ("is_mc", True),
  ("use_it", True)
])

samples["W4Jets"] = OD([
  ("path", ["/hdfs/local/karl/FinalNtuples/W4Jets"]),
  ("nof_files", [1793]), # expected
  ("is_mc", True),
  ("use_it", True)
])

### DY+jets
samples["DYJets"] = OD([
  ("path", ["/hdfs/local/lucia/FinalNtuples/DYJets_ext",
            "/hdfs/local/lucia/FinalNtuples/DYJets_ext_2",
            "/hdfs/local/lucia/FinalNtuples/DYJets_ext_3",
            "/hdfs/local/lucia/FinalNtuples/DYJets_ext_4",
            "/hdfs/local/lucia/FinalNtuples/DYJets_ext_5",
            "/hdfs/local/lucia/FinalNtuples/DYJets_ext_6",
            "/hdfs/local/lucia/FinalNtuples/DYJets_ext_7",
            "/hdfs/local/lucia/FinalNtuples/DYJets_ext_8",
            "/hdfs/local/lucia/FinalNtuples/DYJets_ext_9",
            "/hdfs/local/lucia/FinalNtuples/DYJets_ext_10",
            "/hdfs/local/lucia/FinalNtuples/DYJets_ext_11",
            "/hdfs/local/lucia/FinalNtuples/DYJets_ext_12",
            "/hdfs/local/lucia/FinalNtuples/DYJets_ext_13",
            "/hdfs/local/lucia/FinalNtuples/DYJets_ext_14",
            "/hdfs/local/lucia/FinalNtuples/DYJets_ext_15",
            "/hdfs/local/lucia/FinalNtuples/DYJets_ext_16",
            "/hdfs/local/lucia/FinalNtuples/DYJets_ext_17",
            "/hdfs/local/lucia/FinalNtuples/DYJets_ext_18",
            "/hdfs/local/lucia/FinalNtuples/DYJets_ext_19",
            "/hdfs/local/lucia/FinalNtuples/DYJets_ext_20"]),
  ("nof_files", [8363, #expected 8370
                 7320,
                 7199,
                 6005,
                 6179, #expected 6181
                 6774,
                 5575,
                 6962,
                 6178,
                 5742,
                 6811,
                 4833,
                 5792,
                 6122,
                 5629,
                 6135,
                 6065,
                 5927,
                 6117,
                 3819]),
  ("is_mc", True),
  ("use_it", True)
])

samples["DY1Jets"] = OD([
  ("path", ["/hdfs/local/lucia/FinalNtuples/DY1Jets"]),
  ("nof_files", [6583]), 
  ("is_mc", True),
  ("use_it", True)
])
samples["DY2Jets"] = OD([
  ("path", ["/hdfs/local/lucia/FinalNtuples/DY2Jets"]),
  ("nof_files", [2020]), 
  ("is_mc", True),
  ("use_it", True)
])
samples["DY3Jets"] = OD([
  ("path", ["/hdfs/local/lucia/FinalNtuples/DY3Jets"]),
  ("nof_files", [576]), 
  ("is_mc", True),
  ("use_it", True)
])
samples["DY4Jets"] = OD([
  ("path", ["/hdfs/local/lucia/FinalNtuples/DY4Jets"]),
  ("nof_files", [423]), 
  ("is_mc", True),
  ("use_it", True)
])
samples["DYJetsLowMass"] = OD([
  ("path", ["/hdfs/local/lucia/FinalNtuples/DYJets_LowMass"]),
  ("nof_files", [8990]), 
  ("is_mc", True),
  ("use_it", True)
])

### EWK W/Z+ 2 jets
samples["EWKWMinus"] = OD([
  ("path", ["/hdfs/local/lucia/FinalNtuples/EWKWMinus"]),
  ("nof_files", [99]), 
  ("is_mc", True),
  ("use_it", True)
])
samples["EWKWPlus"] = OD([
  ("path", ["/hdfs/local/lucia/FinalNtuples/EWKWPlus"]),
  ("nof_files", [100]), 
  ("is_mc", True),
  ("use_it", True)
])
samples["EWKZToLL"] = OD([
  ("path", ["/hdfs/local/lucia/FinalNtuples/EWKZToLL"]),
  ("nof_files", [31]), 
  ("is_mc", True),
  ("use_it", True)
])
samples["EWKZToNuNu"] = OD([
  ("path", ["/hdfs/local/lucia/FinalNtuples/EWKZToNuNu"]),
  ("nof_files", [60]), 
  ("is_mc", True),
  ("use_it", True)
])

### SingleMuon data
samples["SingleMuon_Run2015C"] = OD([
  ("path", ["/hdfs/local/lucia/FinalNtuples/SingleMuon_Run2015C"]),
  ("nof_files", [157]), 
  ("is_mc", False),
  ("use_it", True)
])
samples["SingleMuon_Run2015D"] = OD([
  ("path", ["/hdfs/local/lucia/FinalNtuples/SingleMuon_Run2015D"]),
  ("nof_files", [5113]), #expected 5120 
  ("is_mc", False),
  ("use_it", True)
])




########################################################################

cmd_template = """#!/bin/sh

hostname
date

# set up the environment
pwd
export SCRAM_ARCH=slc6_amd64_gcc493
export BUILD_ARCH=slc6_amd64_gcc493
export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
source /home/lucia/setup_cms.sh
cd /home/lucia/ZtautauXsection/CMSSW_7_6_3_patch2/src
pwd
cmsenv

# copy the *.root files to /scratch
mkdir -p {{ scratch_dir_output }}
mkdir -p {{ hdfs_outputpath_sample }}
{% for hdfs_path,scratch_path in cmd_meta %}
cp {{ hdfs_path }} {{ scratch_path }} {% endfor %}

date

# create the input file list
ls -l {{ scratch_dir }}/*.root | awk '{print $9}' > {{ scratch_dir }}/file_list.txt
# run the analysis
/home/lucia/ZtautauXsection/CMSSW_7_6_3_patch2/src/ZmumuAnalysis/./MuMu_Analysis.exe {{ sample_type }} {{ output_filename }} {{ scratch_dir }}/file_list.txt
date

# post-processing

# move the output file to hdfs
mv {{ output_filename }} {{ hdfs_output_filename }}

# remove the input files from scratch
rm {{ scratch_dir }}/*.root
date
echo "Done"

"""

sbatch_template = """#!/bin/bash
{% for logfile,command in sbatch_meta %}
sbatch --partition={{ sbatch_queue }} --output={{ logfile }} {{ command }}{% endfor %}

"""

def create_cmd(scratch_dir, hdfs_outputpath_sample, hdfs_paths, scratch_paths, sample_type, output_filename, hdfs_output_filename):
  '''
  Creates a bash script from the template defined above
  '''
  cmd_meta = zip(hdfs_paths, scratch_paths)
  return jinja2.Template(cmd_template).render(
    scratch_dir = scratch_dir,
    scratch_dir_output = scratch_dir+"/output",
    cmd_meta = cmd_meta,
    sample_type = sample_type,
    output_filename = output_filename,
    hdfs_outputpath_sample = hdfs_outputpath_sample,
    hdfs_output_filename = hdfs_output_filename)

def create_sbatch(sbatch_logfiles, commands, sbatch_queue = "short"):
  '''
  Create the sbatch script from the template defined above
  '''
  sbatch_meta = zip(sbatch_logfiles, commands)
  return jinja2.Template(sbatch_template).render(sbatch_meta = sbatch_meta, sbatch_queue = sbatch_queue)

def generate_file_ids(nof_files, max_files_per_job):
  '''
  Creates a list of lists, where each element contains file IDs per chunk
  '''
  file_limits = range(1, nof_files, max_files_per_job)
  file_limits.append(nof_files + 1)
  job_ids = [range(file_limits[i], file_limits[i + 1]) for i in range(len(file_limits) - 1)]
  return job_ids

def add_chmodX(fullpath):
  '''
  chmod +x fullpath, basically
  '''
  st = os.stat(fullpath)
  os.chmod(fullpath, st.st_mode | stat.S_IEXEC)

def create_if_not_exists(dir_fullpath):
  if not os.path.exists(dir_fullpath): os.makedirs(dir_fullpath)

if __name__ == '__main__':

  # directory where you keep your *.sh and *.out (sbatch stderr,stdout) files,
  # and sbatch.sh which is the only file you are expected to run by hand
  #NOTE change the following path to something more appropriate
  parent_dir = '/home/lucia/ZtautauXsection/CMSSW_7_6_3_patch2/src/ZmumuAnalysis/launchJobs'
  cfg_storage = os.path.join(parent_dir, 'cfgs')
  log_storage = os.path.join(parent_dir, 'logs')
  create_if_not_exists(cfg_storage)
  create_if_not_exists(log_storage)
  sbatch_location = os.path.join(parent_dir, "sbatch.sh")
  
  # set up a logger (stdout & log.txt in parent_dir)
  # needed when trying to detect jobs w/ missing files
  # e.g. by grepping WARNING from the log file
  log = logging.getLogger('')
  log.setLevel(logging.DEBUG)
  fmt = logging.Formatter('%(asctime)s - %(levelname)s: %(message)s')

  ch = logging.StreamHandler(sys.stdout)
  ch.setFormatter(fmt)
  log.addHandler(ch)

  # remove old log file
  log_file = os.path.join(parent_dir, "log.txt")
  if os.path.exists(log_file) and os.path.isfile(log_file): os.remove(log_file)
  fh = logging.FileHandler(log_file)
  fh.setFormatter(fmt)
  log.addHandler(fh)

  # scratch storage -- temporary working directory
  scratch_storage = '/scratch/lucia'

  # final destination where you keep your *.root files
  hdfs_storage = '/hdfs/local/lucia/PreselectedNtuples'

  # how many files per job? a chunk
  split_by = 100
  # if you set it to True, the script checks which files are
  debug = False
  # input/output root file prefix
  file_prefix = "HTauTauAnalysis"

  # list of bash scripts and sbatch log files (one per chunk)
  job_list, log_list = [], []

  for sample_name, sample_info in samples.items():
    if not sample_info["use_it"]: continue
    if sample_info["is_mc"]:
       sample_type = "mc" 
    else:
       if sample_name == "SingleMuon_Run2015C":
          sample_type = "dataC" 
       if sample_name == "SingleMuon_Run2015D":
          sample_type = "dataD" 
    logging.info("Generating config & bash files for sample %s" % sample_name)

    # create a subdirectory for bash scripts and log files
    cfg_storage_sample = os.path.join(cfg_storage, sample_name)
    log_storage_sample = os.path.join(log_storage, sample_name)
    create_if_not_exists(cfg_storage_sample)
    create_if_not_exists(log_storage_sample)
    bash_filepath = os.path.join(cfg_storage_sample, "cfg") + "_%d.sh"

    # loop over multiple dataset per sample and create chunks for each job
    nof_paths = len(sample_info["nof_files"])
    job_ids = []
    for i in range(nof_paths):
      job_ids = job_ids + generate_file_ids(sample_info["nof_files"][i], split_by) 
    ntuple_paths = [os.path.join(x, file_prefix) + "_%d.root" for x in sample_info["path"]]

    # for debugging purposes
    total_nof_missing, total_nof_expected = 0, 0
    ntuple_paths_idx = -1
    for i in range(len(job_ids)):
      # select a chunk
      job_id = job_ids[i]
      # if we are about to assemble file list from a new dataset (DY!) switch there
      if job_id[0] == 1: ntuple_paths_idx += 1

      # set subdirectories and output filename in /scratch; in /hdfs for a chunk
      scratch_storage_sample = os.path.join(scratch_storage, sample_name, "chunk_%d" % i)
      #scratch_filepath_sample = os.path.join("output", file_prefix + "_%d.root" % i)
      scratch_filepath_sample = os.path.join(scratch_storage_sample, "output", file_prefix + "_%d.root" % i)
      hdfs_outputpath_sample = os.path.join(hdfs_storage, sample_name)
      hdfs_filepath_sample   = os.path.join(hdfs_storage, sample_name, file_prefix + "_%d.root" %i)

      # create the sample list to copy
      sample_filelist = [ntuple_paths[ntuple_paths_idx] % x for x in job_id]
      good_samples = sample_filelist
      if debug:
        # check if all the sample files are actually there
        missing_files = []
        good_samples = []
        for sample_file in sample_filelist:
          if not os.path.exists(sample_file):
            missing_files.append(ntpath.basename(sample_file)[len(file_prefix) + 1:-len(".root")])
          else:
            good_samples.append(sample_file)

        nof_missing = len(missing_files)
        total_nof_missing += nof_missing
        total_nof_expected += len(job_id)
        if nof_missing > 0:
          logging.warning("Encountered missing files in chunk #%d of sample %s: %s" % (i, sample_name, str(missing_files)))

      # do not write the *.sh file if none of the files are present
      if len(good_samples) == 0:
        logging.error("Skipping chunk #%d of sample %s because none of the files are present" % (i, sample_name))
        continue

      # create conjugate list from original files residing in /hdfs
      scratch_filelist = [os.path.join(scratch_storage_sample, ntpath.basename(x)) for x in good_samples]

      # log and bash files (needed by sbatch)
      log_filepath_sample = os.path.join(log_storage_sample, "log_%d") % i + "-%a.out"
      bash_filepath_sample = bash_filepath % i
      bash_contents = create_cmd(scratch_storage_sample, hdfs_outputpath_sample, good_samples, scratch_filelist,sample_type, scratch_filepath_sample, hdfs_filepath_sample)
      
      with codecs.open(bash_filepath_sample, 'w', 'utf-8') as f:
        f.write(bash_contents)

      job_list.append(bash_filepath_sample)
      log_list.append(log_filepath_sample)
      add_chmodX(bash_filepath_sample)

    if total_nof_expected > 0:
      logging.debug("Total number of files missing (expected): %d (%d) %.2f" % \
        (total_nof_missing, total_nof_expected, float(total_nof_missing) / total_nof_expected))

  logging.info("Creating sbatch file")
  sbatch_contents = create_sbatch(log_list, job_list)
  with codecs.open(sbatch_location, 'w', 'utf-8') as f:
    f.write(sbatch_contents)
  add_chmodX(sbatch_location)
  logging.info("Done")
