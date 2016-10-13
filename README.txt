git clone https://github.com/lperrini/ZmumuPreselection.git

### main analyser is ZmumuAnalysis/MuMu_Analysis.cc
### which is connected to the ZmumuAnalysis/mutau_Tree.h
### where all the branches of the preselected tree
### are defined and assigned

### how to run the preselection

cd ZmumuAnalysis
source Make.sh MuMu_Analysis.cc

### create the script for the main submission

python genjobs.py
cd launchjobs
sh sbatch.sh

### all the output root files will be stored in
### /hdfs/local/lucia/PreselectedNtuples
### please change this path in genjobs.py before running it!

### now you have to merge all your files, please edit and then run
source MergeFiles.sh

### from this merging are excluded the TT/*.root outputs and the DYJets/*.root 
### these samples need to be merged asking for multiple files given their big size
### at this purpose you can use

perl mergeRootFilesTT.pl 10 /hdfs/local/lucia/PreselectedNtuples/TT/
perl mergeRootFilesDY.pl 50 /hdfs/local/lucia/PreselectedNtuples/DYJets/

### please note that the numbers '10' and '50' are important as they will give 15 files for TT and 25 files for DYJets and this setup is
### respected in the submit_*py cfg files used to run the final analysis on top of these files

### NB: these preselection is run on quasar, while the final selection and the final fits are run on lxplus
### for this reason, I report here the next steps to copy the preselected files to lxplus, and then to the eos
### area of Christian (files are heavy).

### copy files to local lxplus
rsync -avz --progress *.root lperrini@lxplus.cern.ch://afs/cern.ch/work/l/lperrini/ZtautauAnalysis_newJEC/CMSSW_7_6_3_patch2/src/UserCode/ztautau_fwk/data/ntuples/preselected/

### please check carefully these copies - they can fail compromizing the effectiveness of the final selection.

### once on lxplus, do 
cmsStage -f file.root /store/user/veelken/Lucia/...

### this eos area will become your input dir to run the final selection!!!
 

