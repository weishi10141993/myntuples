# Produce Ntuple from DUNE FD MC files

The following instruction is used to produce ROOT n-tuples from FD MC files (mcc11): [FD Beamsim Requests](https://dune-data.fnal.gov/mc/mcc11/index.html).

## DUNE FNAL machines (dunegpvm*) environment setup

[First time only]

```
cd /dune/app/users/weishi                                               # Replace with your username for all commands below
mkdir FDEff (first time only)
cd FDEff

source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
unsetup mrb
setup mrb v4_04_06 
setup dunetpc v09_22_02 -q e19:debug
[optional if run interactively] setup_fnal_security                     # A FNAL grid proxy to submit jobs and access data in dCache via xrootd or ifdh.

mrb newDev
# The prompt ask you to run this:
source /dune/app/users/<your_username>/inspect/localProducts_larsoft_${LARSOFT_VERSION}_debug_${COMPILER}/setup
# For example, mine is: source /dune/app/users/weishi/FDEff/localProducts_larsoft_v09_22_02_debug_e19/setup

cd srcs
git clone https://github.com/weishi10141993/myntuples.git               # First time only, checkout the analysis code from GitHub

mrb uc                                                                  # Tell mrb to update CMakeLists.txt with the latest version numbers of the products.
cd ${MRB_BUILDDIR}                                                      # Go to your build directory
mrb z
mrbsetenv                                                               # Create the bookkeeping files needed to compile programs.
mrb b                                                                   # Compile the code in ${MRB_SOURCE}
```

To run on FD MC files, this produces a ROOT nTuple:

```
cd /dune/app/users/weishi/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis
lar -c MyEnergyAnalysis.fcl -n -1
# Or run with nohup
nohup lar -c MyEnergyAnalysis.fcl -n -1 >& out_myntuple_nohup.log &             # check status: jobs -l
# 10k evts take about 32 minutes
```

The next time you login a DUNE FNAL machine (dunegpvm*), do the following to set up:

```
source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
unsetup mrb
setup mrb v4_04_06 
setup dunetpc v09_22_02 -q e19:debug
source /dune/app/users/weishi/FDEff/localProducts_larsoft_v09_22_02_debug_e19/setup
mrbsetenv
cd /dune/app/users/weishi/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis
```

If changed ```MyEnergyAnalysis_module.cc```, recompile the code (do the setup above first):

```
cd ${MRB_BUILDDIR}                        # Go to your build directory
mrb z                                     # Remove old build directory
mrbsetenv                                 # Create the bookkeeping files needed to compile programs.
# On dunegpvm
mrb b
# On ivy
mrb install                               # Compile the code in ${MRB_SOURCE} and put the results in ${MRB_INSTALL}
```

If added new package in ```srcs``` directory, do ```mrb uc``` and then recompile as above.

To commit changed code changes to remote repository:

```
git commit
git push
```

## Run grid jobs

The instruction is based on the [DUNE computing tutorial](https://wiki.dunescience.org/wiki/DUNE_Computing/Submitting_grid_jobs_May2021#Submit_a_job).

Once the above is compiled and runs without problem interactively, you can start to produce a tarball. First, you need to have a grid setup for the localProducts as the grid job typically runs on a different machine than your working machine,

```
cd /dune/app/users/weishi/FDEff/localProducts_larsoft_v09_22_02_debug_e19
cp setup setup-grid         # make a copy of the setup for grid job
```

then in ```setup-grid```, change ```/dune/app/users/weishi``` to the worker node working directory ```${_CONDOR_JOB_IWD}```.

Now get txt file that lists of input files and work env set up script:

```
cd /dune/app/users/weishi
wget https://raw.githubusercontent.com/weishi10141993/NeutrinoPhysics/main/MCC11FDBeamsim_nu_reco.txt --no-check-certificate
wget https://raw.githubusercontent.com/weishi10141993/NeutrinoPhysics/main/setupFDEffTarBall-grid.sh --no-check-certificate
```

Then make the tarball,
```
tar -czvf FDEff.tar.gz FDEff setupFDEffTarBall-grid.sh MCC11FDBeamsim_nu_reco.txt

# Check the tarball *.tar.gz is indeed created and open with: tar -xf *.tar.gz
```

Now get one of the following grid running scripts
```
wget https://raw.githubusercontent.com/weishi10141993/NeutrinoPhysics/main/run_FDEffTarBall_autogrid.sh --no-check-certificate
# Or this one that allows you to set the number of input files using line number in txt file:
wget https://raw.githubusercontent.com/weishi10141993/NeutrinoPhysics/main/run_FDEffTarBall_grid.sh --no-check-certificate

```

Finally you can submit the job:

this submits N jobs (since we have 9914 files, N=9914 will run 1 files/job),
```
jobsub_submit -G dune -N 9914 --memory=1000MB --disk=1GB --expected-lifetime=30m --cpu=1 --resource-provides=usage_model=DEDICATED,OPPORTUNISTIC,OFFSITE --tar_file_name=dropbox:///dune/app/users/weishi/FDEff.tar.gz --use-cvmfs-dropbox -l '+SingularityImage=\"/cvmfs/singularity.opensciencegrid.org/fermilab/fnal-wn-sl7:latest\"' --append_condor_requirements='(TARGET.HAS_Singularity==true&&TARGET.HAS_CVMFS_dune_opensciencegrid_org==true&&TARGET.HAS_CVMFS_larsoft_opensciencegrid_org==true&&TARGET.CVMFS_dune_opensciencegrid_org_REVISION>=1105&&TARGET.HAS_CVMFS_fifeuser1_opensciencegrid_org==true&&TARGET.HAS_CVMFS_fifeuser2_opensciencegrid_org==true&&TARGET.HAS_CVMFS_fifeuser3_opensciencegrid_org==true&&TARGET.HAS_CVMFS_fifeuser4_opensciencegrid_org==true)' file:///dune/app/users/weishi/run_FDEffTarBall_autogrid.sh
```

This submit 1 job that will run sequentially your specified range of files,

```
jobsub_submit -G dune -N 1 --memory=1000MB --disk=1GB --expected-lifetime=30m --cpu=1 --resource-provides=usage_model=DEDICATED,OPPORTUNISTIC,OFFSITE --tar_file_name=dropbox:///dune/app/users/weishi/FDEff.tar.gz --use-cvmfs-dropbox -l '+SingularityImage=\"/cvmfs/singularity.opensciencegrid.org/fermilab/fnal-wn-sl7:latest\"' --append_condor_requirements='(TARGET.HAS_Singularity==true&&TARGET.HAS_CVMFS_dune_opensciencegrid_org==true&&TARGET.HAS_CVMFS_larsoft_opensciencegrid_org==true&&TARGET.CVMFS_dune_opensciencegrid_org_REVISION>=1105&&TARGET.HAS_CVMFS_fifeuser1_opensciencegrid_org==true&&TARGET.HAS_CVMFS_fifeuser2_opensciencegrid_org==true&&TARGET.HAS_CVMFS_fifeuser3_opensciencegrid_org==true&&TARGET.HAS_CVMFS_fifeuser4_opensciencegrid_org==true)' file:///dune/app/users/weishi/run_FDEffTarBall_grid.sh
```

At the top of ```run_FDEffTarBall_grid.sh```, you can set these variables:
```
Number of input files to run: STARTLINE, ENDLINE
Output directory: OUTDIR
```

Here are some reference settings:

300 events (3 file): ```--memory=502MB --disk=0.1GB --expected-lifetime=30m --cpu=1```

To check job status,

```
jobsub_q --user weishi
jobsub_q 57321037.0@jobsub01.fnal.gov -G dune
# For more options: jobsub_q --help
```

To fetch job output,

```
jobsub_fetchlog 57306776.0@jobsub01.fnal.gov  -G dune
```

Since Feb 15, 2023, [jobsub_lite](https://fifewiki.fnal.gov/wiki/Getting_started_with_jobsub_lite) becomes the new interface to submit and manage jobs.

### Locate files with SAM

To get a list of all files in DUNE dataset, use [SAM](https://dune.github.io/computing-training-basics/03-data-management/index.html):

```
source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
setup dunetpc v09_22_02 -q e19:debug
setup_fnal_security
setup sam_web_client
export SAM_EXPERIMENT=dune
samweb list-file-locations --dim="defname:prodgenie_nu_dune10kt_1x2x6_mcc11_lbl_reco" --schema=root -e dune
```

Or list 100 files in the dataset, do this at command line:

```
for k in `samweb list-files defname: prodgenie_nu_dune10kt_1x2x6_mcc11_lbl_reco with limit 100`; do samweb get-file-access-url --schema root $k; done
```

You can also find 'dimensions' of the sample from [MCC11](https://dune-data.fnal.gov/mc/mcc11/index.html). For example, the following is for the FD nutau reco sample

```
samweb list-files "file_type mc and dune.campaign mcc11 and data_tier 'full-reconstructed' and application reco and version v07_06_02 and file_name nutau_%"
```

This returns a list of files saved here: ```/dune/app/users/weishi/MCC11FDBeamsim_nutau_reco.txt```.

Then to locate full directory of each file:
```
samweb get-file-access-url nutau_dune10kt_1x2x6_12855916_0_20181104T221500_gen_g4_detsim_reco.root --schema=root
```
It returns the full url, e.g.,
```
root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/dune/tape_backed/dunepro/mcc11/protodune/mc/full-reconstructed/08/65/76/12/nutau_dune10kt_1x2x6_12855916_0_20181104T221500_gen_g4_detsim_reco.root
```

You can run ```lar``` directly on the url or open it via ROOT (no need to copy it):
```
lar -c <your fcl file>.fcl -n -1 root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/dune/tape_backed/dunepro/mcc11/protodune/mc/full-reconstructed/08/65/76/12/nutau_dune10kt_1x2x6_12855916_0_20181104T221500_gen_g4_detsim_reco.root

root -l root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/dune/tape_backed/dunepro/mcc11/protodune/mc/full-reconstructed/08/65/76/12/nutau_dune10kt_1x2x6_12855916_0_20181104T221500_gen_g4_detsim_reco.root
```

Another way to locate file:
```
samweb locate-file nutau_dune10kt_1x2x6_12855916_0_20181104T221500_gen_g4_detsim_reco.root
```
It tells you all the locations, e.g.,
```
enstore:/pnfs/dune/tape_backed/dunepro/mcc11/protodune/mc/full-reconstructed/08/65/76/12(14064@vr0423m8)
```

Other SAM functions: https://wiki.dunescience.org/wiki/DUNE_Computing/Main_resources_Jan2021#Data_management:_best_practices

```
source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
setup jobsub_client

# Need to prestage files from tape
samweb prestage-dataset --defname='prodgenie_nu_dune10kt_1x2x6_mcc11_lbl_reco' --parallel=6
```

Some samples are on tape, you can't access unless you prestage it first (may take time): https://dune.github.io/computing-training-basics/07-grid-job-submission/index.html

```
samweb prestage-dataset --defname='protodune-sp_runse6201_reco_v09_09_01_v0'
```

A better way to prestage is to instead do
```
unsetup curl # necessary as of May 2022 because there's an odd interaction with the UPS version of curl, so we need to turn it off
samweb run-project --defname=kherner-may2022tutorial-mc --schema https 'echo %fileurl && curl -L --cert $X509_USER_PROXY --key $X509_USER_PROXY --cacert $X509_USER_PROXY --capath /etc/grid-security/certificates -H "Range: bytes=0-3" %fileurl && echo'
```


## NN group machine environment setup

Note: this setup is not working anymore since ivy is replaced by nnhome (Debian bullseye). (Feb 19, 2023)

On NN group machine (ivy.physics.sunysb.edu, CentOS 6.10), I installed a DUNE software release on the ivy machine using the following setup, you can skip this part and go to [Set up work area on Ivy](#set-up-work-area-on-ivy).

```
#
# I've done this part, so you can skip it
#
mkdir ~/ups
mkdir ~/upstars
cd upstars
wget https://scisoft.fnal.gov/scisoft/bundles/tools/pullProducts
chmod +x ./pullProducts
./pullProducts ~/ups slf6 dune-v08_62_01 e19-py2 prof                                    # Scientific Linux 6
#
# I've done the above, start from settings below
#
```

### Set up work area on Ivy

[First time only]

```
mkdir ~/FDEff                                                                            # First time only
cd ~/FDEff
source /home/wshi/ups/setup      
setup git
setup gitflow
setup mrb
setup dunetpc v08_62_01 -q e19:py2:prof

export MRB_PROJECT=larsoft                                                               # Need to set ${MRB_PROJECT} to the master product
mrb newDev
source /home/<username>/FDEff/localProducts_larsoft_v08_62_01_e19_prof_py2/setup
# For example: source /home/wshi/FDEff/localProducts_larsoft_v08_62_01_e19_prof_py2/setup

cd srcs                            
git clone https://github.com/weishi10141993/myntuples.git                                # First time only, checkout the analysis code from GitHub
cp /home/wshi/dune/product_deps  /home/<username>/FDEff/srcs/myntuples/ups/product_deps  # Change version of larsoft, cetbuildtools, and qualifier list for dunetpc v08_62_01

mrb uc                                                                                   # Tell mrb to update CMakeLists.txt with the latest version numbers of the products.
cd ${MRB_BUILDDIR}                                                                       # Go to your build directory
mrb z
mrbsetenv                                                                                # Create the bookkeeping files needed to compile programs.
mrb install   
```

To run on FD MC files, this produces a TTree in myntuple.root in your work area:

```
cd /home/<your_username>/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis
# For example: cd /home/wshi/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis
lar -c MyEnergyAnalysis_ivy.fcl -n 10 -s /storage/shared/cvilela/DUNE_FD_MC/nu_dune10kt_1x2x6_13422341_0_20181123T225730_gen_g4_detsim_reco.root
```

The next time you login the ivy machine (username@ivy.physics.sunysb.edu), do the following to set up:

```
source /home/wshi/ups/setup
setup mrb
setup dunetpc v08_62_01 -q e19:py2:prof
source /home/<your_username>/FDEff/localProducts_larsoft_v08_62_01_e19_prof_py2/setup
# For example: source /home/wshi/FDEff/localProducts_larsoft_v08_62_01_e19_prof_py2/setup
mrbsetenv

# Go to your work directory and run your study!
cd /home/<your_username>/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis
# For example: cd /home/wshi/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis
```

### Notes for NN group machine (ivy) DUNE users

1. If you want to do "mrb g" on ANY machine, you need to have valid FNAL Kerberos ticket (kinit -f username@FNAL.GOV).

2. FNAL locks down http access to Redmine repositories via http offsite. You will need to be added as a developer so you can get read/write access (e.g., "mrb g dunetpc").

3. SLF6 is used on ivy now, but DUNE officially dropped SL6 builds. So you need to find an older builds and/or source code on SL6 machines. Make sure your builds of LArSoft (Linux64bit+2.6-2.12) match SLF6.

4. CVMFS is not installed on <oak, nngroup>.physics.sunysb.edu machines (Debian GNU/Linux 10 (buster)), you will need to manually install dunetpc releases (and dependencies) on the machine. Refer to this [LArSoft wiki](https://wiki.dunescience.org/wiki/DUNE_LAr_Software_Releases). I tried the instruction but so far no success yet on Debian GNU/Linux 10.

## File access
From DUNE FNAL machines:

* FD MC files: [FD Beamsim Requests](https://dune-data.fnal.gov/mc/mcc11/index.html)

* FD CAFs (no energy deposit details): ```/pnfs/dune/persistent/users/LBL_TDR/CAFs/v4/FD*```

* ND CAFs: ```/pnfs/dune/persistent/users/marshalc/nd_offaxis/v7/CAF```

From NN group machine:

* FD MC files: ```/storage/shared/cvilela/DUNE_FD_MC```

* On-axis ND CAFs to calculate the geometric efficiency correction for ND events: ```/storage/shared/cvilela/CAF/ND_v7```

## ND Geometry Efficiency

For reference, the ND analysis uses [these](https://github.com/DUNE/ND_CAFMaker) to produce CAF files (ntuples).

The ```dumptree.py``` file uses functions in this [repo](https://github.com/cvilelahep/DUNE_ND_GeoEff).

## Instruction for environment setup on SeaWulf

Refer to [Set up work area on SeaWulf](https://github.com/weishi10141993/myntuples/tree/ereco_study#set-up-work-area-on-seawulf).
