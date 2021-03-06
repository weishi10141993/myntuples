# Produce Ntuple from DUNE FD MC files

The following instruction is used to produce ROOT n-tuples from FD MC files (mcc11): [FD Beamsim Requests](https://dune-data.fnal.gov/mc/mcc11/index.html).

## DUNE FNAL machines (dunegpvm*) environment setup

[First time only]

```
cd /dune/app/users/weishi                                               # Replace with your username for all commands below
mkdir FDEff (first time only)
cd FDEff

source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
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

Once the above is compiled and runs without problem interactively, you can start to produce a tarball. First, you need to have a grid setup for the localProducts as the grid job typically runs on a different machine than your working machine,

```
cd /dune/app/users/weishi/FDEff/localProducts_larsoft_v09_22_02_debug_e19
cp setup setup-grid         # make a copy of the setup for grid job
```

then change ```/dune/app/users/weishi``` to ```${INPUT_TAR_DIR_LOCAL}```.

```
cd /dune/app/users/weishi
# Get job set up scripts
wget https://raw.githubusercontent.com/weishi10141993/NeutrinoPhysics/main/setupFDEffTarBall-grid.sh --no-check-certificate
wget https://raw.githubusercontent.com/weishi10141993/NeutrinoPhysics/main/run_FDEffTarBall.sh --no-check-certificate

# Now make the tarball
tar -czvf FDEff.tar.gz FDEff setupFDEffTarBall-grid.sh

# Check the tarball *.tar.gz is indeed created and open with: tar -xf *.tar.gz
```

```
jobsub_submit -G dune -M -N 1 --memory=500MB --disk=0.1GB --expected-lifetime=10m --cpu=1 --resource-provides=usage_model=DEDICATED,OPPORTUNISTIC,OFFSITE --tar_file_name=dropbox:///dune/app/users/weishi/FDEff.tar.gz --use-cvmfs-dropbox -l '+SingularityImage=\"/cvmfs/singularity.opensciencegrid.org/fermilab/fnal-wn-sl7:latest\"' --append_condor_requirements='(TARGET.HAS_Singularity==true&&TARGET.HAS_CVMFS_dune_opensciencegrid_org==true&&TARGET.HAS_CVMFS_larsoft_opensciencegrid_org==true&&TARGET.CVMFS_dune_opensciencegrid_org_REVISION>=1105&&TARGET.HAS_CVMFS_fifeuser1_opensciencegrid_org==true&&TARGET.HAS_CVMFS_fifeuser2_opensciencegrid_org==true&&TARGET.HAS_CVMFS_fifeuser3_opensciencegrid_org==true&&TARGET.HAS_CVMFS_fifeuser4_opensciencegrid_org==true)' file:///dune/app/users/weishi/run_FDEffTarBall.sh
```

Here are some reference settings:

100 events (1 file): ```--memory=338.3MB --disk=0.0GB --expected-lifetime=1m21s --cpu=1```

10k events (100 files): ```--memory=?MB --disk=?GB --expected-lifetime=?h --cpu=?```

Default: ```--memory=2000MB --disk=10GB --expected-lifetime=8h --cpu=1```.

To check job status,

```
jobsub_q --user weishi
# For more options: jobsub_q --help
```

To fetch job output,

```
jobsub_fetchlog --jobid=<id> --unzipdir=<dir>
```

The output is in the scratch area ```/pnfs/dune/scratch/users/weishi/myFDntuples``` as specified in ```run_FDEffTarBall.sh```.

```
source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
setup jobsub_client

# Need to prestage files from tape
samweb prestage-dataset --defname='prodgenie_nu_dune10kt_1x2x6_mcc11_lbl_reco' --parallel=6
```

This instruction is based on the [DUNE computing tutorial](https://wiki.dunescience.org/wiki/DUNE_Computing/Submitting_grid_jobs_May2021#Submit_a_job).

### Use SAM to navigate files

To get a list of all files in DUNE dataset, use SAM:

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

Other SAM functions: https://wiki.dunescience.org/wiki/DUNE_Computing/Main_resources_Jan2021#Data_management:_best_practices


## NN group machine environment setup:

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

## [!!! Under construction !!!] Instruction for environment setup from SeaWulf cluster: seawulf.stonybrook.edu (CentOS 7.9.2009)

Use the following setup to first install a DUNE software release on the SeaWulf machine. Note files older than 1 month will be purged from the SeaWulf scratch directory, so you may want to ask for more quota in your SeaWulf home area and do the following software installation (need more than 20 GB) there. After that, go to [Install CMake binaries](#install-cmake-binaries) and [Set up work area on SeaWulf](#set-up-work-area-on-SeaWulf).

```
cd /gpfs/scratch/<netID>                                                                                                    
mkdir ups
mkdir upstars
cd upstars
wget https://scisoft.fnal.gov/scisoft/bundles/tools/pullProducts
chmod +x ./pullProducts
./pullProducts /gpfs/scratch/<netID>/ups slf7 dune-v09_22_02 e19 debug
```

### Install CMake binaries

```
mkdir ~/CMake
cd ~/CMake
wget https://github.com/Kitware/CMake/releases/download/v3.13.2/cmake-3.13.2-Linux-x86_64.tar.gz
tar xf cmake-3.13.2-Linux-x86_64.tar.gz
export PATH="/gpfs/home/<netID>/CMake/cmake-3.13.2-Linux-x86_64/bin:$PATH"               # Save it in .bashrc if needed
```

### Set up work area on SeaWulf

This setup doesn't work yet. ```mrbsetenv``` reports ```ERROR: Action parsing failed on "unsetuprequired(cmake v3_13_2)"```.

[First time only]

```
mkdir ~/FDEff                                                                            # First time only
cd ~/FDEff
source /gpfs/scratch/<netID>/ups/setup                                                   # You can use mine as I've installed it: source /gpfs/scratch/weishi2/ups/setup
setup git
setup gitflow
setup mrb
setup dunetpc v09_22_02 -q e19:debug

export MRB_PROJECT=larsoft                                                               # Need to set ${MRB_PROJECT} to the master product
mrb newDev
source /gpfs/home/<netID>/FDEff/localProducts_larsoft_v08_62_01_e19_prof_py2/setup
# For example: source /gpfs/home/weishi2/FDEff/localProducts_larsoft_v09_22_02_debug_e19/setup

cd srcs                            
git clone https://github.com/weishi10141993/myntuples.git                                # First time only, checkout the analysis code from GitHub

mrb uc                                                                                   # Tell mrb to update CMakeLists.txt with the latest version numbers of the products.
cd ${MRB_BUILDDIR}                                                                       # Go to your build directory
mrb z
mrbsetenv                                                                               
mrb install   
```
