# DUNE FD Geometry Efficiency

## Instruction for environment setup from DUNE FNAL machines (dunegpvm*)

[First time only]

```
cd /dune/app/users/<your_username>
mkdir FDEff (first time only)
cd FDEff

source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
setup dunetpc v09_22_02 -q e19:debug
[optional] setup_fnal_security                            # A FNAL grid proxy to submit jobs and access data in dCache via xrootd or ifdh.

mrb newDev
source /dune/app/users/<your_username>/inspect/localProducts_larsoft_${LARSOFT_VERSION}_debug_${COMPILER}/setup
#for example: source /dune/app/users/weishi/FDEff/localProducts_larsoft_v09_22_02_debug_e19/setup

cd srcs
git clone https://github.com/weishi10141993/myntuples.git # First time only, checkout the analysis code from GitHub

mrb uc                                                    # Tell mrb to update CMakeLists.txt with the latest version numbers of the products.
cd ${MRB_BUILDDIR}                                        # Go to your build directory
mrb z
mrbsetenv                                                 # Create the bookkeeping files needed to compile programs.
mrb install                                               # Compile the code in ${MRB_SOURCE} and put the results in ${MRB_INSTALL}
```

To run on FD MC files, this produces a ROOT nTuple:

```
cd /dune/app/users/<your_username>/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis
#for example: cd /dune/app/users/weishi/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis
lar -c MyEnergyAnalysis.fcl -n -1
```

The next time you login a DUNE FNAL machine (dunegpvm*), do the following to set up:

```
source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
setup dunetpc v09_22_02 -q e19:debug
source /dune/app/users/<your_username>/inspect/localProducts_larsoft_${LARSOFT_VERSION}_debug_${COMPILER}/setup
#for example: source /dune/app/users/weishi/FDEff/localProducts_larsoft_v09_22_02_debug_e19/setup
mrbsetenv
cd /dune/app/users/<your_username>/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis
```

If changed ```MyEnergyAnalysis_module.cc```, recompile the code (do the setup above first):

```
cd ${MRB_BUILDDIR}                        # Go to your build directory
mrb z                                     # Remove old build directory
mrbsetenv                                 # Create the bookkeeping files needed to compile programs.
mrb install                               # Compile the code in ${MRB_SOURCE} and put the results in ${MRB_INSTALL}
```

If added new package in ```srcs``` directory, do ```mrb uc``` and then recompile as above.

To commit changed code changes to remote repository:

```
git commit
git push
```

## Instruction for environment setup from NN group machine: ivy.physics.sunysb.edu (CentOS 6.10)

I've manually installed a DUNE software release on the ivy.physics.sunysb.edu (CentOS 6.10) machine, so you can skip to: # Set up your work directory.

```
# I've done this part and installed a local release on ivy, so you can skip to: # Set up your work directory
mkdir ~/ups
mkdir ~/upstars
cd upstars
wget https://scisoft.fnal.gov/scisoft/bundles/tools/pullProducts
chmod +x ./pullProducts
./pullProducts ~/ups slf6 dune-v08_62_01 e19-py2 prof                                    # Scientific Linux 6
# I've done the above, start from settings below

# Set up your work directory
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
source /home/wshi/FDEff/localProducts_larsoft_v08_62_01_e19_prof_py2/setup
mrbsetenv

# Go to your work directory and run your study!
cd /home/<your_username>/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis
# For example: cd /home/wshi/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis
```

## Notes/tips for NN group machine (ivy) DUNE users:

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
