# DUNE FD Geometry Efficiency

## Standalone LArSoft Module for FD Geometry Efficiency

From DUNE FNAL machines (dunegpvm*):

```
cd /dune/app/users/<your_username>
mkdir FDEff (first time only)
cd FDEff

source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
setup dunetpc v09_22_02 -q e19:debug
setup_fnal_security

mrb newDev
source /dune/app/users/<your_username>/inspect/localProducts_larsoft_${LARSOFT_VERSION}_debug_${COMPILER}/setup
#for example: source /dune/app/users/weishi/FDEff/localProducts_larsoft_v09_22_02_debug_e19/setup

cd srcs
git clone https://github.com/weishi10141993/myntuples.git # First time only, checkout the analysis code from GitHub

mrb uc                                    # Tell mrb to update CMakeLists.txt with the latest version numbers of the products.
cd ${MRB_BUILDDIR}                        # Go to your build directory
mrb z
mrbsetenv                                 # Create the bookkeeping files needed to compile programs.
mrb install                               # Compile the code in ${MRB_SOURCE} and put the results in ${MRB_INSTALL}
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
setup_fnal_security
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

Fyi, a good instruction on writing a [LArSoft module](https://cdcvs.fnal.gov/redmine/projects/larsoft/wiki/_AnalysisExample_).

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
