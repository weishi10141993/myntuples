# Produce Ntuple from DUNE FD MC files

The following instruction is used to produce ROOT n-tuples from FD MC files (mcc11): [FD Beamsim Requests](https://dune-data.fnal.gov/mc/mcc11/index.html).

## Instruction for environment setup from SeaWulf cluster: seawulf.stonybrook.edu (CentOS 7.9.2009)

The following setup installs DUNE softwares on the SeaWulf machine. After this, go to [Set up work area on SeaWulf](#set-up-work-area-on-SeaWulf).

```
#
# I've done this part, so you can skip it
#
cd /gpfs/projects/McGrewGroup/weishi/DUNE
mkdir ups
mkdir upstars
cd upstars
wget https://scisoft.fnal.gov/scisoft/bundles/tools/pullProducts
chmod +x ./pullProducts
./pullProducts /gpfs/projects/McGrewGroup/weishi/DUNE/ups slf7 dune-v09_22_02 e19 debug
```

## Set up work area on SeaWulf

[First time only]

```
cd /gpfs/projects/McGrewGroup/<usrname>/DUNE
mkdir FDEff                                                                              # First time only
cd FDEff

module unload cmake
module load cmake/3.17.3
# unsetup cmake v3_13_2

source /gpfs/projects/McGrewGroup/weishi/DUNE/ups/setup                                  # source this file under my directory                              
setup git
setup gitflow
setup mrb
setup dunetpc v09_22_02 -q e19:debug

export MRB_PROJECT=larsoft                                                               # Need to set ${MRB_PROJECT} to the master product
mrb newDev
source /gpfs/scratch/weishi2/FDEff/localProducts_larsoft_v09_22_02_debug_e19/setup

cd srcs                            
git clone https://github.com/weishi10141993/myntuples.git -b ereco_study                 # First time only, checkout the analysis code from GitHub

mrb uc                                                                                   # Tell mrb to update CMakeLists.txt with the latest version numbers of the products.
cd ${MRB_BUILDDIR}                                                                       # Go to your build directory
mrb z
mrbsetenv                                                                               
mrb install   
```

To run on FD MC files, this produces a ```TTree``` in ```myntuple.root``` in your work area:

```
cd /gpfs/projects/McGrewGroup/<usrname>/DUNE/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis
lar -c MyEnergyAnalysis.fcl -n 10
```

The next time you login, do the following to set up:

```
source /gpfs/projects/McGrewGroup/weishi/DUNE/ups/setup
setup mrb
setup dunetpc v09_22_02 -q e19:debug
source /gpfs/projects/McGrewGroup/<usrname>/DUNE/FDEff/localProducts_larsoft_v09_22_02_debug_e19/setup
mrbsetenv

# Go to your work directory and run your study!
cd /gpfs/projects/McGrewGroup/<usrname>/DUNE/FDEff/srcs/myntuples/myntuples/MyEnergyAnalysis
```

If the source code ```MyEnergyAnalysis_module.cc``` is changed, recompile the code (do the setup above first):

```
cd ${MRB_BUILDDIR}                        # Go to your build directory
mrb z                                     # Remove old build directory
mrbsetenv                                 # Create the bookkeeping files needed to compile programs
mrb install                               # Compile the code in ${MRB_SOURCE} and put the results in ${MRB_INSTALL}
```

If added new package in ```srcs``` directory, do ```mrb uc``` and then recompile as above.

To commit changed code changes to remote repository:

```
git commit
git push
```
