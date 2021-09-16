#!/bin/bash

DIRECTORY=FDEff
# we cannot rely on "whoami" in a grid job. We have no idea what the local username will be.
# Use the GRID_USER environment variable instead (set automatically by jobsub).
USERNAME=${GRID_USER}
source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
setup dunetpc v09_22_02 -q e19:debug # this also sets up the ifdh used in run_FDEffTarBall.sh
export WORKDIR=${_CONDOR_JOB_IWD} # if we use the RCDS the localProducts area will be placed in $CONDOR_DIR_INPUT
if [ ! -d "$WORKDIR" ]; then
  export WORKDIR=`echo ~`
fi

source ${INPUT_TAR_DIR_LOCAL}/${DIRECTORY}/localProducts*/setup-grid
mrbslp    # setup all products installed in localProducts*
