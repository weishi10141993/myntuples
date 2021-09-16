#!/bin/bash

echo "Running on $(hostname) at ${GLIDEIN_Site}. GLIDEIN_DUNESite = ${GLIDEIN_DUNESite}"

# Set the output location for copyback
OUTDIR=/pnfs/dune/scratch/users/${GRID_USER}/myFDntuples

# Let's rename the output file so it's unique in case we send multiple jobs.
OUTFILE=myntuple_${CLUSTER}_${PROCESS}_$(date -u +%Y%m%dT%H%M%SZ).root

# Make sure we see what we expect
pwd

# Tarball is copied and untarred into a directory on the worker node, accessed via this CONDOR_DIR_INPUT environment variable
ls -l $CONDOR_DIR_INPUT

echo "CONDOR_DIR_INPUT: ${CONDOR_DIR_INPUT}"
echo "INPUT_TAR_DIR_LOCAL: ${INPUT_TAR_DIR_LOCAL}"

if [ -e ${INPUT_TAR_DIR_LOCAL}/${DIRECTORY}/srcs/myntuples/setupFDEffTarBall-grid.sh ]; then
  . ${INPUT_TAR_DIR_LOCAL}/${DIRECTORY}/srcs/myntuples/setupFDEffTarBall-grid.sh
else
  echo "Error, setup script not found. Exiting."
  exit 1
fi

# Go back to the top-level directory since we know that's writable
cd ${_CONDOR_JOB_IWD}

echo "Look at _CONDOR_JOB_IWD"
ls
echo "And _CONDOR_DIR_INPUT"
ls ${_CONDOR_DIR_INPUT}

# Symlink the desired fcl to the current directory
ln -s ${INPUT_TAR_DIR_LOCAL}/${DIRECTORY}/srcs/myntuples/myntuples/MyEnergyAnalysis/MyEnergyAnalysis.fcl .

# Set some other very useful environment variables for xrootd and IFDH
export IFDH_CP_MAXRETRIES=2
export XRD_CONNECTIONRETRY=32
export XRD_REQUESTTIMEOUT=14400
export XRD_REDIRECTLIMIT=255
export XRD_LOADBALANCERTTL=7200
export XRD_STREAMTIMEOUT=14400 # may vary for your job/file type

# Make sure the output directory exists
ifdh ls $OUTDIR 0 # set recursion depth to 0 since we are only checking for the directory; we don't care about the full listing.

if [ $? -ne 0 ]; then
  # if ifdh ls failed, try to make the directory
  ifdh mkdir_p $OUTDIR || { echo "Error creating or checking $OUTDIR"; exit 2; }
fi

echo "Finished checking outdir: $OUTDIR"

# Get the xrootd URL for the input file. Not necessary for SAM inputs when using ifdh_art, etc.
myinfile=$(samweb get-file-access-url --schema=root nu_dune10kt_1x2x6_13009312_0_20181104T221530_gen_g4_detsim_reco.root)

echo "Got xrootd url: $myinfile"

# Now we should be in the work dir if setupFDEffTarBall-grid.sh worked
lar -c MyEnergyAnalysis.fcl -n -1 $myinfile
LAR_RESULT=$?   # check the exit status!!!

if [ $LAR_RESULT -ne 0 ]; then
  echo "lar exited with abnormal status $LAR_RESULT. See error outputs."
  exit $LAR_RESULT
fi

if [ -f myntuple.root ]; then

  mv myntuple.root $OUTFILE

  # and copy our output file back
  ifdh cp -D $OUTFILE $OUTDIR

  # check the exit status to see if the copyback actually worked. Print a message if it did not.
  IFDH_RESULT=$?
  if [ $IFDH_RESULT -ne 0 ]; then
    echo "Error during output copyback. See output logs."
    exit $IFDH_RESULT
  fi
fi

#If we got this far, we succeeded.
echo "Completed successfully."
exit 0
