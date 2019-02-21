#! /bin/bash

# This script prepares the environment
# Allow read&write folders and files by default
umask 002

# Load default modules
module purge
module use /work/ws/nemo/fr_lp1029-IMTEK_SIMULATION-0/modulefiles
module load gromacs/2016.4-gnu-5.2
module load devel/python/3.6.5

# copy filestructure
# cp -r simulation_template example
