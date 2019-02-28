#! /bin/bash

# This script prepares the environment
# Allow read&write folders and files by default
umask 002

module purge
# Get simulation group modules
module use /work/ws/nemo/fr_lp1029-IMTEK_SIMULATION-0/modulefiles
# Get modules from Lukas
module use /work/ws/nemo/fr_le1019/modulefiles
module load smamp
module load gromacs/2016.4-gnu-5.2

# copy filestructure
cp -r .simulation_template 1_charge_cycle
