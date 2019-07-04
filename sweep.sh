#!/bin/bash
source ~/.bashrc

# Calculate cost functions with HORTON
module purge
module load horton/2.1.0b3 2>/dev/null
module load smamp
python sweep_rhoref.py

# Average over cost functions
module purge
module load horton/2.1.0b3 2>/dev/null
module load smamp
python sweep_average.py 

# Fit ESP cost function
module purge
module load horton/2.1.0b3 2>/dev/null
module load smamp
module load gromacs
python sweep_charges.py

# Plot charges
module purge
module load smamp
module load devel/python/3.6.5
python sweep_2d_plot.py
python sweep_plot.py

python collect_charges.py
python plot_charges.py
