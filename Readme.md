# Sarah Folder Structure
Create a consistent, documented pipline and folder structure for self-consistent MD charge optimization.

### Goal
SMAMP (synthetic mimics of antimicrobial peptides) can be used for antimicrobial coating in medical applications.
To understand a class of SMAMP, we simulate them via Molecular Dynamics (MD).
An essential component of MD are point-charges.
These are not available for the SMAMP molecules of interest, so we have to determine them ourselves.
This repo contains scripts and a folderstructure to automate the calculation of these charges.

### Prerequisites
Make sure you have the smamp module:
1. Local installation
```bash
user@machine:~$ pip install --user smamp
```

2. Module installation
Alternatively, you can install it as a module:
```bash
user@machine:~$ mkdir ~/modulefiles
user@machine:~$ cd ~/modulefiles
user@machine:~$ git clone https://github.com/lukaselflein/smamp; cd smamp
user@machine:~$ echo "#%Module1.0" > ~/modulefiles/smamp
user@machine:~$ echo "prepend-path PYTHONPATH $(pwd)" >> ~/modulefiles/smamp
```

Remember to add it to and load from your modulefiles:
```bash
user@machine:~$ module use ~/modulefiles
user@machine:~$ module load smamp
```
3. Use Lukas' local version
```bash
user@machine:~$ module use /home/fr/fr_fr/fr_le1019/modulefiles
user@machine:~$ module load smamp
```

### Organization
The folders are organized in four levels:
1. Lowerest: subfolders for the diverse charge fitting and conversion routines.
2. Snapshot: Each snapshot has its own folder, containing (1).
3. Simulation: For each set of charges, a seperate MD simulation is run. Contains subfolders for the snapshots extracted from the trajectory (2).
4. Current directory: Simulating and fitting charges is iteratated until charges converge. This is the highest level. Subfolders correspond to different simulations (3).
<img src="./.pictures/folder_hierarchy.png" width="400px">

### Content
`Readme.md`: The Readme you are reading right now.
`bin`: Bash and python scripts for fitting, conversion, extracting, and visualization.
`preprocessing.sh`: BASH Commands for preprocessing the input up to and including the DFT calculations.
`postprocessing.sh`: BASH commands to extract densities, calculate and visualize charges.
`.pictures`: Pictures for the Readme.
`.simulation_template`: The full template simulation folder structure.

### Usage
Copy the `simulation_template` folder structure, and rename to a specific simulation name, e.g. `simulation-1.0-a`. 

Then, run the simulations, make sure the data is saved in the appropriate location.
