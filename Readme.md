# Charge Optimization Folderstructure
Create a consistent, documented pipline and folder structure for self-consistent MD charge optimization.

### Goal
SMAMP (synthetic mimics of antimicrobial peptides) can be used for antimicrobial coating in medical applications.
To understand a SMAMP better, we simulate it via Molecular Dynamics (MD).
An essential component of MD are point-charges.
These are not available for the SMAMP molecules of interest, so we have to determine them ourselves.
There are different ways in which we can calculate these charges, and it is not clear which one is the best.
A neccessary criterion for good charges is self-consistency: 
We calculate them from some snapshots of a MD simulation, an then do a new MD simulation with these new charges q*. 
We then re-do the charge calculation, and if the resulting charges q** are similar to q*, this is a good indicator that our charges are meaningful.

This repo contains scripts and a folderstructure to automate the calculation workflow to get these charges.

### Prerequisites
Some tools and often-reused snippets are outsourced into a pip package, [smamp](https://github.com/lukaselflein/smamp).
Make sure you have the smamp module:
1. Local installation
```bash
pip install --user smamp
```

2. Alternatively, you can install it as a module:
```bash
mkdir ~/modulefiles
cd ~/modulefiles
git clone https://github.com/lukaselflein/smamp; cd smamp
echo "#%Module1.0" > ~/modulefiles/smamp
echo "prepend-path PYTHONPATH $(pwd)" >> ~/modulefiles/smamp
```

Remember to add it to and load from your modulefiles:
```bash
module use ~/modulefiles
module load smamp
```
3. Use Lukas' local version
```bash
module use /home/fr/fr_fr/fr_le1019/modulefiles
module load smamp
```

### Organization
The Project is organized in multiple levels:
1. Highest: This is the current directory. Here, you can put the different self-consistent charge cycles, e.g. the first, second, .. iterations. Also, the folder with the calculation scripts `bin` is located on this level.
2. Charge cycle: Here you will find the administrative scripts, for looping over the actual calculations. Also, the different snapshots are in seperate directories on this level.
3. Snapshot: For every timestamp, a new directory is created. It contains folders for the different steps of the charge optimization workflow.
<img src="./.pictures/folder_hierarchy.png" width="400px">

### Usage Example
If you have done all the prerequisite steps, you can start an automated charge optimization.
#### Clone a fresh copy from github
```bash
git clone https://github.com/lukaselflein/charge_optimization_folderstructure
mv charge_optimization_folderstructure/ example/
```
It is probably a good idea to name this folder according to what you want to archieve, e.g., `example`, `convergence study`, etc.

#### Create your first charge cycle
Go into your fresh copy. Inside, copy the template folder:
```bash
cp -r .template_simulation 1_charge_cycle 
cd 1_charge_cycle
```
Now you can work on this first cycle. If you repeat this step after the first cycle is finished, for the second ... iteration `2_charge_cycle`, you have an overview of all charge cycles side-by-side.

#### Copy MD simulation files
We need the trajectory, etc. in our optimization folderstructure. Copy them into `md_simulation`:
```bash
mkdir md_simulation
cp PATH_TO_MD_SIMULATION/example.top md_simulation/
cp PATH_TO_MD_SIMULATION/example.tpr md_simulation/
cp PATH_TO_MD_SIMULATION/example.xtc md_simulation/
p PATH_TO_MD_SIMULATION/example.rtp md_simulation/
```

#### Create the snapshots and subfolder structure
If we want snapshots at time = 100 ps, 200 ps, ... 1000 ps:
```bash
module purge
module load gromacs/2016.4-gnu-5.2
python create_snapshots_from_trajectory.py -tpr md_simulation/example.tpr -top md_simulation/example.top -xtc md_simulation/example.xtc -s 100 -d 100 -e 1000
```

#### Convert UA to AA
```bash
module purge
module load gromacs/2016.4-gnu-5.2
module load matscipy/0.2.0
module load smamp
python loop_convert_UA_to_AA.py
```

#### Submit gpaw optimization, ESP & Rho job
```bash
module purge
module load devel/python/3.6.5
module load smamp
python loop_submit.py
```
Now you will have to wait for the DFT calculations to finish.

### Automation
If you want to automate all of the above steps, you can edit the `preprocessing.sh` script. All the above commands are in there, and you can write your filenames and paths in there.

### Postprocessing
The bash script `postprocessing.sh` should work out of the box, and calculate charges from the converged DFT calculations.
Plots are automatically generated in `images`.

### Content
* `Readme.md`: The Readme you are reading right now.
* `bin`: Bash and python scripts for fitting, conversion, extracting, and visualization.
* `preprocessing.sh`: BASH Commands for preprocessing the input up to and including the DFT calculations.
* `postprocessing.sh`: BASH commands to extract densities, calculate and visualize charges.
* `.pictures`: Pictures for the Readme.
* `.simulation_template`: The full template simulation folder structure.
