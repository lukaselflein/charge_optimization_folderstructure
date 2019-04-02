# Binaries and Scripts

A collection of all code that is reused in the workflow.

Ideally, only this code is referenced in the workflow, either by calling it via a relative path (e.g., "python3 ../../example.py") or by creating symbolic links to the scripts (this might be fragile!).
This way, updates to any of the code will be immediately distributed to all subfolders.

### Content
* `aa2ua_cube.py`: Maps point charges obtained by GPAW and HORTON on the original GROMACS topology.
* `average_cost.py`: Average Cost Functions for Horton to determine Charges for Molecular Dynamics.
* `charges_to_rtp.py`: Transfer Charges from CSV table to .rft file.
* `convert_UA_to_AA.py`: Change structure with implicit Hydrogen to explicitely defined H-atoms.
* `collect_charges.py`: Extract charges from output files and write them into a long-format table.
* `create_snapshots_from_trajectory.py`: Create snapshots from trajectory, initialize folder structure 
* `fitESPconstrained.py`: Fits (united-atom) point charges onto (all-atom) ESP obtained by 
* `gpaw_optimize_and_esp.py`: Minimize the enegy of input structure and extract ESP. 
* `loop_bader.py`: Search folderstructure for DFT output files, calculate bader charges.
* `loop_charges.py`: Calculate HORTON charges for all .cube files in folderstructure. 
* `loop_convert_AA_to_UA.py`: Search for 'all atoms' structures used in DFT, convert them to 'united atoms' format.
* `loop_convert_UA_to_AA.py`: Search for inital united atom structures, convert them to all atom format.
* `loop_cost_functions.py`:  Find files and calculate cost functions.
* `loop_submit.py`:  Search for inital united atom structures, convert them to all atom format.
* `plot_charges.py`: Plot charges with different constraint options from table.
* `submit_gpaw.sh`: NEMO queue submission file for DFT calculations .

