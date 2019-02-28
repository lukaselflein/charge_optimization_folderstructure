# Trajectory Analysis

This is a template for the mid-level folder structure for charge optimization.

### Usage 
To create the full structure, copy the `snapshot_template` folder structure n times for your n snapshots, and give them distinct names (e.g., `100ps`, `200ps`, etc.).
Copy the n snapshots into the corresponding subfolder, e.g.,`100ps/intial_structure/`.
Then, run the mid-level bash script `trajectory_analysis.sh`.

### Content
`snapshot_template`: universal folder structure for fitting charges to snapshot structure
`{100}ps`: a copy of `snapshot_template` filled with concrete data
`charge_comparison`: lists and visualizes charges from the different snapshots
`analyze_snapshots.sh`: runs analytics scripts in subfolders, compares charges
`create_snapshots.py`: initialize empty folder structures from template
