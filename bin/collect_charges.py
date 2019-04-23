""" Extract charges obtained via HORTON and Bader.
Copyright 2019 Simulation Lab
University of Freiburg
Author: Lukas Elflein <elfleinl@cs.uni-freiburg.de>
"""

import os
import pandas as pd


def create_dir(path='./plotting'):
	"""Create new folder for pictures if it does not exist yet."""
	if os.path.isdir(path):
		return path

	os.makedirs(path)
	return path


def collect_bader():
	"""Find charges and put them in one dataframe."""
	# Initialize collection data frame
	coll_df = None
	# Crawl the directory structure
	for subdir, dirs, files in sorted(os.walk('./')):

		# Exclude template folders from search
		if 'template' in subdir or 'exclude' in subdir:
			continue

		# Select the folders with cost function
		if 'bader_charges' in subdir:
			print('Moving to {}'.format(subdir))
			# Extract timestamp
			time = os.path.split(subdir)[0].replace('./', '').replace('_ps_snapshot', '')
			time = int(time)
		
			# Use the first charge file to come across as a template	
			df = pd.read_csv(os.path.join(subdir, 'bader_charges.csv'), sep=r',\s*',
					 engine='python')
			df['timestamp'] = time

			if coll_df is None:
				coll_df = df
			else:
				coll_df = coll_df.append(df)

	# The table still contains redundant hydrogen atoms: 1CD3... 2CB3
	# Delete everything containing '1C' or '2C'
	# print(coll_df[coll_df.atom.str.contains(r'[1-2]C')])
	coll_df = coll_df.drop(coll_df[coll_df.atom.str.contains(r'[1-2]C')].index)

	print('All collected. Transforming wide to long format ...')
	# Transform the wide format into a long format version (for easier plotting)
	coll_df = coll_df.rename({'q': 'bader'}, axis=1)
	coll_df = pd.melt(coll_df, id_vars=['atom', 'residue', 'timestamp'], 
			  value_vars=['bader'])
	coll_df = coll_df.rename({'value': 'charge', 'variable': 'Calculation Variant'},
				 axis=1)
	return coll_df


def collect_horton():
	"""Find charges and put them in one dataframe."""
	# Initialize collection data frame
	coll_df = None
	# Crawl the directory structure
	for subdir, dirs, files in sorted(os.walk('./')):

		# Exclude template folders from search
		if 'template' in subdir or 'exclude' in subdir:
			continue

		# Select the folders with cost function
		if 'horton_cost_function' in subdir:
			print('Moving to {}'.format(subdir))
			# Extract timestamp
			time = os.path.split(subdir)[0].replace('./', '').replace('_ps_snapshot', '')
			time = int(time)
		
			# Use the first charge file to come across as a template	
			df = pd.read_csv(os.path.join(subdir, 'fitted_point_charges.csv'))
			df['timestamp'] = time

			if coll_df is None:
				coll_df = df
			else:
				coll_df = coll_df.append(df)

	print('All collected. Transforming wide to long format ...')
	# Transform the wide format into a long format version (for easier plotting)
	coll_df = coll_df.rename({'q': 'constrained', 'q_unconstrained': 'unconstrained'}, axis=1)
	coll_df = pd.melt(coll_df, id_vars=['atom', 'residue', 'timestamp'], 
			  value_vars=['constrained', 'unconstrained'])
	coll_df = coll_df.rename({'value': 'charge', 'variable': 'Calculation Variant'}, axis=1)
	return coll_df


def extract_init_charges(rtp_path, df):
        """Extract charges from rtp file"""
        atom_names = df.atom.unique()
        residuum_names = df.residue.unique()
        charges = pd.DataFrame()
        with open(rtp_path, 'r') as rtp_file:
                print('Successfully loaded topolgy file {}'.format(rtp_path))
                rtp_text = rtp_file.readlines()

                for line in rtp_text:
                        current_residuum = None
                        # atom names are only unique inside one residuum
                        # Thus, specify which res we are currently in
                        for residuum in residuum_names:
                                if residuum in line:
                                        current_residuum = residuum
                        # Now, we can look up the atom name in the charge table.
                        # First, select the lines with exactly one atom name
                        for atom_name in atom_names:
                                # Select lines with at least one atom name
                                if atom_name in line[0:7]:
                                        second_entry = line[8:18].replace('+', '')
                                        second_entry = second_entry.replace('-', '').strip()
                                        # Select lines with no atom name in second column
                                        if not second_entry in atom_names:
                                                q_value = float(line[24:34].strip(' '))
                                                charges = charges.append({'atom': atom_name,
                                                                          'residue': current_residuum,
                                                                          'q_init': q_value},
                                                                          ignore_index=True)
        return charges


def collect_average():
	"""Put averaged charegs in a dataframe."""
	# Read charges from averaged cost function
	input_path = './horton_charges/fitted_point_charges.csv'
	avg_df = pd.read_csv(input_path)
	# Rename columns for consistency
	avg_df = avg_df.rename({'q': 'averaged cost function'}, axis=1)
	# Transform to long format
	avg_df = pd.melt(avg_df, id_vars=['atom', 'residue'], value_vars=['averaged cost function'])
	avg_df = avg_df.rename({'value': 'charge', 'variable': 'Calculation Variant'}, axis=1)
	return avg_df


def main():
	"""Collect charges and save them to .csv file"""
	# Collect averaged charges
	avg_df = collect_average()

	# Collect all horton charges
	print('Collecting HORTON charges ...')
	horton_df = collect_horton()

	# Collect all bader charges
	print('Collecting Bader charges ...')
	print('Skipping, Bader is broken right now.')
	# bader_df = collect_bader()

	# Paste everything into single dataframe
	print('Combining different charges into one table ... ')
	constr_df = horton_df.loc[horton_df['Calculation Variant'] == 'constrained']
	unconstr_df = horton_df.loc[horton_df['Calculation Variant'] == 'unconstrained']
	collect_df = avg_df
	collect_df = collect_df.append(constr_df, sort=False) 
	collect_df = collect_df.append(unconstr_df, sort=False) 
	# collect_df = collect_df.append(bader_df, sort=False)

	create_dir(path='./plotting')
	collect_df.to_csv('./plotting/all_charges.csv')

if __name__ == '__main__':
	main()
