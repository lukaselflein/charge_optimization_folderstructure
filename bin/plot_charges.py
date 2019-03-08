""" Extract and plot charges obtained via HORTON.
Copyright 2019 Simulation Lab
University of Freiburg
Author: Lukas Elflein <elfleinl@cs.uni-freiburg.de>
"""

import os
import sys
import shutil
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def create_dir(path='./pictures'):
	if os.path.isdir(path):
		return path

	os.makedirs(path)
	return path


def extract_init_charges(rtp_path, df):
        atom_names = df.atom.unique()
        residuum_names = df.residue.unique()
        charges = pd.DataFrame()
        with open(rtp_path, 'r') as rtp_file:
                print('Successfully loaded topolgy file {}'.format(rtp_path))
                rtp_text = rtp_file.readlines()

                current_residuum = None
                for line in rtp_text:
                        # Atom names are only unique inside one residuum
                        # Thus, specify which res we are currently in
                        for residuum in residuum_names:
                                if residuum in line:
                                        current_residuum = residuum
                                        break
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


def generic_scatterplot(df, x_name, y_name):
	"""Default plotting options"""
	fgrid = sns.lmplot(x=x_name, y=y_name, data=df, fit_reg=True)

	# Add annotations
	ax = fgrid.axes[0,0]
	offset = df[x_name].max() * 0.05
	for line in range(0,df.shape[0]):
		text = df['atom'][line]
		x = df[x_name][line]
		y = df[y_name][line]

		# Annotate outliers (far from before=after) only, skip rest
		if abs(x - y) < 0.2:
			continue

		# Add text to the right of the dot
		ax.text(x + offset, y, text,
			horizontalalignment='left', size='medium', 
			color='black', weight='semibold')

	# Add reference function x=y
	sns.lineplot(x=df[x_name], y=df[x_name], ax=ax, color='grey')
	return fgrid, ax

def scatterplot_comparison(df, out_path, x_name='q', y_name='q_init'):
	"""Plot the effect of the first self-consistent cycle"""
	fgrid, ax = generic_scatterplot(df, x_name, y_name)	
	title = 'Comparison before/after self-consistent first cycle'	
	fgrid.fig.suptitle(title)
	fgrid.fig.subplots_adjust(top=0.9)  # Without this, the title is cut off
	plt.ylabel('Charges before optimization')
	plt.xlabel('Charges from Horton')
	plt.savefig(os.path.join(out_path, 'scatter_before_after.png'))
	print('Successfully plotted: {}'.format(title))
	plt.clf()
	return

def scatterplot_constraints(df, out_path, x_name='q', y_name='q_unconstrained'):
	"""Plot the effect of contraints"""
	fgrid, ax = generic_scatterplot(df, x_name, y_name)
	title = 'Effect of constraints on charges'
	fgrid.fig.suptitle(title)
	fgrid.fig.subplots_adjust(top=0.9)  # Without this, the title is cut off
	plt.ylabel('Unconstrained Charges')
	plt.xlabel('Constrained Charges q')
	plt.savefig(os.path.join(out_path, 'scatter_q_unconstrained.png'))
	print('Successfully plotted: {}'.format(title))
	plt.clf()
	return

def collect_charges():
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
			print('Timestamp is {}'.format(time))
		
			# Use the first charge file to come across as a template	
			df = pd.read_csv(os.path.join(subdir, 'fitted_point_charges.csv'))
			df['timestamp'] = time

			if coll_df is None:
				coll_df = df
			else:
				coll_df = coll_df.append(df)

	# Transform the wide format into a long format version (for easier plotting)
	coll_df = pd.melt(coll_df, id_vars=['atom', 'residue', 'timestamp'], 
			  value_vars=['q', 'q_unconstrained'])
	coll_df = coll_df.rename({'value': 'Charge', 'variable': 'Constraint Type', 'atom': 'Atom'},
				 axis=1)
	return coll_df

def pointplot_errorbars(df, out_path):
	"""Plot stripplot with errorbars"""

	
	fig = plt.figure(figsize=(16,10))

	residue = df.residue.unique()

	# Constrained
	c_df = df.loc[df['Constraint Type'] == 'q']
	pp = sns.pointplot('Charge', 'Atom', data=c_df, scale=1.2, 
			     join=False, err_style='bars', color='black',
			     label='Constrained', zorder=100)
	# unconstrained
	sns.set_palette('pastel')
	uc_df = df.loc[df['Constraint Type'] == 'q_unconstrained']
	sns.pointplot('Charge', 'Atom', hue='residue', data=uc_df, alpha=0.3,
			     join=False, err_style='bars', dodge=0.2, zorder=1, ax=pp.axes, scale = 0.75)


	pp.figure.savefig(os.path.join(out_path, 'pointplot_errorbars.png'))
	return

def main():
	"""TODO """

	# Setup the directory structure
	out_path = create_dir()

	# Collect all individual charges
	collect_df = collect_charges()
	# Pointplot charges with errorbars
	pointplot_errorbars(collect_df, out_path=out_path)


	# Read averaged charges
	input_path = './horton_charges/fitted_point_charges.csv'
	df = pd.read_csv(input_path)

	# Plot constrained vs unconstrained charges
	scatterplot_constraints(df, out_path)

	# Plot old charges vs new ones
	old_charges = extract_init_charges(rtp_path='./md_simulation/n7nh2.rtp', df=df)
	# old_charges = extract_init_charges(rtp_path='modified_charges_sarah_gpaw.rtp', df=df)

	# Put both old and new charges in same table for comparison
	merged_df = old_charges.merge(df, how='inner', on=['atom', 'residue'])
	# Throw away unneccessary information
	merged_df = merged_df[['residue', 'atom', 'q_init', 'q']]
	# Plot the old and new charges 
	scatterplot_comparison(merged_df, out_path)

	# merged_df.to_csv(os.path.join(out_path, 'charge_comparison.csv'))
	print('Done.')


if __name__ == '__main__':
	main()
