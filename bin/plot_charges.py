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
                        # atom names are only unique inside one residuum
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
	# Scatterplot with regression line
	fgrid = sns.lmplot(x=x_name, y=y_name, data=df, fit_reg=True)

	# Add annotations
	ax = fgrid.axes[0,0]
	# Move the labels to the right
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
	plt.ylabel('charges before optimization')
	plt.xlabel('charges from Horton')
	plt.savefig(os.path.join(out_path, 'scatter_before_after.png'))
	print('Successfully plotted: {}'.format(title))
	plt.clf()

def scatterplot_constraints(df, out_path, x_name='q', y_name='q_unconstrained'):
	"""Plot the effect of contraints"""
	fgrid, ax = generic_scatterplot(df, x_name, y_name)
	title = 'Effect of constraints on charges'
	fgrid.fig.suptitle(title)
	fgrid.fig.subplots_adjust(top=0.9)  # Without this, the title is cut off
	plt.ylabel('Unconstrained charges')
	plt.xlabel('Constrained charges q')
	plt.savefig(os.path.join(out_path, 'scatter_q_unconstrained.png'))
	print('Successfully plotted: {}'.format(title))
	plt.clf()

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
	print(coll_df[coll_df.atom.str.contains(r'[1-2]C')])
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


def pointplot_errorbars(df, avg_df, out_path):
	"""Plot stripplot with errorbars"""

	
	fig = plt.figure(figsize=(16,10))

	residue = df.residue.unique()

	# Constrained
	c_df = df.loc[df['Calculation Variant'] == 'q']
	pp = sns.pointplot('charge', 'atom', data=c_df, scale=1.2, 
			     join=False, color='black',
			     label='Constrained', zorder=100,
			     err_style="std")
	# unconstrained
	sns.set_palette('pastel')
	uc = df.loc[df['Calculation Variant'] == 'q_unconstrained']
	sns.pointplot('charge', 'atom', hue='residue', data=uc, alpha=0.3,
		      join=False, err_style='bars', dodge=0.2, zorder=1, ax=pp.axes, 
		      scale = 0.75, type='bar', notch=True)

	# Averaged cost function
	sns.pointplot('q', 'atom', data=avg_df, ax=pp.axes, join=False, color='firebrick', marker=2)

	pp.figure.savefig(os.path.join(out_path, 'pointplot_confidence.png'))
	plt.clf()

def make_edgecolor(ax):
	"""Make boxes transparent with colored edges & whiskers"""
	for i,artist in enumerate(ax.artists):
		# Set the linecolor on the artist to the facecolor, and set the facecolor to None
		col = artist.get_facecolor()
		artist.set_edgecolor(col)
		artist.set_facecolor('None')

		# Each box has 6 associated Line2D objects (to make the whiskers, fliers, etc.)
		# Loop over them here, and use the same colour as above
		for j in range(i*6,i*6+6):
			line = ax.lines[j]
			line.set_color(col)
			line.set_mfc(col)
			line.set_mec(col)

def boxplot(df, out_path):
	""" Boxplot."""
	fig = plt.figure(figsize=(16,10))

	# Individual horton charges
	sns.set_palette(sns.color_palette("Set1", n_colors=4, desat=.9))
	bp = sns.boxplot(x='charge', y='atom', hue='Calculation Variant', data=df, whis=100)
	ax = bp.axes

	split_residues = False	
	if split_residues:
		# Unconstrained
		sns.set_palette('pastel')
		uc = df.loc[df['Calculation Variant'] == 'q_unconstrained']
		bp = sns.boxplot(x='charge', y='atom', hue='residue', data=uc, whis=100)

		# Constrained
		c_df = df.loc[df['Calculation Variant'] == 'q']
		sns.boxplot(x='charge', y='atom', data=c_df, color='black', ax=ax, whis=100)

		# Averaged cost function
		sns.pointplot('q', 'atom', data=avg_df, ax=ax, join=False, color='firebrick')

	make_edgecolor(ax)
	ax.yaxis.grid(True)  # Show horizontal gridlines
	ax.tick_params(grid_alpha=0.5)
	ax.set_xlim(-2, 2)
	bp.figure.savefig(os.path.join(out_path, 'boxplot.png'))
	plt.clf()

def collect_all():
	"""Return Horton, Bader, Averaged charges as one long-format dataframe."""
	# Read charges from averaged cost function
	input_path = './horton_charges/fitted_point_charges.csv'
	avg_df = pd.read_csv(input_path)
	avg_df = avg_df.rename({'q': 'averaged cost function'}, axis=1)
	avg_df = pd.melt(avg_df, id_vars=['atom', 'residue'], value_vars=['averaged cost function'])
	avg_df = avg_df.rename({'value': 'charge', 'variable': 'Calculation Variant'}, axis=1)

	# Collect all horton charges
	print('Collecting HORTON charges ...')
	horton_df = collect_horton()

	# Collect all bader charges
	print('Collecting Bader charges ...')
	bader_df = collect_bader()

	# Paste everything into single dataframe
	print('Combining different charges into one table ... ')
	constr_df = horton_df.loc[horton_df['Calculation Variant'] == 'constrained']
	unconstr_df = horton_df.loc[horton_df['Calculation Variant'] == 'unconstrained']
	collect_df = avg_df
	collect_df = collect_df.append(constr_df, sort=False) 
	collect_df = collect_df.append(unconstr_df, sort=False) 
	collect_df = collect_df.append(bader_df, sort=False)
	return collect_df

def main():
	"""Execute everything. """
	print('This is {}.'.format(__file__))

	# Setup the directory structure
	out_path = create_dir()

	# Collect Horton, Bader, Averaged charges into one dataframe
	collect_df = collect_all()

	print('Plotting ... ')

	# Boxplot charges to get min/max errorbars
	boxplot(collect_df, out_path=out_path)
	print('Boxplot saved.')

	# Pointplot charges with errorbars to get confidence intervals
	# pointplot_errorbars(collect_df, avg_df, out_path=out_path)
	# print('Pointplot saved.')

	# Plot constrained vs unconstrained charges
	# scatterplot_constraints(collect_df, out_path)

	# Plot old charges vs new ones
	# old_charges = extract_init_charges(rtp_path='./md_simulation/n7nh2.rtp', df=avg_df)
	# Put both old and new charges in same table for comparison
	# merged_df = old_charges.merge(avg_df, how='inner', on=['atom', 'residue'])
	# Throw away unneccessary information
	# merged_df = merged_df[['residue', 'atom', 'q_init', 'q']]
	# Plot the old and new charges 
	# scatterplot_comparison(merged_df, out_path)

	# merged_df.to_csv(os.path.join(out_path, 'charge_comparison.csv'))
	print('Done.')


if __name__ == '__main__':
	main()
