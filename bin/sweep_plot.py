"""Plot long format table of point charges.
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

from smamp.tools import find

def collect_average():
	"""Put averaged charges in a dataframe."""
	# Read charges from averaged cost function
	input_path = './horton_charges/fitted_point_charges.csv'
	avg_df = pd.read_csv(input_path)
	# Rename columns for consistency
	avg_df = avg_df.rename({'q': 'averaged cost function'}, axis=1)
	# Transform to long format
	avg_df = pd.melt(avg_df, id_vars=['atom', 'residue'], value_vars=['averaged cost function'])
	avg_df = avg_df.rename({'value': 'charge', 'variable': 'Calculation Variant'}, axis=1)
	return avg_df


def default_style(func):
	"""A decorator for setting global plotting styling options."""
	def wrapper(*args, **kwargs):
		fig = plt.figure(figsize=(16,10))
		sns.set_context("talk", font_scale=0.9)
		plt.xlim(-2, 2)
		plt.tick_params(grid_alpha=0.2)
		func(*args, **kwargs)
		plt.clf()
	return wrapper


@default_style
def pointplot(df, out_path, variant='constrained'):
	"""Pointplot with different colors for every snapshot"""
	c_df = df.loc[df['Calculation Variant'] == variant]
	pp = sns.pointplot('charge', 'atom', data=c_df, scale=1.0, 
			     join=False, hue='timestamp', ci='sd', dodge=0.1)

	pp.set_title('{} charges of all residues'.format(variant).capitalize())
	pp.axes.grid(True)  # Show horizontal gridlines
	pp.figure.savefig(os.path.join(out_path, 'pointplot_{}.png'.format(variant)))

@default_style
def main():
	"""Execute everything. """
	print('This is {}.'.format(__file__))

	print('Collecting ...')

	collect_df = pd.DataFrame()
	for lnrho in range(-10, 0):
		charge_file = find(path='.', folder_keyword='horton_charges/sweep_rhoref', 
				  file_keyword='avg_charges_{}.csv'.format(lnrho), 
				  nr_occ=1)[0]
		df = pd.read_csv(charge_file)
		df['lnrho'] = lnrho
		collect_df = collect_df.append(df)

	collect_df = pd.melt(collect_df, id_vars=['atom', 'residue', 'lnrho'], value_vars=['q'])
	print(collect_df)

	print('Plotting ... ')
	pp = sns.pointplot('value', 'atom', data=collect_df, scale=1.0, 
		 	    join=False, hue='lnrho', ci='sd', dodge=0.1, 
			    palette=sns.color_palette("coolwarm", 10))

	pp.figure.savefig(os.path.join('pointplot.png'))

	print('Done.')


if __name__ == '__main__':
	main()
