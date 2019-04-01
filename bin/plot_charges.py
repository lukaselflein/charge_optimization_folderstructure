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


@default_style
def boxplot(df, out_path):
	""" Plot charges calculated under different conditions in one boxplot."""
	# sns.set_palette(sns.color_palette("Set1", n_colors=4, desat=.9))
	bp = sns.boxplot(x='charge', y='atom', hue='Calculation Variant', data=df, whis=100)

	ax = bp.axes
	ax.grid(True)  # Show horizontal gridlines
	make_edgecolor(ax)

	bp.set_title('Distribution of charges with different calculation methods')

	bp.figure.savefig(os.path.join(out_path, 'boxplot.png'))


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
def swarmplot(df, out_path, variant='constrained'):
	"""Plot swarmplot with different colors for every snapshot"""
	c_df = df.loc[df['Calculation Variant'] == variant]
	sp = sns.swarmplot('charge', 'atom', data=c_df, hue='timestamp')

	sp.set_title('{} charges of all residues'.format(variant).capitalize())
	sp.axes.grid(True)  # Show horizontal gridlines
	sp.figure.savefig(os.path.join(out_path, 'swarmplot_{}.png'.format(variant)))

@default_style
def average_pointplot(df, out_path):
	"""Pointplot of averages to compare different calculation methods"""
	pp = sns.pointplot('charge', 'atom', data=df, scale=1.0, 
	                   join=False, hue='Calculation Variant', ci=None, dodge=0.2)

	pp.set_title('Averages of charges with different calculation methods')
	pp.axes.grid(True)  # Show horizontal gridlines
	pp.figure.savefig(os.path.join(out_path, 'averages.png'))

def main():
	"""Execute everything. """
	print('This is {}.'.format(__file__))

	# Setup the directory path
	out_path = './plotting'

	# Load Horton, Bader, Averaged charges from table
	print('Reading dataset ...')
	collect_df = pd.read_csv('./plotting/all_charges.csv')

	print('Plotting ... ')

	# Pointplot to compare mean values
	average_pointplot(collect_df, out_path=out_path)
	print('Average plot saved.')

	# Swarmplot charges to resolve timestamp differences
	for variant in ('constrained', 'unconstrained', 'bader'):
		swarmplot(collect_df, out_path=out_path, variant=variant)
		print('Swarmplot of {} saved.'.format(variant))

	# Boxplot charges to get min/max errorbars
	boxplot(collect_df, out_path=out_path)
	print('Boxplot saved.')

	print('Done.')


if __name__ == '__main__':
	main()
