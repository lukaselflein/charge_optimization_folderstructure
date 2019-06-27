"""Plot long format table of point charges.
Copyright 2019 Simulation Lab
University of Freiburg
Author: Lukas Elflein <elfleinl@cs.uni-freiburg.de>
"""

import os
import re
import sys
import shutil
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

from smamp.tools import find


def collect_snapshots(plot_range=[-9, -5]):
	collect_df = pd.DataFrame()
	cost_paths = find(path='.', folder_keyword='4_horton_cost_function/lnrho_sweep', 
			  file_keyword='charge',
			  nr_occ=None)

	for charge_file in cost_paths:
		# Parse parameters from filename
		lnrho, sigma = charge_file[-15:-4].split('_')[-2:]

		# Read file
		df = pd.read_csv(charge_file)

		# Paste the lnrho parameter into the dataframe
		df['lnrhoref'] = int(lnrho)
		df['sigma'] = float(sigma)

		# Also note the snapshot identifier
		timestamp = re.findall(r'\d+', charge_file)[0]
		df['snapshot'] = timestamp
		# df['diff'] = (df.q - df.q_unconstrained).abs()
		df['diff'] = (df.q - df.q_unconstrained)**2

		collect_df = collect_df.append(df)

	collect_df = pd.melt(collect_df, id_vars=['atom', 'residue', 'snapshot', 'lnrhoref', 'sigma'], 
			     value_vars=['diff'])
	return collect_df


def collect_avg():
	collect_df = pd.DataFrame()
	cost_paths = find(path='.', folder_keyword='horton_charges/sweep_rhoref',
			  file_keyword='charges',
			  nr_occ=None)

	for charge_file in cost_paths:
		#print(charge_file)
		# Parse parameters from filename
		lnrho, sigma = charge_file[-15:-4].split('_')[-2:]
		#print(lnrho, sigma)

		# Read file
		df = pd.read_csv(charge_file)

		# Paste the lnrho parameter into the dataframe
		df['lnrhoref'] = int(lnrho)
		df['sigma'] = float(sigma)

		# Also note the snapshot identifier
		df['diff'] = (df.q - df.q_unconstrained)**2

		collect_df = collect_df.append(df)

	collect_df = pd.melt(collect_df, id_vars=['atom', 'residue', 'lnrhoref', 'sigma'], 
			     value_vars=['diff'])
	return collect_df

def molten_to_2d(df):
	df = df.groupby(['lnrhoref', 'sigma']).mean().apply(np.sqrt)
	df = df.unstack('sigma')
	df.columns = df.columns.droplevel()
	df = df.transpose()
	df = df.sort_index(ascending=False)
	return df

def plot_heatmap(df, outname):
	fig = plt.figure(figsize=(16,10))
	sns.set_context("talk", font_scale=0.9)
	pp = sns.heatmap(data=df, cmap='coolwarm') 
	pp.set_title('RRMSD between constrained and unconstrained charges [e].')
	pp.figure.savefig(outname)

def main():
	"""Execute everything. """
	print('This is {}.'.format(__file__))


	print('Collecting averages')
	df = collect_avg()
	df = molten_to_2d(df)
	plot_heatmap(df, outname='avg_heatmap.png')
	print(df)	

	# Individual snapshots
	print('Collecting snapshots ...')
	df = collect_snapshots(plot_range=range(-9, 0))
	df = molten_to_2d(df)
	plot_heatmap(df, outname='snapshot_heatmap.png')
	print(df)	

	print('Done.')


if __name__ == '__main__':
	main()
