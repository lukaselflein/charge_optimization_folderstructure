""" Extract charges obtained via HORTON and Bader.
Copyright 2019 Simulation Lab
University of Freiburg
Author: Lukas Elflein <elfleinl@cs.uni-freiburg.de>
"""

import os
import re
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from smamp.tools import find
from io import StringIO

def create_dir(path='./plotting'):
	"""Create new folder for pictures if it does not exist yet."""
	if os.path.isdir(path):
		return path

	os.makedirs(path)
	return path


def collect_bfgs_energies():
	"""Find BFGS energies and put them in one dataframe."""
	# Initialize collection data frame
	coll_df = None

	# Crawl the directory structure
	bfgs_list = find(path='./', folder_keyword='dft', file_keyword='BFGS', nr_occ=None)

	for energy_file in bfgs_list:
		print('Moving to {}'.format(energy_file))
		topdir = os.path.split(os.path.split(energy_file)[0])[0]
		time = topdir.replace('./', '').replace('_ps_snapshot', '')
		time = int(time)

		with open(energy_file) as data:
			string = data.read()
			string = re.sub(r"\[.+\]", "", string)
			string = re.sub(r"\*", "", string)
			string = re.sub('BFGSLineSearch:  ', '', string)

		sio = StringIO(string)
		df = pd.read_csv(sio, sep='\s+', engine='python', skiprows=[1])

		df['Timestamp'] = time

		if coll_df is None:
			coll_df = df
		else:
			coll_df = coll_df.append(df)

	return coll_df


def energy_plot(df, out_path, log=False):
	"""Lineplot of energy convergence."""
	if log:
		fig, ax = plt.subplots()
		ax.set_yscale('symlog')
		lp = sns.lineplot('Step', 'Energy', data=df, hue='Timestamp', legend='full', ax=ax)

	else:
		lp = sns.lineplot('Step', 'Energy', data=df, hue='Timestamp', legend='full')

	lp.set_title('Energy convergence of the different snapshots')
	lp.figure.savefig(os.path.join(out_path, 'energy_convergence.png'))
	plt.cla()


def force_plot(df, out_path):
	"""Lineplot of Force convergence."""
	lp = sns.lineplot('Step', 'fmax', data=df, hue='Timestamp', markers=True, legend='full')

	lp.set_title('Force convergence of the different snapshots')
	lp.figure.savefig(os.path.join(out_path, 'force_convergence.png'))
	plt.cla()


def main():
	"""Collect energies and save them to .csv file"""

	create_dir(path='./plotting')

	print('Collecting BFGS files ...')
	collect_df = collect_bfgs_energies()

	if collect_df is None:
		print('No BFGS file was found. Done.')
		exit()

	# Save long format dataset
	print('Saving dataframe ...')
	collect_df.to_csv('./plotting/energy_convergence.csv')

	print('Plotting ...')
	energy_plot(collect_df, './plotting')
	force_plot(collect_df, './plotting')
	print('Done.')

if __name__ == '__main__':
	main()
