"""Plot long format table of point charges.
Copyright 2019 Simulation Lab
University of Freiburg
Author: Lukas Elflein <elfleinl@cs.uni-freiburg.de>
"""

import os
import re
import sys
import shutil
import itertools
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from smamp.tools import find

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


def collect_averages():
   collect_df = pd.DataFrame()
   
   cost_paths = find(path='.', folder_keyword='horton_charges/sweep_rhoref', 
           file_keyword='charge',
           nr_occ=None)

   for charge_file in cost_paths:
      # Parse parameters from filename
      lnrho, sigma = charge_file[-15:-4].split('_')[-2:]
      if sigma == '0.8':
         charge_file = find(path='.', folder_keyword='horton_charges/sweep_rhoref', 
                 file_keyword='charges_{}_{}.csv'.format(lnrho, sigma), 
                 nr_occ=1)[0]
         df = pd.read_csv(charge_file)
         df['lnrho'] = int(lnrho)
         collect_df = collect_df.append(df)

   collect_df = pd.melt(collect_df, id_vars=['atom', 'residue', 'lnrho'], value_vars=['q'])
   return collect_df


@default_style
def plot_averages(df):
   pp = sns.pointplot('value', 'atom', data=df, scale=1.0, 
              join=False, hue='lnrho', ci='sd', dodge=0.1, 
             palette=sns.color_palette("coolwarm", 10))

   pp.set_title('Averaged-Costfunction charges of all residues')
   pp.axes.grid(True)  # Show horizontal gridlines
   pp.figure.savefig('plotting/pointplot.png')



def collect_snapshots():
   cost_paths = find(path='.', folder_keyword='4_horton_cost_function/lnrho_sweep', 
           file_keyword='charges', 
           nr_occ=None)
   collect_df = pd.DataFrame()
   for charge_file in cost_paths:
      # Parse parameters from filename
      lnrho, sigma = charge_file[-15:-4].split('_')[-2:]
      if sigma == '0.8':
         df = pd.read_csv(charge_file)

         # Paste the lnrho parameter into the dataframe
         df['lnrho'] = lnrho

         # Also note the snapshot identifier
         timestamp = re.findall(r'\d+', charge_file)[0]
         df['snapshot'] = timestamp

         collect_df = collect_df.append(df)

   collect_df = pd.melt(collect_df, id_vars=['atom', 'residue', 'lnrho', 'snapshot'], 
              value_vars=['q'])
   return collect_df


@default_style
def plot_snapshots(df):
   pp = sns.pointplot('value', 'atom', data=df, scale=1.0, 
              join=False, hue='lnrho', ci='sd', dodge=0.1, 
             palette=sns.color_palette("coolwarm", 10))

   pp.set_title('Individual snapshot charges of all residues')
   pp.axes.grid(True)  # Show horizontal gridlines
   pp.figure.savefig('plotting/snapshots_pointplot.png')


@default_style
def plot_joint(avg_df, snap_df):
   pp = sns.pointplot('value', 'atom', data=snap_df, scale=0.8, 
              join=False, hue='lnrho', ci='sd', dodge=0.1, 
             palette=sns.color_palette("coolwarm", 10))
   pp = sns.pointplot('value', 'atom', data=avg_df, scale=1.0, 
              join=False, hue='lnrho', ci='sd', dodge=0.1, 
             palette=sns.color_palette("YlGn", 10),
             markers='+')

   pp.set_title('Charges of all residues')
   pp.axes.grid(True)  # Show horizontal gridlines
   pp.figure.savefig('plotting/joint_pointplot.png')

@default_style
def swarmplot(df):
   sp = sns.swarmplot('value', 'atom', data=df, hue='lnrho',
             palette=sns.color_palette("coolwarm", 10))

   sp.set_title('Individual snapshot charges of all residues')
   sp.axes.grid(True)  # Show horizontal gridlines
   sp.figure.savefig('plotting/swarmplot.png')

def plot_variance(df):
   dispersions = []
   lnrhos = []
   for lnrho in df.lnrho.unique():
      lnrhos += [int(lnrho)]
      dispersion = 0
      rho_df = df.loc[df['lnrho'] == lnrho]
      for residue in df.residue.unique():
         res_df = rho_df.loc[rho_df['residue'] == residue]
         for atom in res_df.atom.unique():
            atom_df = res_df.loc[res_df['atom'] == atom]
            dispersion += abs(atom_df.value.max() - atom_df.value.min())
      dispersions += [dispersion]

   x, y = zip(*sorted(zip(lnrhos, dispersions)))

   fig = plt.figure(figsize=(16,10))
   sns.set_context("talk", font_scale=0.9)
   plt.plot(x, y, marker='o')
   plt.xlabel('rhoref')
   plt.ylabel('cumlatative dispersion of charge [e]')
   plt.savefig('plotting/dispersion.png')
   plt.clf()

def main():
   """Execute everything. """
   print('This is {}.'.format(__file__))

   # Averages
   print('Collecting averages ...')
   average_df = collect_averages()
   print('Plotting averages ...')
   plot_averages(average_df)

   # Individual snapshots
   print('Collecting snapshots ...')
   snapshot_df = collect_snapshots()
   plot_variance(snapshot_df)
   print('Plotting snapshots ...')
   plot_snapshots(snapshot_df)
   swarmplot(snapshot_df)
   # plot_joint(average_df, snapshot_df)

   print('Done.')


if __name__ == '__main__':
   main()
