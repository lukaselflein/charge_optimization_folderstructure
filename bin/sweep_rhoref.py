"""Vary the rhoref parameter to find a sane value.
Copyright 2019 Simulation Lab
University of Freiburg
Author: Lukas Elflein <elfleinl@cs.uni-freiburg.de>
"""

from __future__ import print_function

import os
import shutil
#import multiprocessing
#import sys
import random

from loop_cost_functions import calc_cost_function 
from smamp.tools import cd
from smamp.tools import check_existence
def testprint(*args, **kwargs):
   return 'args: {}, kwargs: {}'.format(args, kwargs)

def get_tasks(path_to_subdir):
   """Vary the lnrho weighting parameter, create folder and execute."""
   sweep_dir = 'lnrho_sweep'
   if os.path.exists(sweep_dir):
      # print('Removing old dir')
      # shutil.rmtree(sweep_dir)
      pass
   else:
      print('making dir')
      os.mkdir(sweep_dir)
      print('dir made.')

   tasks = []
   skipped = 0
   for sigma in [0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4]:
      #for lnrho in [-9, -8, -7, -6, -5, -4, -3, -2]:
      for lnrho in [-7, -6, -5.5, -5, -4.75, -4.5, -4.25, -4, -3.5]:
         output_name = os.path.join(sweep_dir, 'cost_{}_{}.h5'.format(lnrho, sigma))
         if os.path.exists(output_name):
            # print('{} exists. Do not include in worklist.'.format(output_name))
            skipped += 1
            continue
         
         else:
            tasks += [(path_to_subdir, lnrho, sigma, output_name)]
   print('{} files found and skipped.'.format(skipped))
   return tasks

def calculate_tasks(tasks):
   print('{} items in worklist.'.format(len(tasks)))
   random.shuffle(tasks)
   for task in tasks:
      if os.path.exists(task[-1]):
         #print('{} exists. Skipping ahead.'.format(task[-1]))
         continue

      print('starting {} {} in {}'.format(task[1], task[2], task[0]))
      with cd(task[0]):
         calc_cost_function(*task)
      

def main():
   """ Execute everything."""
   print('This is {}.'.format(__file__))
   print('Current working dir: {}'.format(os.getcwd()))
   
   tasks = []

   # Crawl the directory structure
   for subdir, dirs, files in sorted(os.walk('.')):

      # Exclude template folders from search
      if 'template' in subdir or 'exclude' in subdir or 'lnrho_sweep' in subdir:
         continue

      # Select the folder to calculate in
      if 'horton_cost_function' in subdir:
         print('Moving to {}'.format(subdir))
         with cd(subdir):
            subdir_tasks = get_tasks(subdir)
            tasks += subdir_tasks
   
   calculate_tasks(tasks)
   print('Done.')

if __name__ == '__main__':
   main()
