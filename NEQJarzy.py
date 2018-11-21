#!/usr/bin/env python3

# coding: utf-8
# # $\Gamma_{NEQ}(x)$ calculation routine
# * target: calculate $\Gamma_{NEQ}(t)=\frac{1}{k_BT}\int_0^tdt'\left<\delta f_c(t)\delta f_c(t')\right>$
# ## Data import

import sys
import numpy as np
import scipy, scipy.integrate as inte
from scipy.ndimage.filters import gaussian_filter
import argparse
import glob
import pandas as pd


def main():

  parser = argparse.ArgumentParser(description='Integrates a constraint force file via trapezoid rule, and performs a friction correction based on Jarzynskis fast growth estimator. First column: reaction coordinate in nm calculated via t * vel. Second column: mean work profile. Third column: friction factors. Fourth column: Jarzynski estimator free energy profile. ATTENTION: Use with python3 or higher!')
  parser.add_argument('-i', metavar='<xvg force file>', type=str, help='xvg constraint force files prefix as given by Gromacs mdrun -pf option before running number. I.e. *[000-XX]')
  parser.add_argument('-s', metavar='<xvg force file>', type=str, help='xvg constraint force files suffix as given by Gromacs mdrun -pf option after running number. I.e. [000-XXX]*.xvg')
  parser.add_argument('-o', metavar='<combined results>', type=str, help='file to write x, dG(x), friction coefficeint by integration (time resolved), and the friction-corrected dG(x).')
  parser.add_argument('-ofrict', metavar='<combined friction results>', type=str, help='file to write x, ACF, friction coefficeint by integration (time resolved), gauss filtered friction coefficient, and slide window averaged friction.')
  parser.add_argument('-vel', metavar='<pull velocity>', type=float, help='pull velocity in nm/ps for converting simulation time t into distance x')
  parser.add_argument('-T', metavar='<temperature>', type=float, help='temperature in K')
  parser.add_argument('-N', metavar='<Number of trajectories>', type=int, help='number of trajectories to average F(x)')
  parser.add_argument('-sigma', metavar='<gauss blurr>', type=int, help='sigma value for Gauss filter for displaying Gamma(x) (recommended: 4 per 100 data points)')

  args = parser.parse_args()

  start_name = args.i
  suffix = args.s
  N = args.N
  vel = args.vel
  RT = 0.0083144598*args.T

  sys.stdout.write("reading data...\n")
  file_names = glob.glob("{}*{}.xvg".format(start_name,suffix))

# ## read in initial data to get length of necessary array

  test_file = pd.read_csv(file_names[0],sep='\s+',header=None,skiprows=17,dtype=float)
  length_data = len(test_file[0].values)
  full_force_set = np.zeros((N,length_data))
  x = np.zeros(length_data)
  t = np.zeros(length_data)
  t = test_file[0].values
  x = test_file[0].values * vel

# ## read in data

  for i in range(0,N):
      current_file_name = file_names[i]
      sys.stdout.write("reading file {}\n".format(current_file_name))
      input_file_data = pd.read_csv(current_file_name,sep='\s+',header=None,skiprows=17,dtype=float)
      full_force_set[i,:] = input_file_data[1].values

# ## preprocessing
# * integrate force to obtain W and <W>
  full_work_set = np.zeros((N,length_data))

  for i in range(0,N):
      full_work_set[i,1:] = scipy.integrate.cumtrapz(full_force_set[i,:],x)

  av_work = np.zeros(length_data)
  for i in range(length_data):
      av_work[i] = np.mean(full_work_set[:,i])

# ## evaluation
# * calculate Jarzynskis estimator: 
  sys.stdout.write("Calculating Jarzynskis fast growth estimator...\n")

  G_Jarz = np.zeros(length_data)
  Wdiss_Jarz = np.zeros(length_data)
  for i in range(length_data):
      G_Jarz[i] = (-1)*RT*np.log( np.mean( np.exp((-1)*full_work_set[:,i]/RT) ) )

  Wdiss_Jarz = av_work - G_Jarz

# * Gauss filter:
  sys.stdout.write("applying Gauss filter...\n") 
  blurr = args.sigma
  blurred_Jarz = np.zeros(length_data)
  blurred_Jarz = gaussian_filter(Wdiss_Jarz, sigma=blurr)

# * differentiate smoothened result to recalculate Gamma values in the time domain
  Gamma = np.zeros(length_data)
  Gamma[1:] = np.diff(blurred_Jarz)/(vel**2)
  

  sys.stdout.write("writing output...\n")
  distresult = open(args.o,"w")
  frictresult = open(args.ofrict,"w")

  distresult.write("#x  <W>  Gamma   wdiss   G_Jarz\n")
  for i in range(length_data):
     distresult.write("{:15.8f} {:20.8f} {:20.8f} {:20.8f} {:20.8f}\n".format(x[i],av_work[i],Gamma[i],Wdiss_Jarz[i],G_Jarz[i]))

  frictresult.write("#x   Gamma   gauss_filtered_Wdiss  \n")
  for i in range(length_data):
     frictresult.write("{:15.8f} {:20.8f} {:20.8f}\n".format(x[i],Gamma[i],blurred_Jarz[i]))

  distresult.close()
  frictresult.close()

  sys.stdout.write("Done!\n")

if __name__ == "__main__":
        main()
