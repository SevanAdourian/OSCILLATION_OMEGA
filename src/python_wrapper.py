#! /usr/bin/env python
# run all kernel files

import os
import numpy as np
import cmath

path_to_bin = '/Users/sevan/Documents/PhD/Research/OMEGA_PHD/scripts/FORTRAN/OMEGA_OSCILLATION_NEW/build/oscillation_omega'
path_to_1d  = '../data/make_s20rts/isoprem808.md'
value_perturb = 0.01
dist_between_knots = 80. # in km
R_EARTH = 6371.
R_CMB = 3480.
pts = np.arange(R_CMB, R_EARTH, dist_between_knots)
# npts = int(np.ceil((R_EARTH-R_CMB)/dist_between_knots))
kernel_rho_file='../build/kernel_rho.txt'
kernel_topo_file='../build/kernel_topo.txt'
output_oscillation='../build/out1'
# print(npts)

# Add flag for kernel computation or not
kern_rho = '.true.'
kern_topo = '.false.'
# Read the 1d model
r = np.loadtxt(path_to_1d,skiprows=3,usecols=(0),unpack=True)
r = r/1.e3

if os.path.exists(kernel_rho_file):
    os.remove(kernel_rho_file)
    print("Removing kernel file")

if os.path.exists(kernel_topo_file):
    os.remove(kernel_topo_file)
    print("Removing kernel file")

if os.path.exists(output_oscillation):
    os.remove(output_oscillation)
    print("Removing output file")

for pt in pts:
# for idx in range(1, len(r)):
    # layer_depth = nl*dist_between_knots
    idx = np.max(np.where((np.abs(r - pt)) == np.abs(r - pt).min()))
    # idx2 = np.max(np.where((np.abs(r - pt+dist_between_knots)) == np.abs(r - pt).min()))
    # print(idx,r[idx])
    print(pt, idx, r[idx])
    cmd = path_to_bin+' '+str(idx+1)+' '+str(value_perturb)+' '+kern_rho+' '+kern_topo+' >> '+output_oscillation
    print(cmd)
    os.system(cmd)

# for disc_perturb in range(2, 13):
# # for idx in range(1, len(r)):
#     cmd = path_to_bin+' '+str(disc_perturb)+' '+str(value_perturb)+' '+kern_rho+' '+kern_topo+' >> '+output_oscillation
#     print(cmd)
#     os.system(cmd)

