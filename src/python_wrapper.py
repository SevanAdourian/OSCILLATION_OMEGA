#! /usr/bin/env python
# run all kernel files

import os
import numpy as np
import cmath
import yaml
import sys

def load_config(conf_file):
    with open(conf_file) as f:
        conf = yaml.load(f)
    return conf
    
# def read_1d_model(conf):
    # Read the 1d model
    # r = np.loadtxt(conf.1d_model_file,skiprows=3,usecols=(0),unpack=True)
    # conf.r = r/1.e3

def cleanup_dir(conf):
    if (kern_rho == '.true.' and os.path.exists(kernel_rho_file)):
        os.remove(kernel_rho_file)
        print("Removing kernel file")

    if (kern_topo=='.true' and os.path.exists(kernel_topo_file)):
        os.remove(kernel_topo_file)
        print("Removing kernel file")

    if os.path.exists(output_oscillation):
        os.remove(output_oscillation)
        print("Removing output file")

    return

def main():
    for conf_file in sys.argv[1:]:
        conf = load_config(conf_file)
        # cleanup_dir(conf)
        # read_1d
        # write_1d_model(conf)
        # write_3d_rho(conf)
        # write_3d_topo(conf)
        # run_oscillation_omega(conf)

if __name__ == '__main__':
    main()

        
# if kern_rho == '.true.':
#     for ii in range(2,nb_disc-1):
#         cmd = path_to_bin+' '+str(ii)+' '+str(value_perturb)+' '+\
#             kern_rho+' '+kern_topo+' >> '+output_oscillation
#         print(cmd)
#         os.system(cmd)

# if kern_topo == '.true.':
#     for disc_perturb in range(2, 13):
#         cmd = path_to_bin+' '+str(disc_perturb)+' '+str(value_perturb)+' '+\
#             kern_rho+' '+kern_topo+' >> '+output_oscillation
#         print(cmd)
#         os.system(cmd)

