#! /usr/bin/env python
# run all kernel files

import os
import numpy as np
import cmath
import yaml
import sys
from itertools import islice

def load_config(conf_file):
    with open(conf_file) as f:
        conf = yaml.safe_load(f)
        
    return conf

def cleanup_dir(conf_dict):
    if (conf_dict['kern_rho'] == '.true.' and os.path.exists(conf_dict['kernel_rho_file'])):
        os.remove(conf_dict['kernel_rho_file'])
        print("Removing density kernel file")

    if (conf_dict['kern_topo'] == '.true' and os.path.exists(conf_dict['kernel_topo_file'])):
        os.remove(conf_dict['kernel_topo_file'])
        print("Removing topography kernel file")

    if os.path.exists(conf_dict['output_oscillation']):
        os.remove(conf_dict['output_oscillation'])
        print("Removing output file")

def read_1d_model_base(conf_dict):
    # Read the 1d model
    r = np.loadtxt(conf_dict['1d_model_base_file'], skiprows=3, unpack=True)
    conf_dict['1d_model'] = r.transpose()
    conf_dict['nodes'] = np.loadtxt(conf_dict['kern_nodes_file'], dtype=int, unpack=True)

def write_1d_model(conf_dict, ind_kernel_node):
    ind_to_duplicate = conf_dict['nodes'][ind_kernel_node]
    dup_1d = conf_dict['1d_model'].copy()
    to_insert = np.array(dup_1d[ind_to_duplicate,:])
    if (conf_dict['1d_model'][ind_to_duplicate][0] !=
        conf_dict['1d_model'][ind_to_duplicate+1][0]):
        dup_1d = np.insert(dup_1d, ind_to_duplicate, to_insert, 0)
    conf_dict['new_1d_model'] = dup_1d

    # Write new to a new file
    with open(conf_dict['1d_model_base_file'], "r") as base_1d:
        head = list(islice(base_1d, 3))
        head_replace = list(head)
        head_replace = head[0]+head[1]+str(int(head_replace[2][0:3])+1)+head[2][3:-1]+'\n'
        head = ''.join(head_replace)

    # always remember, use files in a with statement
    with open(conf_dict['1d_model_for_run'], "w") as new_1d:
        for item in head:
            new_1d.write(item)
            
        np.savetxt(new_1d, dup_1d, fmt="%.2f", delimiter=' ', newline='\n')

def write_3d_rho(conf_dict, ind_kernel_node):
    len_new_1d = len(conf_dict['new_1d_model'])
    rho_3d = np.zeros([len_new_1d,2], dtype=float)
    ind_start = conf_dict['nodes'][ind_kernel_node] - 1
    ind_end = conf_dict['nodes'][ind_kernel_node+1] - 1
    nb_ind = ind_end - ind_start
    delta_rho_single = np.array([conf_dict['rho_perturb'], 0])
    delta_rho_array = np.tile(delta_rho_single, (nb_ind, 1))
    rho_3d[ind_start:ind_end,:] = delta_rho_array
    conf_dict['rho_3d'] = rho_3d
    #
    np.savetxt(conf_dict['3d_delta_rho'], rho_3d, delimiter='\t')
    
def write_3d_topo(conf_dict, ind_discontinuity):
    len_disc = conf_dict['nb_disc'] + 1
    topo_3d = np.zeros([len_disc,2], dtype=float)
    delta_topo_single = np.array([conf_dict['topo_perturb'], 0])
    topo_3d[ind_discontinuity,:] = delta_topo_single
    # 
    np.savetxt(conf_dict['3d_delta_topo'], topo_3d, delimiter='\t')

def run_oscillation_omega(conf_dict):
    cmd = conf_dict['path_to_bin']+' '+conf_dict['1d_model_for_run']+' '+\
        conf_dict['3d_delta_rho']+' '+conf_dict['3d_delta_topo']+' '+\
        conf_dict['kern_rho']+' '+conf_dict['kern_topo']+' >> '+\
        conf_dict['output_oscillation']
    print(cmd)
        # os.system(cmd)

    
def main():
    for config_file in sys.argv[1:]:
        conf = load_config(config_file)
        # cleanup_dir(conf)
        read_1d_model_base(conf)
        ind_to_duplicate = int(2)
        write_1d_model(conf, ind_to_duplicate)
        write_3d_rho(conf, ind_to_duplicate)
        write_3d_topo(conf, ind_to_duplicate)
        run_oscillation_omega(conf)
    return conf

if __name__ == '__main__':
    test = main()
