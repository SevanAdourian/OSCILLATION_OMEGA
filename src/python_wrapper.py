#! /usr/bin/env python
# run all kernel files

import os
import numpy as np
import cmath
import yaml
import sys
from itertools import islice

import shutil

def load_config(conf_file):
    with open(conf_file) as f:
        conf = yaml.safe_load(f)
        
    return conf

def read_output_fortran(filename):
    with open(filename) as f:
        output = yaml.safe_load(f)

    return output

def cleanup_dir(conf_dict):
    if (conf_dict['kern_rho'] == '.true.' and
        os.path.exists(conf_dict['kernel_rho_file'])):
        os.remove(conf_dict['kernel_rho_file'])
        print("Removing density kernel file")

    if (conf_dict['kern_topo'] == '.true' and
        os.path.exists(conf_dict['kernel_topo_file'])):
        os.remove(conf_dict['kernel_topo_file'])
        print("Removing topography kernel file")

    if (conf_dict['kern_topo'] == '.false.' and
        conf_dict['kern_rho'] == '.false.' and
        os.path.exists(conf_dict['output_fortran'])):
        os.remove(conf_dict['output_fortran'])
        print("Removing normal run file")

    if (os.path.exists(conf_dict['output_oscillation'])):
        os.remove(conf_dict['output_oscillation'])
        print("Removing output file")

def read_1d_model_base(conf_dict):
    # Read the 1d model
    r = np.loadtxt(conf_dict['1d_model_base_file'], skiprows=3, unpack=True)
    conf_dict['1d_model'] = r.transpose()
    conf_dict['nodes'] = np.loadtxt(conf_dict['kern_nodes_file'], dtype=int, unpack=True)

def write_1d_model(conf_dict, ind_kernel_node):

    ind_start = -1
    ind_end = -1
    
    if (conf_dict['kern_rho'] == '.true.'):
        ind_to_duplicate_start = conf_dict['nodes'][ind_kernel_node] - 1
        ind_start = ind_to_duplicate_start
        # ind_to_duplicate_start2 = conf_dict['nodes'][ind_kernel_node]
        ind_to_duplicate_end   = conf_dict['nodes'][ind_kernel_node+1] - 1
        ind_end = ind_to_duplicate_end
        # ind_to_duplicate_end2   = conf_dict['nodes'][ind_kernel_node+1] - 1
        dup_1d = conf_dict['1d_model'].copy()
        to_insert_start = np.array(dup_1d[ind_to_duplicate_start,:])
        to_insert_end   = np.array(dup_1d[ind_to_duplicate_end,:])
        if (conf_dict['1d_model'][ind_to_duplicate_start-1][0] !=
            conf_dict['1d_model'][ind_to_duplicate_start][0]):
            dup_1d = np.insert(dup_1d, ind_start+1, to_insert_start, 0)
            ind_start = ind_start + 1
            ind_end = ind_end + 1

        if (conf_dict['1d_model'][ind_to_duplicate_end-1][0] !=
            conf_dict['1d_model'][ind_to_duplicate_end][0]):
            dup_1d = np.insert(dup_1d, ind_end, to_insert_end, 0)
            ind_end = ind_end + 1
            
        conf_dict['new_1d_model'] = dup_1d

    else:
        dup_1d = conf_dict['1d_model'].copy()
        conf_dict['new_1d_model'] = dup_1d

    if (os.path.exists(conf_dict['1d_model_for_run'])):
        os.remove(conf_dict['1d_model_for_run'])
        print("Removing 1d model")
    
        # Write new to a new file
    with open(conf_dict['1d_model_base_file'], "r") as base_1d:                    
        head = list(islice(base_1d, 3))
        head_replace = list(head)
        
        head_replace = head[0]+head[1]+str(len(dup_1d))+head[2][3:-1]+'\n'
        head = ''.join(head_replace)

    with open(conf_dict['1d_model_for_run'], "w") as new_1d:
        for item in head:
            new_1d.write(item)
            
        np.savetxt(new_1d, dup_1d, fmt="%.2f", delimiter=' ', newline='\n')

    # if (conf_dict['kern_rho'] == '.true.'):
    #     ind_debug = np.reshape(np.arange(1,len(dup_1d)+1),(len(dup_1d),1))
    #     r_debug = np.reshape(dup_1d[:,0],((len(dup_1d),1)))
    #     concat_debug = np.concatenate((ind_debug, r_debug), axis=1)
    #     oned_for_run_debug = '1d_debug'+str(ind_kernel_node)+'.txt'
    #     np.savetxt(oned_for_run_debug,
    #                concat_debug[ind_to_duplicate_start-5:ind_to_duplicate_end+5,:],
    #                delimiter='\t', fmt='%3i %.4f')

    # else:
    #     oned_for_run_debug = '1d_debug_normal.txt'
    #     np.savetxt(oned_for_run_debug,
    #                dup_1d[:,0],  delimiter='\t', fmt='%.4f')
        

    return ind_start, ind_end
        
# def write_3d_rho(conf_dict, ind_kernel_node):
def write_3d_rho(conf_dict, ind_start, ind_end, ii):

    if (conf_dict['kern_rho'] == '.true.'):
        len_1d = len(conf_dict['new_1d_model'])
    else:
        len_1d = len(conf_dict['1d_model'])
        
    rho_3d = np.zeros([len_1d,2], dtype=float)

    print(ind_start, ind_end)
    nb_ind = ind_end - ind_start
    
    delta_rho_single = np.array([conf_dict['rho_perturb'], 0])
    delta_rho_array = np.tile(delta_rho_single, (nb_ind, 1))
    rho_3d[ind_start:ind_end,:] = delta_rho_array
    
    conf_dict['rho_3d'] = rho_3d
    #
    # save_rho_debug = conf_dict['3d_delta_rho'][:-4]+str(ii)+'.dat'
    # print(save_rho_debug)
    np.savetxt(conf_dict['3d_delta_rho'], rho_3d, delimiter='\t')
    # rho_debug = np.concatenate((np.reshape(np.arange(1,len(conf_dict['rho_3d'])+1),
    #                                                  (len(conf_dict['rho_3d']),1)),
    #                             conf_dict['rho_3d']), axis=1)
    # np.savetxt(save_rho_debug, rho_debug[ind_start-5:ind_end+5,:],
    #            delimiter='\t', fmt='%3i\t%.2f\t%.2f')
    
def write_3d_topo(conf_dict, ind_discontinuity):
    len_disc = conf_dict['nb_disc'] + 2
    topo_3d = np.zeros([len_disc,2], dtype=float)
    delta_topo_single = np.array([conf_dict['topo_perturb'], 0])
    topo_3d[ind_discontinuity,:] = delta_topo_single
    #
    # save_topo_debug = conf_dict['3d_delta_topo'][:-4]+str(ind_discontinuity)+'.dat'
    np.savetxt(conf_dict['3d_delta_topo'], topo_3d, delimiter='\t')
    # np.savetxt(save_topo_debug, topo_3d, delimiter='\t')

def run_oscillation_omega(conf_dict):
    cmd = conf_dict['path_to_bin']+' '+conf_dict['1d_model_for_run']+' '+\
        conf_dict['3d_delta_rho']+' '+conf_dict['3d_delta_topo']+' '+\
        conf_dict['kern_rho']+' '+conf_dict['kern_topo']+' >> '+\
        conf_dict['output_oscillation']
    print(cmd)
    os.system(cmd)

def compute_kernel_rho(conf_dict):
    for i in range(0,len(conf_dict['nodes'])-1):
    # for i in range(0,5):
        lay = conf_dict['nodes'][i]
        print(str(lay))
        i1, i2 = write_1d_model(conf_dict, i)
        print(i1, i2)
        write_3d_rho(conf_dict, i1, i2, i)
        write_3d_topo(conf_dict, 0)
        # print('Computing kernel for layer '+str(conf_dict['new_1d_model'][lay][0])+
        #       '/'+str(max(conf_dict['new_1d_model'][:][0])))
        run_oscillation_omega(conf_dict)


def generate_map(conf_dict):
    # delete output file if exists
    cleanup_dir(conf_dict)

    ii = 0
    # find index at which we don't perturb rho anymore
    depth_llsvp = (3480. + conf_dict['height_LLSVP']) * 1000
    i1 = 331 - 1 # Take into account python indexing
    i2 = np.argmin(np.abs(conf_dict['1d_model'][:,0] - depth_llsvp)) + 1

    print(i1, i2)
    # generate array of density perturbations
    density_array = np.linspace(0.01, conf_dict['max_perturbation_rho']/100,
                              num=41, endpoint=True)
    # generate array of topography perturbations at CMB
    topo_array    = np.linspace(0, conf_dict['max_perturbation_topo'],
                                num=21, endpoint=True)
    # topo_array    = np.linspace(0, 0, num=5, endpoint=True)

    output_matrix = np.empty((0,3), float)
    # loop on density perturbations
    for rho_perturb in density_array:
        # loop on topography perturbations
        for topo_perturb in topo_array:
            # overwrite value_perturbation
            print(rho_perturb, topo_perturb)
            i_dum, idum2 = write_1d_model(conf_dict, 0)
            conf_dict['rho_perturb'] = rho_perturb
            conf_dict['topo_perturb'] = topo_perturb
            write_3d_rho(conf_dict, i1, i2, ii)
            write_3d_topo(conf_dict, 1)
            run_oscillation_omega(conf_dict)
            output_f = read_output_fortran(conf_dict['output_fortran'])
            output_matrix = np.append(output_matrix,
                                      np.array([[rho_perturb*100, topo_perturb,
                                                output_f['period']]]), axis=0)
            # print(output_matrix)
            ii += 1
            cleanup_dir(conf_dict)

            
    return output_matrix

def main():
    for config_file in sys.argv[1:]:
        conf = load_config(config_file)
        cleanup_dir(conf)
        read_1d_model_base(conf)
        if (conf['kern_rho'] == '.true.'):
            test2 = compute_kernel_rho(conf)
        elif (conf['kern_topo'] == '.false.' and conf['kern_rho'] == '.false.'):
            file_for_map = 'map_oscillation.dat'
            if (os.path.exists(file_for_map)):
                os.remove(file_for_map)

            output = generate_map(conf)
            np.savetxt(file_for_map, output, fmt='%.2f\t%.2f\t%.2f', delimiter='\t')

    return conf

if __name__ == '__main__':
    test = main()
