#!/usr/bin/env python3

import argparse
from ase import Atoms
from ase.io import read, write
from ase.io.trajectory import Trajectory
from os.path import splitext
import io, os
import numpy as np

'''
usage python3 SplitUpMOLPROCOC_ccf12.py -i molecule_geometries.traj -o foldername
'''


#parse command line arguments
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser._action_groups.pop()
required = parser.add_argument_group("required arguments")
required.add_argument("-i", "--input",   type=str,   help="input traj",  required=True)
optional = parser.add_argument_group("optional arguments")
optional.add_argument("-o", "--output",   type=str,   help="output folder name",  default="inp_molpro")
args = parser.parse_args()

#get current working directory and make a scratch 
#directory
path = os.getcwd()
path = path + '/' + args.output
if not os.path.exists(path): os.makedirs(path)
filename, extension = splitext(args.input)


# read input file
traj = Trajectory(args.input)


# define function to translate atomic number to atomic representations
def atomic_numbers_to_labels(Z):
	labels = []
	for z in Z:
		if z == 1:
			labels.append('H')
		elif z == 6: 
			labels.append('C')
		elif z == 8:
			labels.append('O')
		else:
			print("UNSUPPORTED ATOMIC NUMBER", z)
			quit()
	return labels

def calculate_COC(pos, num):
    R = np.zeros(3)
    for i in range(pos.shape[0]):
        R += num[i] * pos[i,:]
    R /= sum(num)
    return R



# get atomic numbers and atomic positions for every structure
for i in range(len(traj)):
    atoms = traj[i]
    string = filename + '_ccf12' + '_%04d' % (i,) +'.inp'
    completeName = os.path.join(path, string)


    
# calculate centre of charge
    pos = atoms.get_positions()
    num = atoms.get_atomic_numbers()
    R = calculate_COC(pos, num)
    pos = pos - R




# write necessary header for Molpro input file
    with open(completeName, "w") as f:
        f.write('geometry={\n')
        f.write(str(len(num)) + '\n')
        f.write('anh freq\n')
        labels = atomic_numbers_to_labels(atoms.get_atomic_numbers())
        for a in range(len(labels)):
            x, y, z = pos[a,:]
            l = labels[a]
            f.write(' ' + l + '    ' + str(x) + '    ' + str(y) + '   ' +  str(z) +"\n")
        f.write('}\n')
        f.write('basis=avtz-f12\n')
        f.write('memory,250,m\n')
        f.write('gthresh,orbital=1.d-8\n')
        f.write('gthresh,energy=1.d-11\n')
        f.write('df-hf\n')
        f.write('force\n')
        f.write('df-ccsd(t)-f12b,ansatz=3C(fix)\n')
        f.write('force\n')
        f.write('expec,type=relax,dm\n')

