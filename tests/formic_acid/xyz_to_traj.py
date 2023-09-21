#!/usr/bin/env python3

#imports
import argparse
from ase import Atoms
from ase.visualize import view
from ase.io import read, write
from ase.io.trajectory import Trajectory
'''
Usage: python xyz_to_traj.py -i file.xyz
'''
#parse command line arguments
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser._action_groups.pop()
required = parser.add_argument_group("required arguments")
required.add_argument("-i", "--input",   type=str,   help="input xyz",  required=True)

from os.path import splitext

args = parser.parse_args()
filename, extension = splitext(args.input)

#read input file 
atoms = read(args.input, index=':')
traj = atoms
new_traj = Trajectory(filename + '.traj', 'w', atoms)
for i in traj:
    new_traj.write(i)

