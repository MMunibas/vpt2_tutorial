#!/usr/bin/env python3

# imports
import argparse
from ase import Atoms

from ase.io import read, write
from ase.io.trajectory import Trajectory
from os.path import splitext
import io, os
'''
Usage: python Traj2Xyz.py -i file.xyz
'''
#parse command line arguments
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser._action_groups.pop()
required = parser.add_argument_group("required arguments")
required.add_argument("-i", "--input",   type=str,   help="input xyz",  required=True)

args = parser.parse_args()
filename, extension = splitext(args.input)
traj = Trajectory(args.input)


#output file name
outFileName = filename + '.xyz'

write(outFileName, traj)

