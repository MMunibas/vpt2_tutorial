#!/usr/bin/env python3

#imports
import argparse
from ase import Atoms
import numpy as np
from ase.visualize import view
from ase.io import read, write
from ase.optimize import *
from ase import units
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.langevin import Langevin
from ase.io.trajectory import Trajectory

#parse command line arguments
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser._action_groups.pop()
required = parser.add_argument_group("required arguments")
required.add_argument("-i", "--input",   type=str,   help="input xyz",  required=True)
optional = parser.add_argument_group("optional arguments")
optional.add_argument("--charge",  type=float, help="total charge", default=0.0)
optional.add_argument("--fmax",  type=float, help="maximal force", default=0.0001)
from os.path import splitext

args = parser.parse_args()
filename, extension = splitext(args.input)

#read input file 
atoms = read(args.input, index=':')
traj = atoms
new_traj = Trajectory(filename + '.traj', 'w', atoms)
for i in traj:
    new_traj.write(i)

