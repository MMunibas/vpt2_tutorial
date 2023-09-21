#!/usr/bin/env python3

# imports
import argparse
from ase import Atoms
from ase.visualize import view
from ase.io import read, write
from ase.calculators.mopac import MOPAC
from ase.optimize import *
from ase import units
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.langevin import Langevin
from ase.io.trajectory import Trajectory
from os.path import splitext
import io, os

#parse command line arguments
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser._action_groups.pop()
required = parser.add_argument_group("required arguments")
required.add_argument("-i", "--input",   type=str,   help="input xyz",  required=True)
optional = parser.add_argument_group("optional arguments")
optional.add_argument("--label",   type=str,   help="prefix of calculator files",  default="calc/mopac-opt")
optional.add_argument("--charge",  type=float, help="total charge", default=0.0)
optional.add_argument("--temperature", type=float, help="Set momenta corresponding to a temperature T", default=300)
optional.add_argument("--timestep",  type=float, help="timestep for Langevin algorithm", default=1)
optional.add_argument("--friction",  type=float, help="friction coeff for Langevin algorithm", default=0.02)
optional.add_argument("--interval",  type=float, help="interval", default=50)

args = parser.parse_args()
filename, extension = splitext(args.input)

traj = Trajectory(args.input)


#output file name
outFileName = filename + '.xyz'

write(outFileName, traj)

